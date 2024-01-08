//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef HRSCC_CHARACTERISTIC_SPLIT_HH
#define HRSCC_CHARACTERISTIC_SPLIT_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#include <cctk.h>

#include <utils.hh>

namespace hrscc {

//! performs a characterist-wise split of a system of equations
/*!
 *  This will reconstruct the fluxes by decomposing the system of equations
 *  into its pseudo-characteristic variables and using the given flux-splitter
 *  to reconstruct the fluxes of each one. At the end the fluxes are
 *  transformed back to the original variables.
 *
 *  \tparam claw_t conservation law that we want to split
 *  \tparam flux_splitter_t flux-splitter to use for the reconstruction
 */
template<typename claw_t, typename flux_splitter_t>
class CharacteristicSplit {
    public:
        typedef claw_t claw;
        typedef flux_splitter_t flux_splitter;
        typedef typename flux_splitter::reconstructor reconstructor;

        //! width of the stencil used for the reconstruction
        enum {width = flux_splitter::width};
        //! number of points needed for the reconstruction
        enum {size = 2*width};

        //! if true the method requires the eigenvalues at the interface
        static bool const requires_interface_eigenvalues =
               flux_splitter::requires_interface_speed;
        //! if true the method requires all the eigenvalues
        static bool const requires_all_eigenvalues =
               flux_splitter::requires_speed
            || flux_splitter::requires_max_speed
            || flux_splitter::requires_min_speed;
        //! if true the method requires the eigenvectors at the interface
        static bool const requires_interface_eigenvectors = true;
        //! if true the method requires all the eigenvectors
        static bool const requires_all_eigenvectors = false;

        CharacteristicSplit() {
            for(int i = 0; i < claw::nequations; ++i) {
                for(int j = 0; j < size; ++j) {
                    _M_eta[i][j] = 0;
                    _M_phi[i][j] = 0;
                }
            }
        }

        //! compute the component of the flux vector at the given point
        /*!
         *  The input arrays should be centered around \e i so that, at the end,
         *  we reconstruct in 1/2.
         *
         *  Superfluous arrays can be put to NULL.
         */
        void compute_fluxes(
                //! [in] grid spacing
                CCTK_REAL delta,
                //! [in] total size of the grid, used to deduce the
                //! displacement of the data in the memory.
                int lsh,
                //! [in] The maximum propagation speed of a signal
                CCTK_REAL cbound,
                //! [in] conserved variables around the reconstruction point,
                //! it should be stored in memory as a claw::nequations x lsh
                //! array
                CCTK_REAL const * u,
                //! [in] fluxes around the reconstruction point, it should be
                //! stored in memory as a claw::nequations x lsh array
                CCTK_REAL const * f,
                //! [in] eigenvalues at the reconstruction point
                CCTK_REAL const * interface_eigenvalue,
                //! [in] eigenvalues on the stencil, it should be stored in
                //! memory as a claw::nequations x lsh array
                CCTK_REAL const * eigenvalue,
                //! [in] left eigenvectors at the reconstruction point, stored
                //! in row-major order
                CCTK_REAL const * interface_left_eigenvector,
                //! [in] right eigenvectors at the reconstruction point, stored
                //! in row-major order
                CCTK_REAL const * interface_right_eigenvector,
                //! [in] left eigenvectors on the stencil, they should be stored
                //! in memory as a claw::nequations x claw::nequations x lsh
                //! array
                CCTK_REAL const *, // left_eigenvector,
                //! [in] right eigenvectors on the stencil, they should be
                //! stored in memory as a claw::nequations x claw::nequations
                //! x lsh array
                CCTK_REAL const *, // right_eigenvector,
                //! [out] reconstructed fluxes
                CCTK_REAL * flux
                ) const {
            CCTK_REAL max_speed = - std::numeric_limits<CCTK_REAL>::max();
            CCTK_REAL min_speed =   std::numeric_limits<CCTK_REAL>::max();

            if(flux_splitter::requires_max_speed) {
                for(int v = 0; v < claw::nequations; ++v) {
                    for(int i = - width + 1; i < width + 1; ++i) {
                        max_speed = std::max(max_speed, eigenvalue[v*lsh+i]);
                    }
                }
            }
            if(flux_splitter::requires_min_speed) {
                for(int v = 0; v < claw::nequations; ++v) {
                    for(int i = - width + 1; i < width + 1; ++i) {
                        min_speed = std::min(min_speed, eigenvalue[v*lsh+i]);
                    }
                }
            }

            for(int i = 1 - width; i < width + 1; ++i) {
                utils::gemv<CCTK_REAL, false, claw::nequations,
                    claw::nequations>::eval(1, &interface_left_eigenvector[0],
                            claw::nequations, &u[i], lsh, 0,
                            &_M_eta[0][i+width-1], size);
                utils::gemv<CCTK_REAL, false, claw::nequations,
                    claw::nequations>::eval(1, &interface_left_eigenvector[0],
                            claw::nequations, &f[i], lsh, 0,
                            &_M_phi[0][i+width-1], size);
            }

            for(int c = 0; c < claw::nequations; ++c) {
                _M_flux[c] = _M_flux_split.local_flux(delta, cbound,
                        &_M_eta[c][width-1], &_M_phi[c][width-1],
                        interface_eigenvalue[c], max_speed, min_speed,
                        &eigenvalue[c*lsh]);
#ifdef HRSCC_CHECK_FOR_NANS
                {
                    using namespace std;
                    assert(!isnan(_M_flux[c]) && !isinf(_M_flux[c]));
                }
#endif
            }

            utils::gemv<CCTK_REAL, false, claw::nequations, claw::nequations>::
                simple(interface_right_eigenvector, _M_flux, flux);
        }
    private:
        flux_splitter _M_flux_split;

        mutable CCTK_REAL _M_eta[claw::nequations][size];
        mutable CCTK_REAL _M_phi[claw::nequations][size];

        mutable CCTK_REAL _M_flux[claw::nequations];
};

} // namespace

#endif
