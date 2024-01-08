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


#ifndef HRSCC_COMPONENT_SPLIT_HH
#define HRSCC_COMPONENT_SPLIT_HH

#include <algorithm>
#include <cmath>
#include <limits>

#include <cctk.h>

namespace hrscc {

//! performs a component-wise split of a system of equations
/*!
 *  This will reconstruct the fluxes by decomposing the system of equations
 *  into its components and using the given flux-splitter to reconstruct the
 *  fluxes of each one.
 *
 *  \tparam claw_t conservation law that we want to split
 *  \tparam flux_splitter_t flux-splitter to use for the reconstruction
 *
 *  Note that flux-splitters which work using the wave speeds are only
 *  supported in the scalar case
 */
template<typename claw_t, typename flux_splitter_t>
class ComponentSplit {
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
        static bool const requires_interface_eigenvectors = false;
        //! if true the method requires all the eigenvectors
        static bool const requires_all_eigenvectors = false;

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
                CCTK_REAL const *,// interface_left_eigenvector,
                //! [in] right eigenvectors at the reconstruction point, stored
                //! in row-major order
                CCTK_REAL const *,// interface_right_eigenvector,
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

            for(int eq = 0; eq < claw::nequations; ++eq) {
                flux[eq] = _M_flux_split.local_flux(delta, cbound,
                    &u[eq*lsh], &f[eq*lsh], interface_eigenvalue[eq],
                    max_speed, min_speed, &eigenvalue[eq*lsh]);
            }
        }
    private:
        flux_splitter _M_flux_split;
};

} // namespace

#endif
