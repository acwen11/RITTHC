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


#ifndef HRSCC_LAX_FRIEDRICHS_HH
#define HRSCC_LAX_FRIEDRICHS_HH

#include <cstdlib>

#include <cctk.h>

#include <hrscc_config_par.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! Class for Lax-Friedrichs flux-splitting
/*!
 *  \tparam reconstructor_t the used reconstructor
 *  \tparam local use the local (as opposed to global) Lax-Friedrichs method
 */
template<typename reconstructor_t, bool local>
class LaxFriedrichsFS {
    public:
        typedef reconstructor_t reconstructor;

        //! width of the stencil used for the reconstruction
        enum {width = reconstructor::width};
        //! number of points needed for the reconstruction
        enum {size = 2*width};

        //! we don't need the speed at the reconstruction point
        static bool const requires_interface_speed = false;
        //! if true this means that the method uses the min speed
        static bool const requires_min_speed = local;
        //! if true this means that the method uses the max speed
        static bool const requires_max_speed = local;
        //! if true this means that the method uses the speed over the stencil
        static bool const requires_speed = false;

        //! reconstruct the fluxes in \f$ x_{i+1/2} \f$
        /*!
         *  The input arrays should be centered around \e i so that, at the end,
         *  we reconstruct in 1/2.
         *
         *  This function has to have the same interface on all the
         *  flux-splitters. If more generality is needed we could consider to
         *  switch to a model in which a polymorphic argument table is used.
         */
        CCTK_REAL local_flux(
                //! [in] grid spacing
                CCTK_REAL delta,
                //! [in] The maximum propagation speed of a signal
                CCTK_REAL cbound,
                //! [in] conserved variable around the reconstruction point
                CCTK_REAL const u[size],
                //! [in] fluxes around the reconstruction point
                CCTK_REAL const f[size],
                //! [in] dummy argument
                CCTK_REAL, // interface_speed,
                //! [in] max speed to use in the local variant of the method
                CCTK_REAL max_speed,
                //! [in] min speed to use in the local variant of the method
                CCTK_REAL min_speed,
                //! [in] speeds on the stencil
                CCTK_REAL const * //speed[size]
                ) const {
            CCTK_REAL flux = 0;
            CCTK_REAL * v = &_M_splitted_flux[width-1];

            CCTK_REAL s = cbound;
            if(local) {
                s = std::max(std::abs(min_speed), std::abs(max_speed));
                // Add dissipation close to the weakly hyperbolic limit
                if(s < config::param::speed_eps*cbound) {
                    s = cbound;
                }
            }

            for(int i = - width + 1; i < width + 1; ++i) {
                v[i] = 0.5*(f[i] + s*u[i]);
            }
            flux += _M_reconstructor.
                reconstructor::template reconstruct<policy::minus>(delta, v);

            for(int i = - width + 1; i < width + 1; ++i) {
                v[i] = 0.5*(f[i] - s*u[i]);
            }
            flux += _M_reconstructor.
                reconstructor::template reconstruct<policy::plus>(delta, v);

            return flux;
        }
    private:
        reconstructor _M_reconstructor;
        mutable CCTK_REAL _M_splitted_flux[size];
};

} // namespace

#endif
