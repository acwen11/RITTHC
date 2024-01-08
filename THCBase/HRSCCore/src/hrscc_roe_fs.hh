//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2012, David Radice <david.radice@aei.mpg.de>
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


#ifndef HRSCC_ROE_FS_HH
#define HRSCC_ROE_FS_HH

#include <cmath>
#include <cstdlib>

#include <cctk.h>

#include <utils.hh>

#include <hrscc_config_par.hh>
#include <hrscc_lax_friedrichs_fs.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! Dummy class for trivial (=no) entropy fix
/*!
 *  I hope you know what you are doing if you use this one
 */
class NoEntropyFix {
    public:
        static bool const requires_min_speed = false;
        static bool const requires_max_speed = false;

        CCTK_REAL local_flux(
                CCTK_REAL delta,
                CCTK_REAL cbound,
                CCTK_REAL const * u,
                CCTK_REAL const * f,
                CCTK_REAL interface_speed,
                CCTK_REAL max_speed,
                CCTK_REAL min_speed,
                CCTK_REAL const * speed
                ) const;
};


//! Class for Roe flux-splitting
/*!
 *  \tparam reconstructor_t the used reconstructor
 *  \tparam entropy_fix_t use monotonicity-based rarefaction-wave detection
 */
template<typename reconstructor_t, typename entropy_fix_t>
class RoeFS {
    public:
        typedef reconstructor_t reconstructor;
        typedef entropy_fix_t entropy_fix;

        //! width of the stencil used for the reconstruction
        enum {width = reconstructor::width};
        //! number of points needed for the reconstruction
        enum {size = 2*width};

        //! if true this means that the method uses the speed at the interface
        static bool const requires_interface_speed = true;
        //! if true this means that the method uses the maximum (local) speed
        static bool const requires_max_speed = entropy_fix::requires_max_speed;
        //! if true this means that the method uses the minimum (local) speed
        static bool const requires_min_speed = entropy_fix::requires_min_speed;
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
                //! [in] characteristic speed at the interface
                CCTK_REAL interface_speed,
                //! [in] max speed to use in the local variant of the method
                CCTK_REAL max_speed,
                //! [in] min speed to use in the local variant of the method
                CCTK_REAL min_speed,
                //! [in] speeds on the stencil
                CCTK_REAL const speed[size]
                ) const {
            if(utils::is_same_type<entropy_fix, NoEntropyFix>::value || (
                       std::abs(interface_speed) > config::param::speed_eps*
                            cbound
                    && utils::sequence<CCTK_REAL, size>::
                            is_strictly_monotonic(&u[1-width])
                    && utils::sequence<CCTK_REAL, size>::
                            is_strictly_monotonic(&f[1-width])
                    )) {
                if(interface_speed > 0) {
                    return _M_reconstructor.
                        reconstructor::template reconstruct<policy::minus>(
                                delta, f);
                }
                else {
                    return _M_reconstructor.
                        reconstructor::template reconstruct<policy::plus>(
                                delta, f);
                }
            }
            else {
                return _M_entropy_fix.local_flux(delta, cbound, u, f,
                        interface_speed, max_speed, min_speed, speed);
            }
        }
    private:
        reconstructor _M_reconstructor;
        entropy_fix _M_entropy_fix;
};

} // namespace

#endif
