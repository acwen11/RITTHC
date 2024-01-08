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


#ifndef HRSCC_WENO_LIMITER_HH
#define HRSCC_WENO_LIMITER_HH


#include <stdio.h>
#include <cctk.h>

#include <hrscc_config_par.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! WENO smoothness indicator limiter
/*!
 *  In some applications it is necessary to avoid over-adaptation of the stencil
 *  in WENO methods. See Taylor et al. J. Comput. Phys. 223 (2007) 384-397.
 */
template<int eno_width, int weno_width, policy::weno_limiter_t type>
class WENOLimiter;

//! WENO smoothens indicator limiter specialized to absolute limiting
template<int eno_width, int weno_width>
class WENOLimiter<eno_width, weno_width, policy::absolute> {
    public:
        //! Apply the smoothness indicator limiter
        static void eval(
                //! [in,out] smoothenss indicator coefficients
                CCTK_REAL IS[weno_width]
                ) {
            for(int i = 0; i < weno_width; ++i) {
                IS[i] = IS[i] * static_cast<CCTK_REAL>(IS[i] >=
                        config::param::weno_alpha);
            }
            if(weno_width > eno_width) {
                CCTK_REAL maxIS = IS[0];
                for(int i = 1; i < weno_width; ++i) {
                    maxIS = maxIS > IS[i] ? maxIS : IS[i];
                }
                IS[eno_width] = maxIS;
            }
        }
};

//! WENO smoothens indicator limiter specialized to dummy limiting
template<int eno_width, int weno_width>
class WENOLimiter<eno_width, weno_width, policy::dummy> {
    public:
        //! Apply the smoothness indicator limiter
        static void eval(
                //! [in,out] smoothenss indicator coefficients
                CCTK_REAL IS[weno_width]
                ) {
            if(weno_width > eno_width) {
                CCTK_REAL maxIS = IS[0];
                for(int i = 1; i < weno_width; ++i) {
                    maxIS = maxIS > IS[i] ? maxIS : IS[i];
                }
                IS[eno_width] = maxIS;
            }
        }
};

//! WENO smoothens indicator limiter specialized for unstable limiting
template<int eno_width, int weno_width>
class WENOLimiter<eno_width, weno_width, policy::unstable> {
    public:
        //! Apply the smoothness indicator limiter
        static void eval(
                //! [in,out] smoothenss indicator coefficients
                CCTK_REAL IS[weno_width]
                ) {
            for(int i = 0; i < weno_width; ++i) {
                IS[i] = 0.0;
            }
        }
};


//! WENO smoothens indicator limiter specialized to relative limiting
template<int eno_width, int weno_width>
class WENOLimiter<eno_width, weno_width, policy::relative> {
    public:
        //! Apply the smoothness indicator limiter
        static void eval(
                //! [in,out] smoothenss indicator coefficients
                CCTK_REAL IS[weno_width]
                ) {
            CCTK_REAL maxIS = IS[0];
            CCTK_REAL minIS = IS[0];
            for(int i = 1; i < weno_width; ++i) {
                maxIS = maxIS > IS[i] ? maxIS : IS[i];
                minIS = minIS < IS[i] ? minIS : IS[i];
            }

            CCTK_REAL RIS = maxIS / (config::param::weno_eps + minIS);

            if(RIS < config::param::weno_alpha) {
                maxIS = 0;
                for(int i = 0; i < weno_width; ++i) {
                    IS[i] = 0;
                }
            }
            if(weno_width > eno_width) {
                IS[eno_width] = maxIS;
            }
        }
};

#undef HRSCC_WENO_LIMITER_BASE

} // namespace

#endif
