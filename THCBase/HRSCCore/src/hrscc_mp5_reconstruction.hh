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


#ifndef HRSCC_MP5_RECONSTRUCTION_HH
#define HRSCC_MP5_RECONSTRUCTION_HH

#include <algorithm>

#include <cctk.h>

#include <utils.hh>

#include <hrscc_config_par.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! Returns minmod function
template<typename T>
T minmod(T a, T b) {
    return 0.5*(UTILS_SIGN(a)+UTILS_SIGN(b))*std::min(std::abs(a), std::abs(b));
}

//! Returns median of three numbers
template<typename T>
T median(T a, T b, T c) {
    return a + minmod(b-a, c-a);
}

//! MP5 reconstruction
/*!
 *  Currently this is a "stand-alone" reconstruction class, but in the future
 *  we may split it as done for the WENO reconstruction to allow for different
 *  orders of reconstructions.
 */
class MP5Reconstruction {
    public:
        //! width of the stencil used for the reconstruction
        enum {width = 3};
        //! number of grid points needed for the reconstruction
        enum {total_width = 5};
        //! number of data points needed for a generic reconstruction
        /*!
         *  Note that in the case in which only the one-directional
         *  reconstruction is computed with a non-symmetric stencil the number
         *  of points is only total_width.
         *
         *  In the general case a buffer of size \e size should be used.
         */
        enum {size = 6};

        //! reconstruct an array in \f$ v_{i+1/2}^\pm \f$
        /*!
         *  The input array should be centered around \e i so that, at the end,
         *  we reconstruct in 1/2.
         *
         *  Please mind that the number of data points to the left and to the
         *  right of 0 will change  depending on the chosen reconstruction and
         *  direction.
         */
        template<policy::orientation_t sign>
        CCTK_REAL reconstruct(
                //! [in] grid spacing
                CCTK_REAL,
                //! [in] data to reconstruct
                CCTK_REAL const v[size]
                ) const {
            CCTK_REAL const * u = &v[sign == 1];
            CCTK_REAL const up15 = (2*u[2*sign] - 13*u[sign] + 47*u[0] +
                    27*u[-sign] - 3*u[-2*sign]) / 60.0;

            CCTK_REAL const deltam = u[0] - u[sign];
            CCTK_REAL const deltap = u[-sign] - u[0];
            CCTK_REAL const ump = u[0] + minmod(deltap,
                    config::param::mp5_alpha*deltam);

            if((up15 - u[0])*(up15 - ump) < 0) {
                return up15;
            }
            else {
                CCTK_REAL const dm  = u[2*sign] + u[0]       - 2*u[sign];
                CCTK_REAL const d   = u[sign]   + u[-sign]   - 2*u[0];
                CCTK_REAL const dp  = u[0]      + u[-2*sign] - 2*u[-sign];

                CCTK_REAL const dmp = minmod(minmod(4*d-dp, 4*dp-d),
                        minmod(d, dp));
                CCTK_REAL const dmm = minmod(minmod(4*dm - d, 4*d - dm),
                        minmod(dm, d));

                CCTK_REAL const ulc = u[0] + 0.5*deltam + 4.0/3.0*dmm;
                CCTK_REAL const umd = 0.5*(u[0] + u[-sign]) - 0.5*dmp;

                CCTK_REAL const uul = u[0] + config::param::mp5_alpha*deltam;

                CCTK_REAL const umin = std::max(
                        std::min(u[0], std::min(u[-sign], umd)),
                        std::min(u[0], std::min(uul, ulc)));
                CCTK_REAL const umax = std::min(
                        std::max(u[0], std::max(u[-sign], umd)),
                        std::max(u[0], std::max(uul, ulc)));

                return median(umin, up15, umax);
            }
        }
};

} // namespace

#endif
