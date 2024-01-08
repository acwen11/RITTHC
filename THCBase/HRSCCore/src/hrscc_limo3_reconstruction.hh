//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2013, David Radice <david.radice@aei.mpg.de>
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


#ifndef HRSCC_LIMO3_RECONSTRUCTION_HH
#define HRSCC_LIMO3_RECONSTRUCTION_HH

#include <algorithm>

#include <cctk.h>

#include <utils.hh>

#include <hrscc_config_par.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! LimO3 reconstruction
class LimO3Reconstruction {
    public:
        //! width of the stencil used for the reconstruction
        enum {width = 2};
        //! number of grid points needed for the reconstruction
        enum {total_width = 3};
        //! number of data points needed for a generic reconstruction
        /*!
         *  Note that in the case in which only the one-directional
         *  reconstruction is computed with a non-symmetric stencil the number
         *  of points is only total_width.
         *
         *  In the general case a buffer of size \e size should be used.
         */
        enum {size = 4};

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
                CCTK_REAL delta,
                //! [in] data to reconstruct
                CCTK_REAL const v[size]
                ) const {
            CCTK_REAL const * u = &v[(sign == policy::plus)];

            CCTK_REAL const deltap = u[-sign] - u[0];
            if(std::abs(deltap) <= std::numeric_limits<CCTK_REAL>::epsilon()) {
                return u[0];
            }
            CCTK_REAL const deltam = u[0] - u[sign];
            CCTK_REAL const theta = deltam/deltap;

            CCTK_REAL const eta = (utils::pow<2>(deltam)+utils::pow<2>(deltap))
                        / (utils::pow<2>(config::param::limo3_r*delta));
            CCTK_REAL const chi = std::max(0.0, std::min(1.0,
                        0.5 + (eta - 1.0)/(2*config::param::limo3_eps)));

            CCTK_REAL const P3 = (2.0 + theta)/3.0;
            CCTK_REAL phi = 0;
            if(theta >= 0) {
                phi = std::max(0.0, std::min(P3, std::min(2*theta, 1.6)));
            }
            else {
                phi = std::max(0.0, std::min(P3, -theta/2.0));
            }

            return u[0] + deltap/2.0*(P3 + chi*(phi - P3));
        }
};

} // namespace

#endif
