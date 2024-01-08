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


#ifndef HRSCC_TVD_RECONSTRUCTION_HH
#define HRSCC_TVD_RECONSTRUCTION_HH

#include <cctk.h>

#include <hrscc_limiters.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! Core class for TVD reconstruction
template<typename limiter>
class TVDReconstruction {
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
            CCTK_REAL const dp = u[1] - u[0];
            CCTK_REAL const dm = u[0] - u[-1];
            CCTK_REAL r;
            if(std::abs(dp) > std::numeric_limits<CCTK_REAL>::epsilon()) {
                r = dm / dp;
            }
            else {
                r = 0;
            }
            return u[0] - 0.5*sign*limiter::eval(r)*(u[1] - u[0]);
        }
};

//! Alias for the classical minmod-based reconstruction
typedef TVDReconstruction<limiters::minmod> MinModReconstruction;
//! Alias for the TVD reconstruction with Van Leer's limiter
typedef TVDReconstruction<limiters::vanleer> VanLeerReconstruction;
//! Alias for the TVD reconstruction SuperBee limiter
typedef TVDReconstruction<limiters::superbee> SuperBeeReconstruction;

} // namespace

#endif
