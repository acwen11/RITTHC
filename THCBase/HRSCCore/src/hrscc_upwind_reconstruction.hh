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


#ifndef HRSCC_UPWIND_RECONSTRUCTION_HH
#define HRSCC_UPWIND_RECONSTRUCTION_HH

#include <cctk.h>

#include <hrscc_typedefs.hh>

namespace hrscc {

//! First order upwind reconstruction
/*!
 *  Simple 1st order accurate upwind reconstruction:
 *  \f$
        v^+_{i+1/2} := v_{i+1}, \quad
        v^-_{i+1/2} := v_{i}.
    \f$
 */
class UpwindReconstruction {
    public:
        //! width of the stencil used for the reconstruction
        enum {width = 1};
        //! number of grid points needed for the reconstruction
        enum {total_width = 2};
        //! number of data points needed for a generic reconstruction
        /*!
         *  Note that in the case in which only the one-directional
         *  reconstruction is computed with a non-symmetric stencil the number
         *  of points is only total_width.
         *
         *  In the general case a buffer of size \e size should be used.
         */
        enum {size = 3};

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
            return v[(sign == policy::plus)];
        }
};

} // namespace

#endif
