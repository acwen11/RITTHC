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


#ifndef HRSCC_U5_RECONSTRUCTION_HH
#define HRSCC_U5_RECONSTRUCTION_HH

#include <cctk.h>

#include <hrscc_typedefs.hh>

namespace hrscc {

//! Unlimited 5th order reconstruction
/*!
 *  Currently this is a "stand-alone" reconstruction class, but in the future
 *  we may split it as done for the WENO reconstruction to allow for different
 *  orders of reconstructions.
 */
class U5Reconstruction {
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
            return (2*u[2*sign] - 13*u[sign] + 47*u[0] +
                    27*u[-sign] - 3*u[-2*sign]) / 60.0;
        }
};

} // namespace

#endif
