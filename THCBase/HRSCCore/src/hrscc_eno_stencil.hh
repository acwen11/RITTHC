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


#ifndef HRSCC_ENO_STENCIL_HH
#define HRSCC_ENO_STENCIL_HH

#include <cctk.h>

namespace hrscc {

//! Stencil for ENO reconstruction
/*!
 *  This class stores the weights needed for the reconstruction on a fixed
 *  stencil of width \e eno_width.
 */
template<int eno_width>
class ENOStencil {
    public:
        //! The width of the stencil
        enum {width = eno_width};
        //! The weights of the stencil
        /*!
         *  The last row is used only for symmetric (W)ENO reconstruction:
         *  \f$ a_{kl} \f$.
         */
        static CCTK_REAL const a[width + 1][width];
};

} // namespace

#endif
