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


#ifndef HRSCC_WENO_WEIGHTS_HH
#define HRSCC_WENO_WEIGHTS_HH

#include <cctk.h>

namespace hrscc {

//! Weights for the non-linear smoothness indicadors
template<int eno_width>
class WENOWeights {
    public:
        //! The weights
        /*!
         *  The last row is only used for symmetric WENO reconstruction:
         *  \f$ d_{kml} \f$
         */
        static CCTK_REAL const d[eno_width + 1][eno_width - 1][eno_width];
};

} // namespace

#endif
