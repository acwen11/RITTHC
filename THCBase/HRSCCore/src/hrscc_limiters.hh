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


#ifndef HRSCC_LIMITERS_HH
#define HRSCC_LIMITERS_HH

#include <cctk.h>

namespace hrscc {
//! Limiters used for the TVD reconstruction
namespace limiters {

//! MinMod limiter
class minmod {
    public:
        static CCTK_REAL eval(
                CCTK_REAL r             //!< [in] slope ratios
                );
};

//! SuperBee limiter
class superbee {
    public:
        static CCTK_REAL eval(
                CCTK_REAL r             //!< [in] slope ratios
                );
};

//! Van Leer limiter
class vanleer {
    public:
        static CCTK_REAL eval(
                CCTK_REAL r             //!< [in] slope ratios
                );
};

} // namespace limiters
} // namespace hrscc

#endif
