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


#ifndef HRSCC_MACRO_HH
#define HRSCC_MACRO_HH

#include <limits>

#include <utils.hh>

//! machine precision
#ifndef HRSCC_EPSILON
#define HRSCC_EPSILON                                                         \
    1e6*std::numeric_limits<CCTK_REAL>::epsilon()
#endif

#ifdef HRSCC_DEBUG

#ifndef HRSCC_CHECK_EIGENVECTORS
//! check that the left and right eigenvectors are orthonormal
#define HRSCC_CHECK_EIGENVECTORS
#endif

#ifndef HRSCC_CHECK_FOR_NANS
//! check for NaNs and Infs in the computed RHSs
#define HRSCC_CHECK_FOR_NANS
#endif

#endif

#endif
