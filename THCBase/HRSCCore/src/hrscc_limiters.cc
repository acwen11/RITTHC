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


#include <algorithm>
#include <cmath>

#include <hrscc_config_par.hh>
#include <hrscc_limiters.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {
namespace limiters {

CCTK_REAL minmod::eval(CCTK_REAL r) {
    CCTK_REAL const th = hrscc::config::param::minmod_theta;
    return std::max(0.0, std::min(th*r, std::min(0.5*(1.+r), th)));
}

CCTK_REAL superbee::eval(CCTK_REAL r) {
    return std::max(0.0, std::max(std::min(2*r, 1.0), std::min(r, 2.0)));
}

CCTK_REAL vanleer::eval(CCTK_REAL r) {
    return (r + std::abs(r)) / (1 + std::abs(r));
}

}
}
