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


// Weights from Martin et al. J. of Comput. Phys 220 (2006) 270-289

#include <cmath>

#include <cctk.h>

#include <hrscc_weno_weights.hh>

#define beta  (std::sqrt(13.0/12.0))
#define gamma (std::sqrt(781.0/720.0))

namespace hrscc {

template<>
CCTK_REAL const WENOWeights<2>::d[3][1][2] = {
    {{1.0, -1.0}},
    {{1.0, -1.0}},
    {{1.0, -1.0}}
};

template<>
CCTK_REAL const WENOWeights<3>::d[4][2][3] = {
    // k = 0
    {
        {1.0/2.0, -4.0/2.0, 3.0/2.0},
        {beta, -2*beta, beta}
    },
    // k = 1
    {
        {-1.0/2.0, 0, 1.0/2.0},
        {beta, -2*beta, beta}
    },
    // k = 2
    {
        {-3.0/2.0, 4.0/2.0, -1.0/2.0},
        {beta, -2*beta, beta}
    },
    // k = 3
    {
        {-5.0/2.0, 8.0/2.0, -3.0/2.0},
        {beta, -2*beta, beta}
    }
};

template<>
CCTK_REAL const WENOWeights<4>::d[5][3][4] = {
    // k = 0
    {
        {-2.0/6.0, 9.0/6.0, -18.0/6.0, 11.0/6.0},
        {-beta, 4*beta, -5*beta, 2*beta},
        {-gamma, 3*gamma, -3*gamma, gamma}
    },
    // k = 1
    {
        {1.0/6.0, -6.0/6.0, 3.0/6.0, 2.0/6.0},
        {0, beta, -2*beta, beta},
        {-gamma, 3*gamma, -3*gamma, gamma}
    },
    // k = 2
    {
        {-2.0/6.0, -3.0/6.0, 6.0/6.0, -1.0/6.0},
        {beta, -2*beta, beta, 0},
        {-gamma, 3*gamma, -3*gamma, gamma}
    },
    // k = 3
    {
        {-11.0/6.0, 18.0/6.0, -9.0/6.0, 2.0/6.0},
        {2*beta, -5*beta, 4*beta, -beta},
        {-gamma, 3*gamma, -3*gamma, gamma}
    },
    // k = 4
    {
        {-26.0/6.0, 57.0/6.0, -42.0/6.0, 11.0/6.0},
        {3*beta, -8*beta, 7*beta, -2*beta},
        {-gamma, 3*gamma, -3*gamma, gamma}
    }
};

#undef beta
#undef gamma

} // namespace
