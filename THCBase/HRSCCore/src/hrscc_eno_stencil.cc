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

#include <cctk.h>

#include <hrscc_eno_stencil.hh>

namespace hrscc {

template<>
CCTK_REAL const ENOStencil<2>::a[3][2] = {
    {-0.5, 1.5},
    {0.5,  0.5},
    {1.5, -0.5}
};

template<>
CCTK_REAL const ENOStencil<3>::a[4][3] = {
    {  2.0  / 6.0,  -7.0 / 6.0, 11.0 / 6.0  },
    { -1.0  / 6.0,   5.0 / 6.0,  2.0  / 6.0 },
    {  2.0  / 6.0,   5.0 / 6.0, -1.0  / 6.0 },
    { 11.0  / 6.0,  -7.0 / 6.0,  2.0  / 6.0 }
};

template<>
CCTK_REAL const ENOStencil<4>::a[5][4] = {
    { -6.0 / 24.0,  26.0 / 24.0, -46.0 / 24.0, 50.0 / 24.0 },
    {  2.0 / 24.0, -10.0 / 24.0,  26.0 / 24.0,  6.0 / 24.0 },
    { -2.0 / 24.0,  14.0 / 24.0,  14.0 / 24.0, -2.0 / 24.0 },
    {  6.0 / 24.0,  26.0 / 24.0, -10.0 / 24.0,  2.0 / 24.0 },
    { 50.0 / 24.0, -46.0 / 24.0,  26.0 / 24.0, -6.0 / 24.0 }
};

} // namespace
