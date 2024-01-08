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

#include <hrscc_weno_stencil.hh>

namespace hrscc {

// r = 2
template<>
CCTK_REAL const WENOStencil<2, policy::standard, policy::order>::C[2] =
    {1.0/3.0, 2.0/3.0};

template<>
CCTK_REAL const WENOStencil<2, policy::symmetric, policy::order>::C[3] =
    {1.0/6.0, 4.0/6.0, 1.0/6.0};

// r = 3
template<>
CCTK_REAL const WENOStencil<3, policy::standard, policy::order>::C[3] =
    {1.0/10.0, 6.0/10.0, 3.0/10.0};

template<>
CCTK_REAL const WENOStencil<3, policy::symmetric, policy::order>::C[4] =
    {1.0/20.0, 9.0/20.0, 9.0/20.0, 1.0/20.0};

template<>
CCTK_REAL const WENOStencil<3, policy::symmetric, policy::bandwidth>::C[4] =
    {0.094647545896, 0.428074212384, 0.408289331408, 0.068988910311};

// r = 4
template<>
CCTK_REAL const WENOStencil<4, policy::standard, policy::order>::C[4] =
    {1.0/35.0, 12.0/35.0, 18.0/35.0, 4.0/35.0};

template<>
CCTK_REAL const WENOStencil<4, policy::symmetric, policy::order>::C[5] =
    {1.0/70.0, 16.0/70.0, 36.0/70.0, 16.0/70.0, 1.0/70.0};

template<>
CCTK_REAL const WENOStencil<4, policy::symmetric, policy::bandwidth>::C[5] =
    {0.040195483373, 0.249380000671, 0.480268625626, 0.200977547673,
        0.029178342658};

} // namespace
