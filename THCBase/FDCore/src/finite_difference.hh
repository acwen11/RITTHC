//  FDCore: generic finite-difference operators for Cactus
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


#ifndef FINITE_DIFFERENCE_HH
#define FINITE_DIFFERENCE_HH

#include <cctk.h>
#include <finite_difference.h>

// stencil type
typedef CCTK_REAL (*stencil_t)(CCTK_REAL const * const, int, int);

// Low-level FD interface: the user should not use this directly
template<int order>
class fd {
    public:
        // Width of the stencil
        enum {width = order + 1};

        // Weights of the stencil: so that stencil[p][i] is the weight of
        // at the ith point when the derivative is computed at the pth point
        // of the stencil, for example, if order == 6, stencil[3][i] gives
        // the weights for the centered finite-difference
        static CCTK_REAL const stencil[width][width];

        // Static diff: compute the finite difference using the given stencil
        static CCTK_REAL sdiff(
                // Grid function to differentiate at the diff. point
                CCTK_REAL const * const grid_function_ijk,
                // Point in the stensil in which to compute the derivative
                int dpoint,
                // Vector stride
                int stride) {
            CCTK_REAL d = 0;
            for(int p = 0; p < width; ++p) {
                d += stencil[dpoint][p] * grid_function_ijk[(p-dpoint)*stride];
            }
            return d;
        }
};

#endif
