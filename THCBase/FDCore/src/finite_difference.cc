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


#include <finite_difference.hh>

#include <assert.h>
#include <stdio.h>

///////////////////////////////////////////////////////////////////////////////
// Finite difference stencils (computed with Mathematica)
///////////////////////////////////////////////////////////////////////////////
template<>
CCTK_REAL const fd<1>::stencil[2][2] = {
    {-1.0, 1.0},
    {-1.0, 1.0}
};

template<>
CCTK_REAL const fd<2>::stencil[3][3] = {
    {-1.5,  2.0, -0.5},
    {-0.5,  0,    0.5},
    { 0.5, -2.0,  1.5}
};

template<>
CCTK_REAL const fd<3>::stencil[4][4] = {
    {-1.8333333333333333, 3., -1.5, 0.33333333333333333},
    {-0.33333333333333333, -0.5, 1.0, -0.16666666666666667},
    {0.16666666666666667, -1.0, 0.5, 0.33333333333333333},
    {-0.33333333333333333, 1.5, -3.0, 1.8333333333333333}
};

template<>
CCTK_REAL const fd<4>::stencil[5][5] = {
    {-2.0833333333333335, 4, -3, 1.3333333333333333, -0.25},
    {-0.25,-0.8333333333333334,1.5,-0.5,0.08333333333333333},
    {0.08333333333333333,-0.6666666666666666,0, 0.6666666666666666,
        -0.08333333333333333},
    {-0.08333333333333333,0.5,-1.5,0.8333333333333334,0.25},
    {0.25,-1.3333333333333333,3,-4,2.0833333333333335}
};

template<>
CCTK_REAL const fd<5>::stencil[6][6] = {
    {-2.2833333333333333, 5.0, -5.0, 3.3333333333333333, -1.25, 0.2},
    {-0.2, -1.0833333333333333, 2.0, -1.0, 0.33333333333333333, -0.05},
    {0.05, -0.5, -0.33333333333333333, 1.0, -0.25, 0.033333333333333333},
    {-0.033333333333333333, 0.25, -1.0, 0.33333333333333333, 0.5, -0.05},
    {0.05, -0.33333333333333333, 1.0, -2.0, 1.0833333333333333, 0.2},
    {-0.2, 1.25, -3.3333333333333333, 5.0, -5.0, 2.2833333333333333}
};

template<>
CCTK_REAL const fd<6>::stencil[7][7] = {
    {-2.45,6,-7.5,6.666666666666667,-3.75,1.2, -0.16666666666666666},
    {-0.16666666666666666,-1.2833333333333334,2.5,-1.6666666666666667,
        0.8333333333333334,-0.25,0.03333333333333333},
    {0.03333333333333333,-0.4,-0.5833333333333334,1.3333333333333333,-0.5,
        0.13333333333333333,-0.016666666666666666},
    {-0.016666666666666666,0.15,-0.75,0,0.75,-0.15,0.016666666666666666},
    {0.016666666666666666,-0.13333333333333333,0.5,-1.3333333333333333,
        0.5833333333333334,0.4,-0.03333333333333333},
    {-0.03333333333333333,0.25,-0.8333333333333334,1.6666666666666667,-2.5,
        1.2833333333333334,0.16666666666666666},
    {0.16666666666666666,-1.2,3.75,-6.666666666666667,7.5,-6,2.45}
};

template<>
CCTK_REAL const fd<7>::stencil[8][8] = {
    {-2.5928571428571429, 7.0, -10.5, 11.666666666666667, -8.75, 4.2,
        -1.1666666666666667, 0.14285714285714286},
    {-0.14285714285714286, -1.45, 3.0, -2.5, 1.6666666666666667, -0.75, 0.2,
        -0.023809523809523810},
    {0.023809523809523810, -0.33333333333333333, -0.78333333333333333,
        1.6666666666666667, -0.83333333333333333, 0.33333333333333333,
        -0.083333333333333333, 0.0095238095238095238},
    {-0.0095238095238095238, 0.1, -0.6, -0.25, 1.0, -0.3,
        0.066666666666666667, -0.0071428571428571429},
    {0.0071428571428571429, -0.066666666666666667, 0.3, -1.0, 0.25,
        0.6, -0.1, 0.0095238095238095238},
    {-0.0095238095238095238, 0.083333333333333333, -0.33333333333333333,
        0.83333333333333333, -1.6666666666666667, 0.78333333333333333,
        0.33333333333333333, -0.023809523809523810},
    {0.023809523809523810, -0.2, 0.75, -1.6666666666666667, 2.5, -3.0,
        1.45, 0.14285714285714286},
    {-0.14285714285714286, 1.1666666666666667, -4.2, 8.75,
        -11.666666666666667, 10.5, -7.0, 2.5928571428571429}
};

template<>
CCTK_REAL const fd<8>::stencil[9][9] = {
    {-2.717857142857143,8,-14,18.666666666666668,-17.5,11.2,
        -4.666666666666667,1.1428571428571428,-0.125},
    {-0.125,-1.5928571428571427,3.5,-3.5,2.9166666666666665,-1.75,0.7,
        -0.16666666666666666,0.017857142857142856},
    {0.017857142857142856,-0.2857142857142857,-0.95,2,-1.25,
        0.6666666666666666,-0.25,0.05714285714285714,-0.005952380952380952},
    {-0.005952380952380952,0.07142857142857142,-0.5,-0.45,1.25,-0.5,
        0.16666666666666666,-0.03571428571428571,0.0035714285714285713},
    {0.0035714285714285713,-0.0380952380952381,0.2,-0.8,0,0.8,-0.2,
        0.0380952380952381,-0.0035714285714285713},
    {-0.0035714285714285713,0.03571428571428571,-0.16666666666666666,0.5,
        -1.25,0.45,0.5,-0.07142857142857142,0.005952380952380952},
    {0.005952380952380952,-0.05714285714285714,0.25,-0.6666666666666666,1.25,
        -2,0.95,0.2857142857142857,-0.017857142857142856},
    {-0.017857142857142856,0.16666666666666666,-0.7,1.75,-2.9166666666666665,
            3.5,-3.5,1.5928571428571427,0.125},
    {0.125,-1.1428571428571428,4.666666666666667,-11.2,17.5,
            -18.666666666666668,14,-8,2.717857142857143}
};

///////////////////////////////////////////////////////////////////////////////
// Automatic interface
///////////////////////////////////////////////////////////////////////////////
extern "C" CCTK_REAL adiff_1(
        CCTK_POINTER_TO_CONST _GH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int dir,
        int order) {
    cGH const * const cctkGH = (cGH const * const)_GH;

    int const stride =
        (dir == DIRECTION_X)*(CCTK_GFINDEX3D(cctkGH, 1, 0, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Y)*(CCTK_GFINDEX3D(cctkGH, 0, 1, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Z)*(CCTK_GFINDEX3D(cctkGH, 0, 0, 1)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0));

    int const idx[3] = {i, j, k};
    int dpoint = order/2;
    if(idx[dir] < dpoint) {
        dpoint = idx[dir];
    }
    else if(idx[dir] >= cctkGH->cctk_lsh[dir] - dpoint) {
        dpoint = order + 1 + idx[dir] - cctkGH->cctk_lsh[dir];
    }

    return sdiff_1(grid_function + CCTK_GFINDEX3D(cctkGH, i, j, k), order,
            dpoint, stride);
}

extern "C" CCTK_REAL adiff_x(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order) {
    return adiff_1(cctkGH, grid_function, i, j, k, DIRECTION_X, order);
}


extern "C" CCTK_REAL adiff_y(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order) {
    return adiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Y, order);
}

extern "C" CCTK_REAL adiff_z(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order) {
    return adiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Z, order);
}

///////////////////////////////////////////////////////////////////////////////
// Centered interface
///////////////////////////////////////////////////////////////////////////////
extern "C" CCTK_REAL cdiff_x(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order) {
    return cdiff_1(cctkGH, grid_function, i, j, k, DIRECTION_X, order);
}

extern "C" CCTK_REAL cdiff_y(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order) {
    return cdiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Y, order);
}

extern "C" CCTK_REAL cdiff_z(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order) {
    return cdiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Z, order);
}

extern "C" CCTK_REAL cdiff_1(
        CCTK_POINTER_TO_CONST _GH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int dir,
        int order) {
    cGH const * const cctkGH = (cGH const * const)_GH;
    int const stride =
        (dir == DIRECTION_X)*(CCTK_GFINDEX3D(cctkGH, 1, 0, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Y)*(CCTK_GFINDEX3D(cctkGH, 0, 1, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Z)*(CCTK_GFINDEX3D(cctkGH, 0, 0, 1)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0));
    return sdiff_1(grid_function + CCTK_GFINDEX3D(cctkGH, i, j, k), order,
            order/2, stride);

}

///////////////////////////////////////////////////////////////////////////////
// "Margin" interface
///////////////////////////////////////////////////////////////////////////////
extern "C" CCTK_REAL mdiff_1(
        CCTK_POINTER_TO_CONST _GH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int dir,
        int order,
        int margin) {
    cGH const * const cctkGH = (cGH const * const)_GH;

    int const stride =
        (dir == DIRECTION_X)*(CCTK_GFINDEX3D(cctkGH, 1, 0, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Y)*(CCTK_GFINDEX3D(cctkGH, 0, 1, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Z)*(CCTK_GFINDEX3D(cctkGH, 0, 0, 1)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0));

    int const idx[3] = {i, j, k};
    int const margin_l = margin*(cctkGH->cctk_bbox[2*dir] &&
            idx[dir] >= margin);
    int const margin_r = margin*(cctkGH->cctk_bbox[2*dir+1] &&
            idx[dir] < cctkGH->cctk_lsh[dir] - margin);

    int dpoint = order/2;
    if(idx[dir] - margin_l < dpoint) {
        dpoint = idx[dir] - margin_l;
    }
    else if(idx[dir] + margin_r >= cctkGH->cctk_lsh[dir] - dpoint) {
        dpoint = order - cctkGH->cctk_lsh[dir] + margin_r + 1 + idx[dir];
    }

    return sdiff_1(grid_function + CCTK_GFINDEX3D(cctkGH, i, j, k), order,
            dpoint, stride);
}

extern "C" CCTK_REAL mdiff_x(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order,
        int margin) {
    return mdiff_1(cctkGH, grid_function, i, j, k, DIRECTION_X, order, margin);
}

extern "C" CCTK_REAL mdiff_y(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order,
        int margin) {
    return mdiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Y, order, margin);
}

extern "C" CCTK_REAL mdiff_z(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order,
        int margin) {
    return mdiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Z, order, margin);
}

///////////////////////////////////////////////////////////////////////////////
// Upwind interface
///////////////////////////////////////////////////////////////////////////////
extern "C" CCTK_REAL udiff_1(
        CCTK_POINTER_TO_CONST _GH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int dir,
        int order,
        int shift) {
    cGH const * const cctkGH = (cGH const * const)_GH;

    int const stride =
        (dir == DIRECTION_X)*(CCTK_GFINDEX3D(cctkGH, 1, 0, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Y)*(CCTK_GFINDEX3D(cctkGH, 0, 1, 0)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0)) +
        (dir == DIRECTION_Z)*(CCTK_GFINDEX3D(cctkGH, 0, 0, 1)  -
                              CCTK_GFINDEX3D(cctkGH, 0, 0, 0));

    int dpoint = order/2 + shift;

    return sdiff_1(grid_function + CCTK_GFINDEX3D(cctkGH, i, j, k), order,
            dpoint, stride);
}

extern "C" CCTK_REAL udiff_x(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order,
        int shift) {
    return udiff_1(cctkGH, grid_function, i, j, k, DIRECTION_X, order, shift);
}

extern "C" CCTK_REAL udiff_y(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order,
        int shift) {
    return udiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Y, order, shift);
}

extern "C" CCTK_REAL udiff_z(
        CCTK_POINTER_TO_CONST cctkGH,
        CCTK_REAL const * grid_function,
        int i, int j, int k,
        int order,
        int shift) {
    return udiff_1(cctkGH, grid_function, i, j, k, DIRECTION_Z, order, shift);
}

///////////////////////////////////////////////////////////////////////////////
// Low-level C interface
///////////////////////////////////////////////////////////////////////////////
extern "C" CCTK_REAL sdiff_1(
        CCTK_REAL const * grid_function_ijk,
        int order,
        int dpoint,
        int stride) {
    assert(dpoint >= 0 && dpoint <= order);
    stencil_t const dop[8] = {&fd<1>::sdiff, &fd<2>::sdiff, &fd<3>::sdiff,
        &fd<4>::sdiff, &fd<5>::sdiff, &fd<6>::sdiff, &fd<7>::sdiff,
        &fd<8>::sdiff};
    return (*dop[order-1])(grid_function_ijk, dpoint, stride);
}
