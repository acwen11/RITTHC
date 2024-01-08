/*
 *  FDCore: generic finite-difference operators for Cactus
 *  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include <cctk.h>

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************/
/* Handle for directional derivatives                                        */
/*****************************************************************************/
#define DIRECTION_X 0
#define DIRECTION_Y 1
#define DIRECTION_Z 2

/*****************************************************************************/
/* Finite differencing operators: high-level interface                       */
/*****************************************************************************/
/* Auto. diff: compute the first derivative at the given order of accuracy */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL adiff_x(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order);                         /* order of the derivative */
/* Auto. diff: compute the first derivative at the given order of accuracy */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL adiff_y(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order);                         /* order of the derivative */
/* Auto. diff: compute the first derivative at the given order of accuracy */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL adiff_z(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order);                         /* order of the derivative */
/* Auto. diff: compute the first derivative at the given order of accuracy */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL adiff_1(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int direction,                      /* direction of the derivative */
        int order);                         /* order of the derivative */

/* Centered diff: compute the first derivative at the given order of accuracy */
/* using a centered stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL cdiff_x(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order);                         /* order of the derivative */
/* Centered diff: compute the first derivative at the given order of accuracy */
/* using a centered stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL cdiff_y(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order);                         /* order of the derivative */
/* Centered diff: compute the first derivative at the given order of accuracy */
/* using a centered stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL cdiff_z(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order);                         /* order of the derivative */
/* Centered diff: compute the first derivative at the given order of accuracy */
/* using a centered stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL cdiff_1(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int direction,                      /* direction of the derivative */
        int order);                         /* order of the derivative */

/* Margin diff: compute the first derivative at the given order of accuracy */
/* excluding the first/last "margin" points of the physical grid from */
/* the stencil. */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL mdiff_x(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order,                          /* order of the derivative */
        int margin);                        /* number ofexcluded points */
/* Margin diff: compute the first derivative at the given order of accuracy */
/* excluding the first/last "margin" points of the physical grid from */
/* the stencil. */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL mdiff_y(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order,                          /* order of the derivative */
        int margin);                        /* number ofexcluded points */

/* Margin diff: compute the first derivative at the given order of accuracy */
/* excluding the first/last "margin" points of the physical grid from */
/* the stencil. */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL mdiff_z(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order,                          /* order of the derivative */
        int margin);                        /* number ofexcluded points */
/* Margin diff: compute the first derivative at the given order of accuracy */
/* excluding the first/last "margin" points of the physical grid from */
/* the stencil. */
/* This will take care of using a one-sided stencil when needed */
CCTK_REAL mdiff_1(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int direction,                      /* direction of the derivative */
        int order,                          /* order of the derivative */
        int margin);                        /* number ofexcluded points */

/* Upwind diff: compute the first derivative at the given order of accuracy */
/* using a left or right biased stencil */
/* The shift has to be smaller than the size of the stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL udiff_x(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order,                          /* order of the derivative */
        int shift);                         /* stencil shift */
/* Upwind diff: compute the first derivative at the given order of accuracy */
/* using a left or right biased stencil */
/* The shift has to be smaller than the size of the stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL udiff_y(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order,                          /* order of the derivative */
        int shift);                         /* stencil shift */
/* Upwind diff: compute the first derivative at the given order of accuracy */
/* using a left or right biased stencil */
/* The shift has to be smaller than the size of the stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL udiff_z(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int order,                          /* order of the derivative */
        int shift);                         /* stencil shift */
/* Upwind diff: compute the first derivative at the given order of accuracy */
/* using a centered stencil */
/* The shift has to be smaller than the size of the stencil */
/* This will NOT take care of handling boundary points! */
CCTK_REAL udiff_1(
        CCTK_POINTER_TO_CONST _GH,          /* grid hierarchy */
        CCTK_REAL const * gf,               /* function to differentiate */
        int i, int j, int k,                /* grid point */
        int direction,                      /* direction of the derivative */
        int order,                          /* order of the derivative */
        int shift);                         /* stencil shift */

/*****************************************************************************/
/* Finite differencing operators: low-level interface                        */
/*****************************************************************************/
/* Static diff: compute the finite difference using the given stencil */
CCTK_REAL sdiff_1(
        CCTK_REAL const * gf_ijk,           /* grid func. at the diff. point */
        int order,                          /* which stencil to use */
        int dpoint,                         /* diff. point on the stencil */
        int stride);                        /* stride */

#ifdef __cplusplus
}
#endif

#endif
