#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "finite_difference.h"

#define NPOINTS 100
#define FDORDER 8

#ifndef M_PI
#define M_PI 3.1415926535897932384
#endif

CCTK_REAL func(CCTK_REAL x) {
    //return sin(20*M_PI*x);
    return 23*pow(x,8) + 64*pow(x,6) - 480*pow(x,4) + 720*pow(x,2) - 120;
}

CCTK_REAL dfunc(CCTK_REAL x) {
    //return 20*M_PI*cos(20*M_PI*x);
    return 184*pow(x,7) + 384*pow(x,5) - 1920*pow(x,3) + 1440*x;
}

int main(void) {
    CCTK_REAL * f;
    CCTK_REAL * df;
    cGH * cctkGH = malloc(sizeof(cGH));

    cctkGH->cctk_lsh[0] = 1;
    cctkGH->cctk_lsh[1] = 1;
    cctkGH->cctk_lsh[2] = NPOINTS;

    cctkGH->cctk_delta_space[0] = 1.0;
    cctkGH->cctk_delta_space[1] = 1.0;
    cctkGH->cctk_delta_space[2] = 2.0 / (double)NPOINTS;

    cctkGH->cctk_bbox[0] = 1;
    cctkGH->cctk_bbox[1] = 1;
    cctkGH->cctk_bbox[2] = 1;
    cctkGH->cctk_bbox[3] = 1;
    cctkGH->cctk_bbox[4] = 1;
    cctkGH->cctk_bbox[5] = 1;

    f  = malloc(NPOINTS*sizeof(CCTK_REAL));
    df = malloc(NPOINTS*sizeof(CCTK_REAL));

    for(int i = 0; i < NPOINTS; ++i) {
        *(f + i) = func(-1 + (double)i * cctkGH->cctk_delta_space[2]);
    }
    for(int i = 0; i < NPOINTS; ++i) {
        *(df + i) = mdiff_z(cctkGH, f, 0, 0, i, FDORDER, 2) /
            cctkGH->cctk_delta_space[2];
    }

    CCTK_REAL delta;
    CCTK_REAL e_inf = 0;

    CCTK_REAL x;
    FILE * outfile = fopen("test_finite_difference.dat", "w");
    fprintf(outfile, "# 1: x, 2: f, 3: df, 4: FD\n");
    for(int i = 0; i < NPOINTS; ++i) {
        x = -1 + (double)i * cctkGH->cctk_delta_space[2];
        delta = fabs(dfunc(x) - df[i]);
        e_inf = (e_inf > delta) ? (e_inf) : (delta);
        fprintf(outfile, "%f %f %f %f\n", x, f[i], dfunc(x), df[i]);
    }
    fclose(outfile);

    printf("e_inf = %f\n", e_inf);

    free(df);
    free(f);
    free(cctkGH);
}
