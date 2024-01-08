// Dummy cctk.h include file

#include <stdio.h>

#ifndef CCTK_H
#define CCTK_H

#define CCTK_INT int
#define CCTK_REAL double

#define CCTK_POINTER_TO_CONST void const *

#define CCTK_GFINDEX3D(cctkGH,I,J,K)                                           \
        ((I) + cctkGH->cctk_lsh[0]*((J) + cctkGH->cctk_lsh[1]*(K)))

#define CCTK_DELTA_SPACE(I) (cctkGH->cctk_delta_space[(I)])

#define CCTK_WARN_ABORT 0
#define CCTK_WARN(A,B) printf(B);

typedef struct _cGH {
    int cctk_lsh[3];
    int cctk_bbox[6];
    CCTK_REAL cctk_delta_space[3];
    CCTK_REAL cctk_levfac[3];
} cGH;

#endif
