//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
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


#include <cassert>
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"
#include "utils.hh"

using namespace thc::m1;
using namespace std;
using namespace utils;

#define NDIM 3

#define GFINDEX1D(__k, ig, iv) \
    ((iv) + (ig)*5 + (__k)*(5*ngroups*nspecies))

#define PINDEX1D(ig, iv) \
    ((iv) + (ig)*5)

extern "C" void THC_M1_Smoothing(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_Smoothing");
    }

    // Aliases for the cons pointers
    CCTK_REAL ** cons = new CCTK_REAL * [ngroups*nspecies*5];
    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
        int const i4D = CCTK_VectGFIndex3D(cctkGH, 0, 0, 0, ig);
        cons[PINDEX1D(ig, 0)] = &rN[i4D];
        cons[PINDEX1D(ig, 1)] = &rFx[i4D];
        cons[PINDEX1D(ig, 2)] = &rFy[i4D];
        cons[PINDEX1D(ig, 3)] = &rFz[i4D];
        cons[PINDEX1D(ig, 4)] = &rE[i4D];
    }

#pragma omp parallel
    {
        for (int dir = 0; dir < NDIM; ++dir) {
            int index[3];
            int lsh[3];

            // We have to leave as the most internal loop the one on the
            // direction. For this reason we will remap the usual indices
            // i,j,k into different points of index[:].
            int ii, ij, ik;
            switch(dir) {
                case 0:
                    ii = 2;
                    ij = 1;
                    ik = 0;

                    lsh[0] = cctk_lsh[2];
                    lsh[1] = cctk_lsh[1];
                    lsh[2] = cctk_lsh[0];

                    break;
                case 1:
                    ii = 1;
                    ij = 2;
                    ik = 0;

                    lsh[0] = cctk_lsh[2];
                    lsh[1] = cctk_lsh[0];
                    lsh[2] = cctk_lsh[1];

                    break;
                case 2:
                    ii = 1;
                    ij = 0;
                    ik = 2;

                    lsh[0] = cctk_lsh[1];
                    lsh[1] = cctk_lsh[0];
                    lsh[2] = cctk_lsh[2];

                    break;
            }

            // Indices aliases
            int & i = index[ii];
            int & j = index[ij];
            int & k = index[ik];

            // Actual indices
            int __i, __j, __k;

#pragma omp for collapse(2)
            for (__i = THC_M1_NGHOST; __i < lsh[0] - THC_M1_NGHOST; ++__i)
            for (__j = THC_M1_NGHOST; __j < lsh[1] - THC_M1_NGHOST; ++__j) {
                // ----------------------------------------------
                for (__k = THC_M1_NGHOST-1; __k < lsh[2]-THC_M1_NGHOST; ++__k) {
                    index[0] = __i;
                    index[1] = __j;
                    index[2] = __k;
                    int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
										index[2] += 1;
                    int const idxp1 = CCTK_GFINDEX3D(cctkGH, i, j, k);
										index[2] += 1;
                    int const idxp2 = CCTK_GFINDEX3D(cctkGH, i, j, k);
										index[2] -= 3;
                    int const idxm1 = CCTK_GFINDEX3D(cctkGH, i, j, k);
										index[2] += 1;

                    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                        for (int iv = 0; iv < 5; ++iv) {
                            CCTK_REAL const ujm = cons[PINDEX1D(ig, iv)][idxm1];
                            CCTK_REAL const uj = cons[PINDEX1D(ig, iv)][ijk];
                            CCTK_REAL const ujp = cons[PINDEX1D(ig, iv)][idxp1];
                            CCTK_REAL const ujpp = cons[PINDEX1D(ig, iv)][idxp2];

                            CCTK_REAL const dup = ujpp - ujp;
                            CCTK_REAL const duc = ujp - uj;
                            CCTK_REAL const dum = uj - ujm;

                            if (dup*duc < 0 && dum*duc < 0) {
																cons[PINDEX1D(ig, iv)][ijk] = (1.0/2.0) * (uj + (1.0/2.0) * (ujm + ujp));
                            }
                        }
                    }
                }
            }
        }
    }

    delete[] cons;
}
