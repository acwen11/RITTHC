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


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "finite_difference.h"
#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"
#include "utils.hh"

using namespace thc::m1;
using namespace utils;

void smooth_gf(CCTK_REAL* cons_gf, CCTK_REAL* cons_tmp, CCTK_INT* m1_mask,
							const cGH *restrict cctkGH, const int *restrict cctk_lsh) {

    DECLARE_CCTK_PARAMETERS
#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_smoothing,
                k, THC_M1_NGHOST, cctk_lsh[2]-THC_M1_NGHOST,
                j, THC_M1_NGHOST, cctk_lsh[1]-THC_M1_NGHOST,
                i, THC_M1_NGHOST, cctk_lsh[0]-THC_M1_NGHOST) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            // Skip masked points
            if (m1_mask[ijk]) {
                continue;
            }

            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

								// Check x first
								int const im1 = CCTK_VectGFIndex3D(cctkGH, i-1, j, k, ig);
								int const ip1 = CCTK_VectGFIndex3D(cctkGH, i+1, j, k, ig);
								int const ip2 = CCTK_VectGFIndex3D(cctkGH, i+2, j, k, ig);

								const CCTK_REAL uijk = cons_gf[ijk];
								const CCTK_REAL uim1 = cons_gf[im1];
								const CCTK_REAL uip1 = cons_gf[ip1];
								const CCTK_REAL uip2 = cons_gf[ip2];

								const CCTK_REAL dim = uijk - uim1;
								const CCTK_REAL dip = uip1 - uijk;
								const CCTK_REAL dip1 = uip2 - uip1;

								const bool checkerx = (dim * dip < 0) && (dip * dip1 < 0);

								if (!checkerx) {
										cons_tmp[i4D] = uijk;
								}
								else {
										int const jp1 = CCTK_VectGFIndex3D(cctkGH, i, j+1, k, ig);
										int const jm1 = CCTK_VectGFIndex3D(cctkGH, i, j-1, k, ig);
										int const kp1 = CCTK_VectGFIndex3D(cctkGH, i, j, k+1, ig);
										int const km1 = CCTK_VectGFIndex3D(cctkGH, i, j, k-1, ig);

										const CCTK_REAL ujp1 = cons_gf[jp1];
										const CCTK_REAL ujm1 = cons_gf[jm1];
										const CCTK_REAL ukp1 = cons_gf[kp1];
										const CCTK_REAL ukm1 = cons_gf[km1];

										const CCTK_REAL djp = ujp1 - uijk;
										const CCTK_REAL dkp = ukp1 - uijk;
										const CCTK_REAL djm = uijk - ujm1;
										const CCTK_REAL dkm = uijk - ukm1;

										const bool checkeryz = (djp*djm < 0) && (dkp*dkm < 0);

										if (checkerboard) {
												cons_tmp[i4D] = (1.0/2.0) * (uijk + (1.0/6.0) * (uip1 + uim1 + ujp1 + ujm1 + ukp1 + ukm1));
										}
										else {
												cons_tmp[i4D] = uijk;
										}
								}
						}

        } UTILS_ENDLOOP3(thc_m1_smoothing);
    } // pragma omp parallel
}

void copy_tmp(CCTK_REAL* cons_gf, CCTK_REAL* cons_tmp, CCTK_INT* m1_mask,
							const cGH *restrict cctkGH, const int *restrict cctk_lsh) {

    DECLARE_CCTK_PARAMETERS
#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_copy_tmp,
                k, THC_M1_NGHOST, cctk_lsh[2]-THC_M1_NGHOST,
                j, THC_M1_NGHOST, cctk_lsh[1]-THC_M1_NGHOST,
                i, THC_M1_NGHOST, cctk_lsh[0]-THC_M1_NGHOST) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            // Skip masked points
            if (m1_mask[ijk]) {
                continue;
            }

            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

								cons_gf[i4D] = cons_tmp[i4D];
						}

        } UTILS_ENDLOOP3(thc_m1_smoothing);
    } // pragma omp parallel
}

extern "C" void THC_M1_Smoothing(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_Smoothing");
    }

		smooth_gf(rN, rN_tmp, thc_m1_mask, cctkGH, cctk_lsh);

		smooth_gf(rE, rE_tmp, thc_m1_mask, cctkGH, cctk_lsh);

		smooth_gf(rFx, rFx_tmp, thc_m1_mask, cctkGH, cctk_lsh);
		smooth_gf(rFy, rFy_tmp, thc_m1_mask, cctkGH, cctk_lsh);
		smooth_gf(rFz, rFz_tmp, thc_m1_mask, cctkGH, cctk_lsh);

		copy_tmp(rN, rN_tmp, thc_m1_mask, cctkGH, cctk_lsh);

		copy_tmp(rE, rE_tmp, thc_m1_mask, cctkGH, cctk_lsh);

		copy_tmp(rFx, rFx_tmp, thc_m1_mask, cctkGH, cctk_lsh);
		copy_tmp(rFy, rFy_tmp, thc_m1_mask, cctkGH, cctk_lsh);
		copy_tmp(rFz, rFz_tmp, thc_m1_mask, cctkGH, cctk_lsh);
}
