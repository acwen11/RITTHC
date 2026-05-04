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

#include "thc_M1_macro.h"

#include "utils.hh"

// #include "utils_macro.h"

using namespace std;

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)>(Y)?(Y):(X))
#define SQ(X) ((X)*(X))

extern "C" void THC_M1_SetMask(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_SetMask");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_setmask,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            // if (hydro_excision && hydro_excision_mask[ijk]) {
            // Set mask ourselves, using conditions and default parameters from THC_Core
            /* Excise region within the AH */
            if (m1_excision && sf_active[excision_surface]) {
                int const sn = excision_surface;

                /* Adapted from CarpetMask */
                CCTK_REAL const x0 = sf_origin_x[sn];
                CCTK_REAL const y0 = sf_origin_y[sn];
                CCTK_REAL const z0 = sf_origin_z[sn];

                CCTK_REAL const theta0 = sf_origin_theta[sn];
                CCTK_REAL const phi0   = sf_origin_phi  [sn];
                CCTK_REAL const dtheta = sf_delta_theta[sn];
                CCTK_REAL const dphi   = sf_delta_phi  [sn];

                int const ntheta = sf_ntheta[sn];
                int const nphi   = sf_nphi  [sn];

                CCTK_REAL const dx = x[ijk] - x0;
                CCTK_REAL const dy = y[ijk] - y0;
                CCTK_REAL const dz = z[ijk] - z0;
                CCTK_REAL const rad = sqrt(SQ(dx) + SQ(dy) + SQ(dz));

                if(rad < 1.0e-12) {
                    thc_m1_mask[ijk] = 1;
                }
                else {
                    CCTK_REAL theta = acos(MIN(1.0, MAX(-1.0, dz/rad)));
                    assert (! std::isnan(theta));
                    assert (theta >= 0);
                    assert (theta <= M_PI);

                    CCTK_REAL phi =
                        fmod (atan2 (dy, dx) + (2 * M_PI),
                            (2 * M_PI));
                    assert (! std::isnan(phi));
                    assert (phi >= 0);
                    assert (phi < 2 * M_PI);

                    int const a = floor ((theta - theta0) / dtheta + 0.5);
                    assert (a >= 0);
                    assert (a < ntheta);
                    int const b = floor ((phi   - phi0  ) / dphi   + 0.5);
                    assert (b >= 0);
                    assert (b < nphi);

                    CCTK_REAL const dr =
                      sf_radius[a + maxntheta * (b + maxnphi * sn)];
                    if(rad < dr * excision_margin) {
                        thc_m1_mask[ijk] = 1;
                    }
                    else {
                        thc_m1_mask[ijk] = 0;
                    }
								}
            }
            else {
                thc_m1_mask[ijk] = 0;
            }

            if (thc_m1_mask[ijk]) {
                for (int ig = 0; ig < nspecies*ngroups; ++ig) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                    rN[i4D] = 0;
                    rE[i4D] = 0;
                    rFx[i4D] = 0;
                    rFy[i4D] = 0;
                    rFz[i4D] = 0;
                }
						}
        } UTILS_ENDLOOP3(thc_m1_setmask);
    }
}
