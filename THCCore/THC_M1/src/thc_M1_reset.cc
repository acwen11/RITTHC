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


#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

extern "C" void THC_M1_Reset(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_Reset");
    }

    size_t siz = UTILS_GFSIZE(cctkGH)*nspecies*sizeof(CCTK_REAL);
    size_t isiz = UTILS_GFSIZE(cctkGH)*sizeof(CCTK_INT);

		for (int kk=0; kk<cctk_lsh[2]; kk++) {
			for (int jj=0; jj<cctk_lsh[1]; jj++) {
				for (int ii=0; ii<cctk_lsh[0]; ii++) {
					for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4D = CCTK_VectGFIndex3D(cctkGH, ii, jj, kk, ig);
						rE[i4D] = rad_E_floor;
						rN[i4D] = rad_N_floor;
					}
				}
			}
		}

    // std::memset(rN, rad_N_floor, siz);
    // std::memset(rE, rad_E_floor, siz);
    std::memset(rFx, 0, siz);
    std::memset(rFy, 0, siz);
    std::memset(rFz, 0, siz);
    std::memset(thc_m1_mask, 0, isiz);

		*max_rhserror = 0;
}
