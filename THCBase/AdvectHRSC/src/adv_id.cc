//  AdvectHRSC: solves the advection equation using HRSCCore and Cactus
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


#include <cmath>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "adv_id.hh"

#include "utils_macro.h"

#define SQ(X) ((X)*(X))

extern "C" CCTK_REAL adv_id_sine(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z) {
    DECLARE_CCTK_PARAMETERS

    return id_phi_amplitude * std::sin(2*M_PI*(
                id_k[0]*x +
                id_k[1]*y +
                id_k[2]*z));
}

extern "C" CCTK_REAL adv_id_gaussian(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z) {
    DECLARE_CCTK_PARAMETERS

    return id_phi_amplitude * std::exp(-0.5*(
        id_isigma[0]*SQ(x - id_gpos[0]) +
        id_isigma[1]*SQ(y - id_gpos[1]) +
        id_isigma[2]*SQ(z - id_gpos[2])));
}

extern "C" CCTK_REAL adv_id_square(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z) {
    DECLARE_CCTK_PARAMETERS

    return id_phi_amplitude * static_cast<CCTK_REAL>(
        (x > id_lowbound[0]) && (y > id_lowbound[1]) &&
        (z > id_lowbound[2]) && (x < id_upbound[0] ) &&
        (y < id_upbound[1] ) && (z < id_upbound[2] ));
}

extern "C" CCTK_REAL adv_id_random(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z) {
    DECLARE_CCTK_PARAMETERS

    return id_phi_amplitude * UTILS_RANDOM(1.0, 1.0);
}

extern "C" void AdvectHRSC_ID(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("AdvectHRSC_ID");
    }

    int const ncells = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

    CCTK_REAL * velx = &vel[0*ncells];
    CCTK_REAL * vely = &vel[1*ncells];
    CCTK_REAL * velz = &vel[2*ncells];

    CCTK_REAL (*func)(cGH const * const , CCTK_REAL,
            CCTK_REAL, CCTK_REAL) = NULL;

    if(CCTK_Equals(id_phi_type, "sine")) {
        func = &adv_id_sine;
    }
    else if(CCTK_Equals(id_phi_type, "gaussian")) {
        func = &adv_id_gaussian;
    }
    else if(CCTK_Equals(id_phi_type, "square")) {
        func = &adv_id_square;
    }
    else if(CCTK_Equals(id_phi_type, "random")) {
        func = &adv_id_random;
    }
    else {
        CCTK_WARN(CCTK_WARN_ABORT, "Something is wrong here! "
                "Dump a 216 digit number and die...");
    }

#pragma omp parallel
    {
        int ijk;
        UTILS_LOOP3(initial_data,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            phi[ijk] = func(cctkGH, x[ijk], y[ijk], z[ijk]);

            if(CCTK_Equals(id_vel_type, "given")) {
                velx[ijk] = id_vel[0];
                vely[ijk] = id_vel[1];
                velz[ijk] = id_vel[2];
            }
            else if(CCTK_Equals(id_vel_type, "rigid_rotation")) {
                velx[ijk] = -y[ijk];
                vely[ijk] = x[ijk];
                velz[ijk] = 0;
            }
            else {
#pragma omp critical
                CCTK_WARN(CCTK_WARN_ABORT, "Something is wrong here! "
                        "Dump a 216 digit number and die...");
            }
        } UTILS_ENDLOOP3(initial_data);
    }
}
