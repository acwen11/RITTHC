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


#include <cassert>
#include <cmath>
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "adv_id.hh"

#include "utils_macro.h"

#define SQ(X) ((X)*(X))

extern "C" void AdvectHRSC_Error(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("AdvectHRSC_Error");
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

    CCTK_REAL domainspec[21];
    int ierr = GetDomainSpecification(3, &domainspec[0], &domainspec[3],
            &domainspec[6], &domainspec[9], &domainspec[12], &domainspec[15],
            &domainspec[18]);
    assert(!ierr);

    CCTK_REAL const xmin = domainspec[0];
    CCTK_REAL const xmax = domainspec[3];

    CCTK_REAL const ymin = domainspec[1];
    CCTK_REAL const ymax = domainspec[4];

    CCTK_REAL const zmin = domainspec[2];
    CCTK_REAL const zmax = domainspec[5];

    CCTK_REAL shape[3] = {xmax - xmin, ymax - ymin, zmax - zmin};

#pragma omp parallel
    {
        int ijk;
        UTILS_LOOP3(initial_data,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            CCTK_REAL xp = x[ijk] - cctk_time*velx[ijk];
            CCTK_REAL yp = y[ijk] - cctk_time*vely[ijk];
            CCTK_REAL zp = z[ijk] - cctk_time*velz[ijk];

            CCTK_REAL tmp;

            if(shape[0] > std::numeric_limits<CCTK_REAL>::epsilon()) {
                while(xp < xmin) {
                    xp += shape[0];
                }
                xp = shape[0]*std::modf((xp - xmin)/shape[0], &tmp) + xmin;
            }
            else {
                xp = xmin;
            }

            if(shape[1] > std::numeric_limits<CCTK_REAL>::epsilon()) {
                while(yp < ymin) {
                    yp += shape[1];
                }
                yp = shape[1]*std::modf((yp - ymin)/shape[1], &tmp) + ymin;
            }
            else {
                yp = ymin;
            }

            if(shape[2] > std::numeric_limits<CCTK_REAL>::epsilon()) {
                while(zp < zmin) {
                    zp += shape[2];
                }
                zp = shape[2]*std::modf((zp - zmin)/shape[2], &tmp) + zmin;
            }
            else {
                zp = zmin;
            }

            error[ijk] = phi[ijk] - func(cctkGH, xp, yp, zp);
        } UTILS_ENDLOOP3(initial_data);
    }
}
