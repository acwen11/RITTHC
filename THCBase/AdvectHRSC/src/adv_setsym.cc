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


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

extern "C" void AdvectHRSC_SetSym(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("AdvectHRSC_SetSym");
    }

    int const sym_phi[3]   = { 1,  1,  1};
    int const sym_velx[3]  = {-1,  1,  1};
    int const sym_vely[3]  = { 1, -1,  1};
    int const sym_velz[3]  = { 1,  1, -1};

    SetCartSymVN(cctkGH, sym_phi,   "AdvectHRSC::phi");
    SetCartSymVN(cctkGH, sym_velx,  "AdvectHRSC::vel[0]");
    SetCartSymVN(cctkGH, sym_vely,  "AdvectHRSC::vel[1]");
    SetCartSymVN(cctkGH, sym_velz,  "AdvectHRSC::vel[2]");
}
