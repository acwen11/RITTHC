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

extern "C" void AdvectHRSC_MoLRegister(CCTK_ARGUMENTS) {
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("AdvectHRSC_MoLRegister");
    }

    int ierr = 0;
    ierr += MoLRegisterEvolved(CCTK_VarIndex("AdvectHRSC::phi"),
            CCTK_VarIndex("AdvectHRSC::rhs"));
    ierr += MoLRegisterConstrained(CCTK_VarIndex("AdvectHRSC::vel[0]"));
    ierr += MoLRegisterConstrained(CCTK_VarIndex("AdvectHRSC::vel[1]"));
    ierr += MoLRegisterConstrained(CCTK_VarIndex("AdvectHRSC::vel[2]"));

    if(ierr) {
        CCTK_WARN(CCTK_WARN_ABORT, "Could not register with MoL");
    }
}
