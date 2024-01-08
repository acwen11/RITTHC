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


#include <cstdio>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "adv_claw.hh"

namespace {

int index(char const * name) {
    int i = CCTK_VarIndex(name);
    if(i < 0) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not find the index of \"%s\"", name);
    }
    return i;
}

} // namespace

namespace hrscc {

template<> int CLaw<LinearAdvection>::conserved_idx[1] = {0};
template<> int CLaw<LinearAdvection>::primitive_idx[1] = {0};
template<> int CLaw<LinearAdvection>::rhs_idx[1] = {0};
template<> int CLaw<LinearAdvection>::field_idx[3] = {0, 0, 0};
template<> int CLaw<LinearAdvection>::bitmask_idx[0] = {};
template<> int CLaw<LinearAdvection>::num_flux_idx[3] = {0, 0, 0};

template<> CCTK_REAL CLaw<LinearAdvection>::conserved_lbound[1] = {0};

} // namespace

using namespace hrscc;

LinearAdvection::LinearAdvection() {
    DECLARE_CCTK_PARAMETERS

    char vname[BUFSIZ];

    CLaw<LinearAdvection>::conserved_idx[0] = index("AdvectHRSC::phi");
    CLaw<LinearAdvection>::primitive_idx[0] = index("AdvectHRSC::phi");

    CLaw<LinearAdvection>::rhs_idx[0]       = index("AdvectHRSC::rhs");

    CLaw<LinearAdvection>::field_idx[0]     = index("AdvectHRSC::vel[0]");
    CLaw<LinearAdvection>::field_idx[1]     = index("AdvectHRSC::vel[1]");
    CLaw<LinearAdvection>::field_idx[2]     = index("AdvectHRSC::vel[2]");

    if (refluxing) {
        for (int d = 0; d < 3; ++d) {
            std::snprintf(vname, BUFSIZ, "Refluxing::flux[%d]", d);
            CLaw<LinearAdvection>::num_flux_idx[d] = index(vname);
        }
    }
}
