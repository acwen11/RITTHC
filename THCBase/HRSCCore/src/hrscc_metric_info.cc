//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2014, David Radice <dradice@caltech.edu>
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


#include <cctk.h>
#include <cctk_Arguments.h>

#include <hrscc_metric_info.hh>

int hrscc::metric_info::idx_lapse = -1;
int hrscc::metric_info::idx_metric[3] = {-1, -1, -1};

namespace {

int varindex(char const * name) {
    int i = CCTK_VarIndex(name);
    if(i < 0) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not find the index of \"%s\"", name);
    }
    return i;
}

void get_metric_indices(cGH const * cctkGH) {
    hrscc::metric_info::idx_lapse     = varindex("ADMBase::alp");
    hrscc::metric_info::idx_metric[0] = varindex("ADMBase::gxx");
    hrscc::metric_info::idx_metric[1] = varindex("ADMBase::gyy");
    hrscc::metric_info::idx_metric[2] = varindex("ADMBase::gzz");
}

}

extern "C" void HRSCC_GetMetric(CCTK_ARGUMENTS) {
    get_metric_indices(cctkGH);
}
