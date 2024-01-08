//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2013, David Radice <dradice@caltech.edu>
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

#include <hrscc_gridinfo.hh>

int hrscc::gridinfo::idx_coordinates[4] = {-1, -1, -1, -1};
int hrscc::gridinfo::idx_weight = -1;

namespace {

int varindex(char const * name) {
    int i = CCTK_VarIndex(name);
    if(i < 0) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not find the index of \"%s\"", name);
    }
    return i;
}

void get_grid_indices(cGH const * cctkGH) {
    hrscc::gridinfo::idx_coordinates[0] = varindex("Grid::x");
    hrscc::gridinfo::idx_coordinates[1] = varindex("Grid::y");
    hrscc::gridinfo::idx_coordinates[2] = varindex("Grid::z");
    hrscc::gridinfo::idx_coordinates[3] = varindex("Grid::r");
    hrscc::gridinfo::idx_weight         = varindex("CarpetReduce::weight");
}

}

extern "C" void HRSCC_GetGridCoordinates(CCTK_ARGUMENTS) {
    get_grid_indices(cctkGH);
}
