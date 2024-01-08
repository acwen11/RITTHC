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


#ifndef HRSCC_GRID_HH
#define HRSCC_GRID_HH

namespace hrscc {

//! contains grid indices for coordinates and weights
class gridinfo {
    public:
        //! Cactus grid indices for Grid::x, Grid::y, Grid::z and Grid::r
        static int idx_coordinates[4];
        //! Cactus grid index for CarpetReduce::weight
        static int idx_weight;
};

} // namespace hrscc

#endif
