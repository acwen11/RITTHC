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


#ifndef ADV_TRAITS_HH
#define ADV_TRAITS_HH

class LinearAdvection;

namespace hrscc {

template<>
class traits<LinearAdvection> {
    public:
        enum {nequations = 1};
        enum {nexternal = 3};
        enum {nbitmasks = 0};
        static bool const pure = true;
};

} // namespace

#endif
