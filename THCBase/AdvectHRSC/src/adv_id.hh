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


#ifndef ADV_ID_HH
#define ADV_ID_HH

extern "C" CCTK_REAL adv_id_sine(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z);

extern "C" CCTK_REAL adv_id_gaussian(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z);

extern "C" CCTK_REAL adv_id_square(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z);

extern "C" CCTK_REAL adv_id_random(cGH const * const cctkGH, CCTK_REAL x,
        CCTK_REAL y, CCTK_REAL z);

#endif
