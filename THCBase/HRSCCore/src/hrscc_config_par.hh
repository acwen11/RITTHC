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


#ifndef HRSCC_CONFIG_PAR_HH
#define HRSCC_CONFIG_PAR_HH

#include <cctk.h>

// We cannot use Cactus parameters as template parameters, because the final
// objects will be created outside of this thorn, were the parameters are not
// available. For this reason we store all of them in a global data structure.
// All the parameters are represented as int so that we can loop over them
// using Boost.PP to instantiate all the classes.

namespace hrscc {
namespace config {

struct method {
    static int const FD  = 0;
    static int const FV  = 1;
    static int const siz = 2;
};

struct pplim {
    static int const no  = 0;
    static int const yes = 1;
    static int const siz = 2;
};

struct refluxing {
    static int const no  = 0;
    static int const yes = 1;
    static int const siz = 2;
};

struct reconstruction {
    static int const LimO3    = 0;
    static int const MinMod   = 1;
    static int const MP5      = 2;
    static int const SuperBee = 3;
    static int const VanLeer  = 4;
    static int const WENO3    = 5;
    static int const WENO5    = 6;
    static int const WENO7    = 7;
    static int const siz      = 8;
};

struct riemann_solver {
    static int const GLF  = 0;
    static int const LLF  = 1;
    static int const HLLE = 2;
    static int const siz  = 3;
};

struct flux_split {
    static int const GLF = 0;
    static int const LLF = 1;
    static int const RF  = 2;
    static int const siz = 3;
};

struct system_split {
    static int const characteristics = 0;
    static int const components      = 1;
    static int const siz             = 2;
};

// Here is where we store all the parameters
struct param {
    static int method_i;
    static int pplim_i;
    static int refluxing_i;
    static int reconstruction_i;
    static int riemann_solver_i;
    static int flux_split_i;
    static int system_split_i;

    static CCTK_REAL maxspeed;
    static bool      cartesian;
    static CCTK_REAL speed_eps;
    static CCTK_REAL pplim_alpha;
    static CCTK_REAL minmod_theta;
    static CCTK_REAL limo3_eps;
    static CCTK_REAL limo3_r;
    static CCTK_REAL mp5_alpha;
    static CCTK_REAL weno_alpha;
    static CCTK_REAL weno_eps;
};

} // namespace config

} // namespace hrscc

#endif
