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

#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "hrscc_config_par.hh"
#include "hrscc_config.hh"

using namespace hrscc;

int config::param::method_i;
int config::param::pplim_i;
int config::param::refluxing_i;
int config::param::reconstruction_i;
int config::param::riemann_solver_i;
int config::param::flux_split_i;
int config::param::system_split_i;

CCTK_REAL config::param::maxspeed;
bool      config::param::cartesian;
CCTK_REAL config::param::speed_eps;
CCTK_REAL config::param::minmod_theta;
CCTK_REAL config::param::pplim_alpha;
CCTK_REAL config::param::limo3_eps;
CCTK_REAL config::param::limo3_r;
CCTK_REAL config::param::mp5_alpha;
CCTK_REAL config::param::weno_alpha;
CCTK_REAL config::param::weno_eps;

namespace {

void param_store(cGH const * cctkGH) {
    DECLARE_CCTK_PARAMETERS

    std::ostringstream msg;
    if(CCTK_Equals(scheme, "FD")) {
        config::param::method_i = config::method::FD;
        msg << "Using FD. ";
    }
    else if(CCTK_Equals(scheme, "FV")) {
        config::param::method_i = config::method::FV;
        msg << "Using FV. ";
    }
    else {
        CCTK_ERROR("");
    }

    if(pplim) {
        config::param::pplim_i = config::pplim::yes;
    }
    else {
        config::param::pplim_i = config::pplim::no;
    }

    if(refluxing) {
        config::param::refluxing_i = config::refluxing::yes;
    }
    else {
        config::param::refluxing_i = config::refluxing::no;
    }

    if(CCTK_Equals(reconstruction, "LimO3")) {
        config::param::reconstruction_i = config::reconstruction::LimO3;
    }
    else if(CCTK_Equals(reconstruction, "MinMod")) {
        config::param::reconstruction_i = config::reconstruction::MinMod;
    }
    else if(CCTK_Equals(reconstruction, "MP5")) {
        config::param::reconstruction_i = config::reconstruction::MP5;
    }
    else if(CCTK_Equals(reconstruction, "SuperBee")) {
        config::param::reconstruction_i = config::reconstruction::SuperBee;
    }
    else if(CCTK_Equals(reconstruction, "VanLeer")) {
        config::param::reconstruction_i = config::reconstruction::VanLeer;
    }
    else if(CCTK_Equals(reconstruction, "WENO3")) {
        config::param::reconstruction_i = config::reconstruction::WENO3;
    }
    else if(CCTK_Equals(reconstruction, "WENO5")) {
        config::param::reconstruction_i = config::reconstruction::WENO5;
    }
    else if(CCTK_Equals(reconstruction, "WENO7")) {
        config::param::reconstruction_i = config::reconstruction::WENO7;
    }
    else {
        CCTK_ERROR("");
    }

    if(CCTK_Equals(riemann_solver, "GLF")) {
        config::param::riemann_solver_i = config::riemann_solver::GLF;
    }
    else if(CCTK_Equals(riemann_solver, "LLF")) {
        config::param::riemann_solver_i = config::riemann_solver::LLF;
    }
    else if(CCTK_Equals(riemann_solver, "HLLE")) {
        config::param::riemann_solver_i = config::riemann_solver::HLLE;
    }
    else {
        CCTK_ERROR("");
    }

    if(CCTK_Equals(flux_split, "GLF")) {
        config::param::flux_split_i = config::flux_split::GLF;
    }
    else if(CCTK_Equals(flux_split, "LLF")) {
        config::param::flux_split_i = config::flux_split::LLF;
    }
    else if(CCTK_Equals(flux_split, "RF")) {
        config::param::flux_split_i = config::flux_split::RF;
    }
    else {
        CCTK_ERROR("");
    }

    if(CCTK_Equals(system_split, "characteristics")) {
        config::param::system_split_i = config::system_split::characteristics;
    }
    else if(CCTK_Equals(system_split, "components")) {
        config::param::system_split_i = config::system_split::components;
    }
    else {
        CCTK_ERROR("");
    }

    config::param::maxspeed     = maxspeed;
    config::param::cartesian    = cartesian;
    config::param::speed_eps    = speed_eps;
    config::param::pplim_alpha  = pplim_alpha;
    config::param::minmod_theta = minmod_theta;
    config::param::limo3_eps    = limo3_eps;
    config::param::limo3_r      = limo3_r;
    config::param::mp5_alpha    = mp5_alpha;
    config::param::weno_alpha   = weno_alpha;
    config::param::weno_eps     = weno_eps;

    msg << "Method ID = " << config::get_metid() << ".";
    std::string msg_str = msg.str();
    CCTK_INFO(msg_str.c_str());
}

}

extern "C" void HRSCC_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_PARAMETERS

    if(CCTK_Equals(flux_split, "RF")) {
        if(!CCTK_Equals(system_split, "characteristics")) {
            CCTK_PARAMWARN("The Roe flux-split requires the use of the "
                    "characteristics-based system-split!");
        }
    }
    if(refluxing && !CCTK_Equals(scheme, "FV")) {
        CCTK_PARAMWARN("Refluxing is only supported by the FV methods!");
    }
}

extern "C" void HRSCC_ParamStore(CCTK_ARGUMENTS) {
    param_store(cctkGH);
}
