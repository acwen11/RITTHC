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


#include <cassert>

#include "hrscc_config.hh"

namespace hrscc {
namespace config {

int get_fd_metid() {
    assert(HRSCC_CONFIG_NUMBER_FD_METHODS == pplim::siz*system_split::siz*
            flux_split::siz*reconstruction::siz);
    return param::reconstruction_i +
        reconstruction::siz*(param::flux_split_i +
            flux_split::siz*(param::system_split_i +
                system_split::siz*(param::pplim_i)));
}

int get_fv_metid() {
    assert(HRSCC_CONFIG_NUMBER_FV_METHODS == riemann_solver::siz*
            reconstruction::siz*pplim::siz*refluxing::siz);
    return param::reconstruction_i +
        reconstruction::siz*(param::riemann_solver_i +
            riemann_solver::siz*(param::pplim_i +
                pplim::siz*(param::refluxing_i)));
}

int get_metid() {
    if(param::method_i == method::FD) {
        return get_fd_metid();
    }
    else if(param::method_i == method::FV) {
        return get_fv_metid();
    }
    else {
        CCTK_ERROR("This is a bug");
    }
}

} // namespace config
} // namespace hrscc
