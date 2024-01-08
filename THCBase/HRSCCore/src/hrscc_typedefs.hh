//  HRSCCore: HRSC methods for Cactus
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


#ifndef HRSCC_TYPEDEFS_HH
#define HRSCC_TYPEDEFS_HH

namespace hrscc {
//! contains enumerators used for policy specification
namespace policy {

//! working direction in the dimensional split
enum direction_t {x = 0, y = 1, z = 2};

//! reconstruct \f$ v_{i+1/2}^- \f$ or \f$ v_{i+1/2}^+ \f$
enum orientation_t {minus = -1, plus = 1};

//! type of WENO reconstruction: standard or symmetric
enum weno_stencil_t {standard = 0, symmetric = 1};
//! optimize reconstruction for order or bandwith
enum weno_optim_t {order, bandwidth};

//! type of limiting for the smoothness indicators in the WENO method
enum weno_limiter_t {
    //! Jiang and Shu original limiter
    absolute,
    //! no limiting (default behaviour)
    dummy,
    //! Taylor et al. J. Comput. Phys 223 (2007) 384-397 limiter
    relative,
    //! set IS == 0 (this makes the method unstable!)
    unstable
};

} // namespace policy
} // namespace hrscc

#endif
