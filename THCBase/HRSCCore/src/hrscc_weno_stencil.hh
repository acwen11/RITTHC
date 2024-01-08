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


#ifndef HRSCC_WENO_STENCIL_HH
#define HRSCC_WENO_STENCIL_HH

#include <cctk.h>

#include <hrscc_typedefs.hh>

namespace hrscc {

//! Stencil for WENO reconstruction
/*!
 *  This class stores the weights needed for the WENO reconstruction
 *  \tparam eno_width width of the ENO stencil
 *  \tparam weno_type type of the WENO method: symmetric or upwind biased
 *  \tparam weno_optim type of WENO reconstruction: order or bandwith optimized
 */
template<
    int eno_width,
    policy::weno_stencil_t weno_type,
    policy::weno_optim_t weno_optim
>
class WENOStencil {
    public:
        //! The width of the stencil (number of ENO stencils)
        enum {width = eno_width + weno_type};
        //! The weights: \f$ C_{k} \f$.
        static CCTK_REAL const C[eno_width + weno_type];
};

} // namespace

#endif
