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


#ifndef HRSCC_CONFIG_HH
#define HRSCC_CONFIG_HH

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include <cctk.h>
#include <hrscc.hh>

// Here we enumerate all possible schemes.
// First the basic components of each scheme are defined, then the full schemes
// are enumerated by taking into account all possible combinations of the
// components.
// Finally we use Boost.PP to loop over the range of all the schemes and
// instantiate all of them.

namespace hrscc {
namespace config {

///////////////////////////////////////////////////////////////////////////////
// select_reconstruction
///////////////////////////////////////////////////////////////////////////////
template<int reconstruction_i>
struct select_reconstruction;

template<>
struct select_reconstruction<0> {
    typedef LimO3Reconstruction method;
};

template<>
struct select_reconstruction<1> {
    typedef MinModReconstruction method;
};

template<>
struct select_reconstruction<2> {
    typedef MP5Reconstruction method;
};

template<>
struct select_reconstruction<3> {
    typedef SuperBeeReconstruction method;
};

template<>
struct select_reconstruction<4> {
    typedef VanLeerReconstruction method;
};

template<>
struct select_reconstruction<5> {
    typedef WENO3Reconstruction method;
};

template<>
struct select_reconstruction<6> {
    typedef WENO5Reconstruction method;
};

template<>
struct select_reconstruction<7> {
    typedef WENO7Reconstruction method;
};

///////////////////////////////////////////////////////////////////////////////
// select_flux_split
///////////////////////////////////////////////////////////////////////////////
template<typename reconstructor_t, int flux_split_i>
struct select_flux_split;

template<typename reconstructor_t>
struct select_flux_split<reconstructor_t, 0> {
    typedef LaxFriedrichsFS<reconstructor_t, false> method;
};

template<typename reconstructor_t>
struct select_flux_split<reconstructor_t, 1> {
    typedef LaxFriedrichsFS<reconstructor_t, true> method;
};

template<typename reconstructor_t>
struct select_flux_split<reconstructor_t, 2> {
    typedef LaxFriedrichsFS<reconstructor_t, true> efix;
    typedef RoeFS<reconstructor_t, efix> method;
};

///////////////////////////////////////////////////////////////////////////////
// select_system_split
///////////////////////////////////////////////////////////////////////////////
template<typename claw_t, typename flux_split_t, int system_split_i>
struct select_system_split;

template<typename claw_t, typename flux_split_t>
struct select_system_split<claw_t, flux_split_t, 0> {
    typedef CharacteristicSplit<claw_t, flux_split_t> method;
};

template<typename claw_t, typename flux_split_t>
struct select_system_split<claw_t, flux_split_t, 1> {
    typedef ComponentSplit<claw_t, flux_split_t> method;
};

///////////////////////////////////////////////////////////////////////////////
// select_riemann_solver
///////////////////////////////////////////////////////////////////////////////
template<typename claw_t, int riemann_solver_i>
struct select_riemann_solver;

template<typename claw_t>
struct select_riemann_solver<claw_t, 0> {
    typedef LaxFriedrichsRS<claw_t, false> method;
};

template<typename claw_t>
struct select_riemann_solver<claw_t, 1> {
    typedef LaxFriedrichsRS<claw_t, true> method;
};

template<typename claw_t>
struct select_riemann_solver<claw_t, 2> {
    typedef HLLERS<claw_t> method;
};

///////////////////////////////////////////////////////////////////////////////
// select_fd
///////////////////////////////////////////////////////////////////////////////
#define HRSCC_CONFIG_NUMBER_FD_METHODS    96

template<typename claw_t, int metid>
struct select_fd {
    static int const pplim_i = metid/(reconstruction::siz*flux_split::siz*
            system_split::siz);
    static int const system_split_i = (metid - pplim_i*reconstruction::siz*
            flux_split::siz*system_split::siz)/(reconstruction::siz*
            flux_split::siz);
    static int const flux_split_i = (metid - pplim_i*reconstruction::siz*
            flux_split::siz*system_split::siz - system_split_i*
            reconstruction::siz*flux_split::siz)/reconstruction::siz;
    static int const reconstruction_i = metid - pplim_i*reconstruction::siz*
            flux_split::siz*system_split::siz - system_split_i*
            reconstruction::siz*flux_split::siz - flux_split_i*
            reconstruction::siz;

    typedef typename select_reconstruction<reconstruction_i>::method
        reconstruction_t;
    typedef typename select_flux_split<reconstruction_t, flux_split_i>::method
        flux_split_t;
    typedef typename select_system_split<claw_t, flux_split_t, system_split_i>
        ::method system_split_t;
    static bool const pplimiter = pplim_i == 1;

    typedef FiniteDifference<claw_t, system_split_t, pplimiter> method;
};

int get_fd_metid();

///////////////////////////////////////////////////////////////////////////////
// select_fv
///////////////////////////////////////////////////////////////////////////////
#define HRSCC_CONFIG_NUMBER_FV_METHODS    96

template<typename claw_t, int metid>
struct select_fv {
    static int const refluxing_i = metid/(reconstruction::siz*
            riemann_solver::siz*pplim::siz);
    static int const pplim_i = (metid - refluxing_i*reconstruction::siz*
            riemann_solver::siz*pplim::siz)/(reconstruction::siz*
            riemann_solver::siz);
    static int const riemann_solver_i = (metid - refluxing_i*reconstruction::siz*
            riemann_solver::siz*pplim::siz - pplim_i*reconstruction::siz*
            riemann_solver::siz)/reconstruction::siz;
    static int const reconstruction_i = metid - refluxing_i*reconstruction::siz*
            riemann_solver::siz*pplim::siz - pplim_i*reconstruction::siz*
            riemann_solver::siz - riemann_solver_i*reconstruction::siz;

    typedef typename select_reconstruction<reconstruction_i>::method
        reconstruction_t;
    typedef typename select_riemann_solver<claw_t, riemann_solver_i>::method
        riemann_solver_t;
    static bool const pplimiter = pplim_i == 1;
    static bool const refluxing = refluxing_i == 1;

    typedef FiniteVolume<claw_t, reconstruction_t, riemann_solver_t,
            pplimiter, refluxing> method;
};

int get_fv_metid();

int get_metid();

} // namespace config

///////////////////////////////////////////////////////////////////////////////
// Callable functions
///////////////////////////////////////////////////////////////////////////////
#define HRSCC_CONFIG_ITS_ME

template<typename claw_t>
void compute_rhs(cGH const * cctkGH) {
#define HRSCC_CONFIG_ACTION compute_rhs()
#include "hrscc_config_guts.hh"
#undef HRSCC_CONFIG_ACTION
}

template<typename claw_t>
void cons_to_all(cGH const * cctkGH) {
#define HRSCC_CONFIG_ACTION cons_to_all()
#include "hrscc_config_guts.hh"
#undef HRSCC_CONFIG_ACTION
}

template<typename claw_t>
void prim_to_all(cGH const * cctkGH) {
#define HRSCC_CONFIG_ACTION prim_to_all()
#include "hrscc_config_guts.hh"
#undef HRSCC_CONFIG_ACTION
}

#undef HRSCC_CONFIG_ITS_ME

} // namespace hrscc

#endif
