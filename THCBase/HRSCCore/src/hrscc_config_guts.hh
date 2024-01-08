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


#ifndef HRSCC_CONFIG_ITS_ME
#error "This file should not be included directly by the end user!"
#endif

int const metid = config::get_metid();

#ifndef HRSCC_DISABLE_FD
if(config::param::method_i == config::method::FD) {
    switch(metid) {
#define HRSCC_CONFIG_GUTS_FD_CASE(z, N, unused)                                \
        case N:                                                                \
            typename config::select_fd<claw_t, N>::                            \
                method(cctkGH).HRSCC_CONFIG_ACTION;                            \
            break;

#ifdef HRSCC_FD_ENABLE_ONLY
        HRSCC_CONFIG_GUTS_FD_CASE(z, HRSCC_FD_ENABLE_ONLY, z);
#else
        BOOST_PP_REPEAT_FROM_TO(0, HRSCC_CONFIG_NUMBER_FD_METHODS,
                HRSCC_CONFIG_GUTS_FD_CASE, ~);
#endif

#undef HRSCC_CONFIG_GUTS_FD_CASE

        default:
            CCTK_ERROR("Unknown FD method");
    }
}
#endif // HRSCC_DISABLE_FD

#ifndef HRSCC_DISABLE_FV
#ifndef HRSCC_DISABLE_FD
else if(config::param::method_i == config::method::FV) {
#else
if(config::param::method_i == config::method::FV) {
#endif // HRSCC_DISABLE_FD
    switch(metid) {
#define HRSCC_CONFIG_GUTS_FV_CASE(z, N, unused)                                \
        case N:                                                                \
            typename config::select_fv<claw_t, N>::                            \
                method(cctkGH).HRSCC_CONFIG_ACTION;                            \
            break;

#ifdef HRSCC_FV_ENABLE_ONLY
        HRSCC_CONFIG_GUTS_FV_CASE(z, HRSCC_FV_ENABLE_ONLY, z);
#else
        BOOST_PP_REPEAT_FROM_TO(0, HRSCC_CONFIG_NUMBER_FV_METHODS,
                HRSCC_CONFIG_GUTS_FV_CASE, ~);
#endif

#undef HRSCC_CONFIG_GUTS_FV_CASE

        default:
            CCTK_ERROR("Unknown FV method");
    }
}
#endif // HRSCC_DISABLE_FV

else {
    CCTK_ERROR("Unkown numerical method");
}
