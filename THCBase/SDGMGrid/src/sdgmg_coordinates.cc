//  SDGMGrid: SDGM grid for Cactus
//  Copyright (C) 2011, Erik Schnetter and David Radice
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


#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_WarnLevel.h>

#include <hrscc_gni_grid.hh>
#include <hrscc_gll_element.hh>

static CCTK_REAL sqr(CCTK_REAL const x) { return x*x; }

// Setup the coordinate system
extern "C" void SDGMG_SetupCoordinates(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INFO("Setting up the coordinates");

    // Setup local coordinates of collocation points
    switch(sdgm_order) {
#define SDGMG_SETUP_COORDINATES(z, N, unused)                                  \
        case N:                                                                \
            {                                                                  \
                hrscc::GNIGrid<hrscc::GLLElement<N> > grid(cctkGH);            \
                grid.setup_coordinates();                                      \
            }                                                                  \
            break;

        BOOST_PP_REPEAT_FROM_TO(0, BOOST_PP_ADD(HRSCC_GLL_ELEMENT_MAX_ORDER, 1),
                SDGMG_SETUP_COORDINATES, ~);

#undef SDGMG_SETUP_COORDINATES

        default:
            CCTK_WARN(CCTK_WARN_ABORT, "The requested GLL elements are not "
                    "implemented");
    }

    // Transform local to global coordinates
    int const m = MultiPatch_GetMap(cctkGH);

    int const np = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
    std::vector<CCTK_INT> ml(np);
    std::vector<CCTK_REAL> xl(np), yl(np), zl(np);
    for (int k=0; k<cctk_lsh[2]; ++k) {
      for (int j=0; j<cctk_lsh[1]; ++j) {
        for (int i=0; i<cctk_lsh[0]; ++i) {
          int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
          int const indn = i + cctk_lsh[0] * (j + cctk_lsh[1] * k);
          xl[indn] = x[ind3d];
          yl[indn] = y[ind3d];
          zl[indn] = z[ind3d];
          ml[indn] = m;
        }
      }
    }

    CCTK_REAL const* const locals[3] = {&xl[0], &yl[0], &zl[0]};
    std::vector<CCTK_REAL> xg(np), yg(np), zg(np);
    CCTK_REAL* const globals[3] = {&xg[0], &yg[0], &zg[0]};
    MultiPatch_LocalToGlobal(cctkGH, 3, np, &ml[0], locals, (void*)globals,
                             NULL, NULL, NULL, NULL, NULL, NULL);

    for (int k=0; k<cctk_lsh[2]; ++k) {
      for (int j=0; j<cctk_lsh[1]; ++j) {
        for (int i=0; i<cctk_lsh[0]; ++i) {
          int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
          int const indn = i + cctk_lsh[0] * (j + cctk_lsh[1] * k);
          x[ind3d] = xg[indn];
          y[ind3d] = yg[indn];
          z[ind3d] = zg[indn];
          r[ind3d] = sqrt(sqr(x[ind3d]) + sqr(y[ind3d]) + sqr(z[ind3d]));
        }
      }
    }
}

// Setup the weights
extern "C" void SDGMG_SetupWeights(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INFO("Setting up reduction weights");

    switch(sdgm_order) {
#define SDGMG_SETUP_WEIGHTS(z, N, unused)                                      \
        case N:                                                                \
            {                                                                  \
                hrscc::GNIGrid<hrscc::GLLElement<N> > grid(cctkGH);            \
                grid.setup_weights();                                          \
            }                                                                  \
            break;

        BOOST_PP_REPEAT_FROM_TO(0, BOOST_PP_ADD(HRSCC_GLL_ELEMENT_MAX_ORDER, 1),
                SDGMG_SETUP_WEIGHTS, ~);

#undef SDGMG_SETUP_WEIGHTS

        default:
            CCTK_WARN(CCTK_WARN_ABORT, "The requested GLL elements are not "
                    "implemented");
    }
}
