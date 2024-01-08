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


#ifndef HRSCC_LAX_FRIEDRICHS_RS_HH
#define HRSCC_LAX_FRIEDRICHS_RS_HH

#include <algorithm>
#include <cmath>

#include <cctk.h>

#include <hrscc_config_par.hh>

namespace hrscc {

//! Lax-Friedrichs Riemann Solver
/*!
 *  \tparam claw_t conservation law
 *  \tparam local use the local (as opposed to global) Lax-Friedrichs method
 */
template<typename claw_t, bool local>
class LaxFriedrichsRS {
    public:
        //! conservation law type
        typedef claw_t claw;
        //! number of equations
        enum {nequations = claw::nequations};
        //! if true the method requires the eigenvalues at the interface
        static bool const requires_eigenvalues = local;
        //! if true the method requires the eigenvectors at the interface
        static bool const requires_eigenvectors = false;

        //! solves the given local Riemann problem
        void compute_fluxes(
                    //! [in] Maximum possible propagation speed
                    CCTK_REAL const cbound,
                    //! [in] conserved variables at the interface point
                    //! should be stored as a 2 x claw::nequations array
                    CCTK_REAL const u[2*nequations],
                    //! [in] primitive variables at the interface point
                    //! should be stored as a 2 x claw::nequations array
                    CCTK_REAL const *, // p
                    //! [in] fluxes at the interface point
                    //! should be stored as a 2 x claw::nequations array
                    CCTK_REAL const f[2*nequations],
                    //! [in] eigenvalues at the location of the interface.
                    //! should be stored as a 2 x claw::nequations array
                    CCTK_REAL const eigenvalue[2*nequations],
                    //! [in] left_eigenvectors: not used
                    CCTK_REAL const *, // left_eigenvector,
                    //! [in] right_eigenvectors: not used
                    CCTK_REAL const *, // right_eigenvector,
                    //! [out] computed fluxes
                    CCTK_REAL flux[nequations]) {
            CCTK_REAL const * ul = &u[0];
            CCTK_REAL const * ur = &u[nequations];
            CCTK_REAL const * fl = &f[0];
            CCTK_REAL const * fr = &f[nequations];

            CCTK_REAL s = cbound;
            if(local) {
                s = std::abs(eigenvalue[0]);
                for(int i = 1; i < 2*nequations; ++i) {
                    s = std::max(s, std::abs(eigenvalue[i]));
                }
                // Add dissipation close to the weakly hyperbolic limit
                if(s < config::param::speed_eps*cbound) {
                    s = cbound;
                }
            }

            for(int c = 0; c < nequations; ++c) {
                flux[c] = 0.5*((fr[c] + fl[c]) - s*(ur[c] - ul[c]));
#ifdef HRSCC_CHECK_FOR_NANS
                {
                  using namespace std;
                  assert(!isnan(flux[c]) && !isinf(flux[c]));
                }
#endif
            }
        }
};

} // namespace

#endif
