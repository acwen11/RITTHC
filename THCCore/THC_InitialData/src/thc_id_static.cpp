//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
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


#include <math.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "IllinoisGRMHD_headers.h"

#include "utils_macro.h"

using namespace std;
// Let's just redefine this IGM func
void compute_P_cold(igm_eos_parameters eos, CCTK_REAL rho_in,
                              CCTK_REAL &P_cold) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  if(rho_in==0) {
    P_cold   = 0.0;
    return;
  }
  int polytropic_index      = find_polytropic_K_and_Gamma_index(eos,rho_in);
  CCTK_REAL K_ppoly_tab     = eos.K_ppoly_tab[polytropic_index];
  CCTK_REAL Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];
  CCTK_REAL eps_integ_const = eos.eps_integ_const[polytropic_index];

	P_cold = K_ppoly_tab*std::pow(rho_in,Gamma_ppoly_tab);
}
  
extern "C" void THC_ID_Static(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

		igm_eos_parameters eos;
  	initialize_igm_eos_parameters_from_input(igm_eos_key,cctk_time,eos);

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_Static");
    }

    bool const set_Y_e = CCTK_Equals(initial_Y_e, "THC_Initial");

    int const siz = UTILS_GFSIZE(cctkGH);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_static,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            rho[ijk]  = static_rho;
            velx[ijk] = static_velx;
            vely[ijk] = static_vely;
            velz[ijk] = static_velz;
            eps[ijk]  = static_eps;

						// Compute Pcold
						CCTK_REAL static_press;
						compute_P_cold(eos, static_rho, static_press);
						press[ijk] = static_press;

            if (set_Y_e) {
                Y_e[ijk] = static_ye;
            }

            for(int e = 0; e < ntracers; ++e) {
                tracer[ijk + e*siz] = 0;
            }
        } UTILS_ENDLOOP3(thc_id_static);
    }
}
