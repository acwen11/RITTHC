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


#ifndef HRSCC_FD_RECONSTRUCTION_HH
#define HRSCC_FD_RECONSTRUCTION_HH

#include <cctk.h>

#include <utils.hh>

#include <hrscc_typedefs.hh>

namespace hrscc {

//! Class for the derivatives reconstruction
/*!
 *  \tparam reconstructor_t the used reconstructor
 */
template<typename reconstructor_t>
class FDReconstruction {
    public:
        typedef reconstructor_t reconstructor;

        //! width of the stencil used for the reconstruction
        enum {width = reconstructor::width};

        //! reconstruct the derivative of a given grid function
        /*!
         *  \tparam dir direction of the derivative
         */
        template<policy::direction_t dir, policy::orientation_t sign>
        void diff(
                //! [in] Cactus grid hierarchy
                cGH const * const cctkGH,
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [out] derivative of the target grid function
                CCTK_REAL * diff_grid_function
                ) const {
#pragma omp parallel
            {
                // Indices aliases
                int lsh[3] = {0, 0, 0};
                int index[3] = {0, 0, 0};
                int nghost[3] = {0, 0, 0};
                int ii, ij, ik;
                switch(dir) {
                    case 0:
                        ii = 2;
                        ij = 1;
                        ik = 0;

                        lsh[0] = cctkGH->cctk_lsh[2];
                        lsh[1] = cctkGH->cctk_lsh[1];
                        lsh[2] = cctkGH->cctk_lsh[0];
                        nghost[0] = cctkGH->cctk_nghostzones[2];
                        nghost[1] = cctkGH->cctk_nghostzones[1];
                        nghost[2] = cctkGH->cctk_nghostzones[0];

                        break;
                    case 1:
                        ii = 1;
                        ij = 2;
                        ik = 0;

                        lsh[0] = cctkGH->cctk_lsh[2];
                        lsh[1] = cctkGH->cctk_lsh[0];
                        lsh[2] = cctkGH->cctk_lsh[1];
                        nghost[0] = cctkGH->cctk_nghostzones[2];
                        nghost[1] = cctkGH->cctk_nghostzones[0];
                        nghost[2] = cctkGH->cctk_nghostzones[1];

                        break;
                    case 2:
                        ii = 1;
                        ij = 0;
                        ik = 2;

                        lsh[0] = cctkGH->cctk_lsh[1];
                        lsh[1] = cctkGH->cctk_lsh[0];
                        lsh[2] = cctkGH->cctk_lsh[2];
                        nghost[0] = cctkGH->cctk_nghostzones[1];
                        nghost[1] = cctkGH->cctk_nghostzones[0];
                        nghost[2] = cctkGH->cctk_nghostzones[2];

                        break;
                    default:
                        abort();
                }
                assert(nghost[2] >= width);

                // Indices aliases
                int & i = index[ii];
                int & j = index[ij];
                int & k = index[ik];
                int ijk;

                // Actual indices
                int __i, __j, __k;

                // Some scratch space
                // We need to declare these as volatile to prevent the Intel
                // compiler to wrongly optimise the following loops away
                // Tested with icpc-12.1.5.339 and icpc-12.1.3
                CCTK_REAL volatile * u = NULL;
                CCTK_REAL volatile * Ru = NULL;
                try{
                    u  = new CCTK_REAL[cctkGH->cctk_lsh[dir]];
                    Ru = new CCTK_REAL[cctkGH->cctk_lsh[dir]];
                }
                catch(std::bad_alloc & e) {
#pragma omp critical
                    CCTK_WARN(CCTK_WARN_ABORT, "Out of memory!");
                }

                CCTK_REAL const idelta = cctkGH->cctk_levfac[dir] /
                    cctkGH->cctk_delta_space[dir];
                CCTK_REAL const delta = 1.0/idelta;

                // Reconstructor
                reconstructor R;

#pragma omp for collapse(2)
                for(__i = nghost[0]; __i < lsh[0] - nghost[0]; ++__i)
                for(__j = nghost[1]; __j < lsh[1] - nghost[1]; ++__j) {
                    // 1st pass: store the data in a contiguous area of memory
                    for(__k = 0; __k < lsh[2]; ++__k) {
                        index[0] = __i;
                        index[1] = __j;
                        index[2] = __k;
                        ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

                        u[__k]  = grid_function[ijk];
                        Ru[__k] = 0;
                    }

                    // 2nd pass: reconstruct u
                    for(__k = nghost[2] - 1; __k < lsh[2] - nghost[2]; ++__k) {
                        index[2] = __k;
                        // Cast the volatiliness away
                        CCTK_REAL const * uk =
                            const_cast<CCTK_REAL const *>(&u[__k]);
                        Ru[__k] = R.reconstructor::template
                            reconstruct<sign>(delta, uk);
                    }

                    // 3rd pass: compute the derivative
                    for(__k = nghost[2]; __k < lsh[2] - nghost[2]; ++__k) {
                        index[2] = __k;
                        ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

                        diff_grid_function[ijk] = idelta*(Ru[__k] - Ru[__k-1]);
                    }
                }

                delete[] u;
                delete[] Ru;
            }
        }
};

} // namespace

#endif
