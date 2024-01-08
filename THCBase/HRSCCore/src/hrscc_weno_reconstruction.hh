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


#ifndef HRSCC_WENO_RECONSTRUCTION_HH
#define HRSCC_WENO_RECONSTRUCTION_HH

#include <cctk.h>

#include <utils.hh>

#include <hrscc_config_par.hh>
#include <hrscc_eno_stencil.hh>
#include <hrscc_weno_limiter.hh>
#include <hrscc_weno_stencil.hh>
#include <hrscc_weno_weights.hh>

namespace hrscc {

//! Core class for WENO reconstruction
/*!
 *  \tparam eno_stencil_t  type of stencil for the base ENO method
 *  \tparam weno_stencil_t type of WENO stencil to use on top of the ENO method
 *  \tparam weno_weights_t type of weights to use in the final stage
 *  \tparam weno_limiter_t type of smoothness indicator limiter to use
 */
template<
    typename eno_stencil_t,
    typename weno_limiter_t,
    typename weno_stencil_t,
    typename weno_weights_t
>
class WENOReconstruction {
    public:
        typedef eno_stencil_t eno_stencil;
        typedef weno_stencil_t weno_stencil;
        typedef weno_weights_t weno_weights;
        typedef weno_limiter_t weno_limiter;

        //! width of the stencil used for the basic interpolation
        enum {width = eno_stencil::width};
        //! number of stencils to use for the final interpolation
        enum {nstencil = weno_stencil::width};
        //! number of grid points needed for the reconstruction
        enum {total_width = width + nstencil - 1};
        //! number of data points needed for a generic reconstruction
        /*!
         *  Note that in the case in which only the one-directional
         *  reconstruction is computed with a non-symmetric stencil the number
         *  of points is only total_width.
         *
         *  In the general case a buffer of size \e size should be used.
         */
        enum {size = 2*width};

        //! reconstruct an array in \f$ v_{i+1/2}^\pm \f$
        /*!
         *  The input array should be centered around \e i so that, at the end,
         *  we reconstruct in 1/2.
         *
         *  Please mind that the number of data points to the left and to the
         *  right of 0 will change  depending on the chosen reconstruction and
         *  direction.
         */
        template<policy::orientation_t sign>
        CCTK_REAL reconstruct(
                //! [in] grid spacing
                CCTK_REAL,
                //! [in] data to reconstruct
                CCTK_REAL const v[size]
                ) const {
            CCTK_REAL const * u = &v[(width-1)*sign + (sign == 1)];

            for(int k = 0; k < nstencil; ++k) {
                _M_IS[k] = 0;
                for(int m = 0; m < width-1; ++m) {
                    _M_IS[k] += utils::pow<2>(utils::dot<CCTK_REAL, width>::
                            eval(&weno_weights::d[k][m][0], 1,
                                &u[-k*sign], -sign));
                }
            }

            weno_limiter::eval(&_M_IS[0]);

            CCTK_REAL sum_alpha = 0;
            for(int k = 0; k < nstencil; ++k) {
                _M_alpha[k] = weno_stencil::C[k] / utils::pow<2>(
                        config::param::weno_eps + _M_IS[k]);
                sum_alpha += _M_alpha[k];
            }

            CCTK_REAL value = 0;
            CCTK_REAL omega_k = 0;
            for(int k = 0; k < nstencil; ++k) {
                omega_k = _M_alpha[k] / sum_alpha;

                value += omega_k * utils::dot<CCTK_REAL, width>::eval(
                        &eno_stencil::a[k][0], 1, &u[-k*sign], -sign);
            }

            return value;
        }

    private:
        mutable CCTK_REAL _M_alpha[nstencil];
        mutable CCTK_REAL _M_IS[nstencil];
};

//! Alias for the classical WENO3 reconstruction
typedef WENOReconstruction<
    ENOStencil<2>,
    WENOLimiter<2, 2, policy::dummy>,
    WENOStencil<2, policy::standard, policy::order>,
    WENOWeights<2> > WENO3Reconstruction;

//! Alias for the classical WENO5 reconstruction
typedef WENOReconstruction<
    ENOStencil<3>,
    WENOLimiter<3, 3, policy::dummy>,
    WENOStencil<3, policy::standard, policy::order>,
    WENOWeights<3> > WENO5Reconstruction;

//! Alias for the classical WENO7 reconstruction
typedef WENOReconstruction<
    ENOStencil<4>,
    WENOLimiter<4, 4, policy::dummy>,
    WENOStencil<4, policy::standard, policy::order>,
    WENOWeights<4> > WENO7Reconstruction;

//! Alias for the WENO3B reconstruction
typedef WENOReconstruction<
    ENOStencil<3>,
    WENOLimiter<3, 4, policy::relative>,
    WENOStencil<3, policy::symmetric, policy::bandwidth>,
    WENOWeights<3> > WENO3BReconstruction;

//! Alias for the WENO4B reconstruction
typedef WENOReconstruction<
    ENOStencil<4>,
    WENOLimiter<4, 5, policy::relative>,
    WENOStencil<4, policy::symmetric, policy::bandwidth>,
    WENOWeights<4> > WENO4BReconstruction;

} // namespace

#endif
