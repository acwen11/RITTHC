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


#ifndef HRSCC_FINITE_DIFFERENCE_HH
#define HRSCC_FINITE_DIFFERENCE_HH

#include <cassert>
#include <cmath>

#include <cctk.h>

#include <finite_difference.h>

#include <hrscc_claw_solver.hh>
#include <hrscc_config_par.hh>
#include <hrscc_macro.hh>
#include <hrscc_observer.hh>
#include <hrscc_traits.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

template<typename law_t, typename system_splitter_t, bool pplimiter>
class FiniteDifference;

//! traits of the FiniteDifference class
template<typename law_t, typename system_splitter_t, bool pplimiter>
class traits<FiniteDifference<law_t, system_splitter_t, pplimiter> > {
        public:
            //! number of ghost zones to use
            enum {nghostzones = system_splitter_t::width};
            //! is this method positivity preserving
            static bool const pospres = pplimiter;
            //! this method does not support refluxing
            static bool const refluxing = false;
};


//! class implementing the Finite Difference HRSC method
/*!
 *  \tparam law_t conservation law that we are solving
 *  \tparam system_splitter_t decomposition technique used for systems of c.laws
 *  \tparam pplimiter use positivity-preserving limiter
 *             (see Hu et al. JCP 242 (2013) 169-180)
 */
template<typename law_t, typename system_splitter_t, bool pplimiter = false>
class FiniteDifference: public CLawSolver<law_t, FiniteDifference<law_t,
        system_splitter_t, pplimiter> > {
    public:
        typedef law_t law;
        typedef CLaw<law> claw;
        typedef system_splitter_t system_splitter;
        typedef typename system_splitter::flux_splitter flux_splitter;
        typedef typename system_splitter::reconstructor reconstructor;
        typedef FiniteDifference<law, system_splitter, pplimiter> method;
        typedef CLawSolver<law, method> clawsolver;

        //! number of ghost zones to use
        enum {nghostzones = system_splitter::width};
        //! is this method positivity preserving
        static bool const pospres = pplimiter;

        //! the constructor
        FiniteDifference(
                //! [in] Cactus grid hierarchy
                cGH const * const cctkGH
                ): clawsolver(cctkGH) {}

        //! computes the derivative of a given grid function
        /*!
         *  \tparam dir direction of the derivative
         *
         *  The derivative is computed with an high-order centered
         *  finite-differencing formula, using the same number of ghost regions
         *  as the scheme used to compute the fluxes.
         */
        template<policy::direction_t dir>
        void diff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [out] derivative of the given grid function
                CCTK_REAL * diff_grid_function
                ) const {
#pragma omp parallel
            {
                int const order = nghostzones*2;
                CCTK_REAL const idelta =
                    clawsolver::_M_cctkGH->cctk_levfac[dir] /
                    clawsolver::_M_cctkGH->cctk_delta_space[dir];
                UTILS_LOOP3(diff,
                        k, clawsolver::_M_cctkGH->cctk_nghostzones[2],
                           clawsolver::_M_cctkGH->cctk_lsh[2] -
                           clawsolver::_M_cctkGH->cctk_nghostzones[2],
                        j, clawsolver::_M_cctkGH->cctk_nghostzones[1],
                           clawsolver::_M_cctkGH->cctk_lsh[1] -
                           clawsolver::_M_cctkGH->cctk_nghostzones[1],
                        i, clawsolver::_M_cctkGH->cctk_nghostzones[0],
                           clawsolver::_M_cctkGH->cctk_lsh[0] -
                           clawsolver::_M_cctkGH->cctk_nghostzones[0]) {
                    int const ijk = CCTK_GFINDEX3D(clawsolver::_M_cctkGH,
                            i, j, k);
                    diff_grid_function[ijk] = cdiff_1(clawsolver::_M_cctkGH,
                            grid_function, i, j, k, dir, order) * idelta;
                } UTILS_ENDLOOP3(diff);
            }
        }
        //! computes the derivative of a given grid function
        /*!
         *  \tparam dir direction of the derivative
         *
         *  The derivative is computed with an high-order centered
         *  finite-differencing formula, using the same number of ghost regions
         *  as the scheme used to compute the fluxes.
         */
        template<policy::direction_t dir>
        CCTK_REAL diff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int const index[3]
                ) const {
            int const order = nghostzones*2;
            CCTK_REAL const idelta =
                clawsolver::_M_cctkGH->cctk_levfac[dir] /
                clawsolver::_M_cctkGH->cctk_delta_space[dir];
            return cdiff_1(clawsolver::_M_cctkGH, grid_function, index[0],
                    index[1], index[2], dir, order)*idelta;
        }

        //! actually implements the computation of the RHS
        void compute_rhs() {
            this->_M_compute_rhs<policy::x>();
            this->_M_compute_rhs<policy::y>();
            this->_M_compute_rhs<policy::z>();
        }
    private:
        template<policy::direction_t dir>
        void _M_compute_rhs() {
#pragma omp parallel
            {
                int index[3];
                int lsh[3];
                int nghost[3];
                bool shift[3];
#ifdef HRSCC_CHECK_EIGENVECTORS
                CCTK_REAL eig_err;
#endif

                // We have to leave as the most internal loop the one on the
                // direction. For this reason we will remap the usual indices
                // i,j,k into different points of index[:].
                int ii, ij, ik;
                switch(dir) {
                    case 0:
                        ii = 2;
                        ij = 1;
                        ik = 0;

                        lsh[0] = clawsolver::_M_cctkGH->cctk_lsh[2];
                        lsh[1] = clawsolver::_M_cctkGH->cctk_lsh[1];
                        lsh[2] = clawsolver::_M_cctkGH->cctk_lsh[0];
                        nghost[0] = clawsolver::_M_cctkGH->cctk_nghostzones[2];
                        nghost[1] = clawsolver::_M_cctkGH->cctk_nghostzones[1];
                        nghost[2] = clawsolver::_M_cctkGH->cctk_nghostzones[0];

                        shift[0] = true;
                        shift[1] = false;
                        shift[2] = false;

                        break;
                    case 1:
                        ii = 1;
                        ij = 2;
                        ik = 0;

                        lsh[0] = clawsolver::_M_cctkGH->cctk_lsh[2];
                        lsh[1] = clawsolver::_M_cctkGH->cctk_lsh[0];
                        lsh[2] = clawsolver::_M_cctkGH->cctk_lsh[1];
                        nghost[0] = clawsolver::_M_cctkGH->cctk_nghostzones[2];
                        nghost[1] = clawsolver::_M_cctkGH->cctk_nghostzones[0];
                        nghost[2] = clawsolver::_M_cctkGH->cctk_nghostzones[1];

                        shift[0] = false;
                        shift[1] = true;
                        shift[2] = false;

                        break;
                    case 2:
                        ii = 1;
                        ij = 0;
                        ik = 2;

                        lsh[0] = clawsolver::_M_cctkGH->cctk_lsh[1];
                        lsh[1] = clawsolver::_M_cctkGH->cctk_lsh[0];
                        lsh[2] = clawsolver::_M_cctkGH->cctk_lsh[2];
                        nghost[0] = clawsolver::_M_cctkGH->cctk_nghostzones[1];
                        nghost[1] = clawsolver::_M_cctkGH->cctk_nghostzones[0];
                        nghost[2] = clawsolver::_M_cctkGH->cctk_nghostzones[2];

                        shift[0] = false;
                        shift[1] = false;
                        shift[2] = true;

                        break;
                }
                assert(nghost[2] >= nghostzones);

                // Indices aliases
                int & i = index[ii];
                int & j = index[ij];
                int & k = index[ik];
                int ijk;

                // Actual indices
                int __i, __j, __k;

                CCTK_REAL const tstep =
                    clawsolver::_M_cctkGH->cctk_delta_time /
                    clawsolver::_M_cctkGH->cctk_timefac;
                CCTK_REAL const idelta =
                    clawsolver::_M_cctkGH->cctk_levfac[dir] /
                    clawsolver::_M_cctkGH->cctk_delta_space[dir];
                CCTK_REAL const delta = 1.0/idelta;
                CCTK_REAL const CFL   = tstep*idelta;
                CCTK_REAL const a2CFL = 2*config::param::pplim_alpha*CFL;

                Observer<claw> observer(
                        clawsolver::_M_cctkGH,
                        clawsolver::_M_coordinates,
                        clawsolver::_M_grid_variable,
                        clawsolver::_M_bitmask);
                system_splitter ssplitter;

                // Some scratch space
                CCTK_REAL * u = NULL;
                CCTK_REAL * f = NULL;
                CCTK_REAL * smooth = NULL;
                CCTK_REAL * eigenvalue = NULL;
                CCTK_REAL * left_eigenvector = NULL;
                CCTK_REAL * right_eigenvector = NULL;
                CCTK_REAL * flux_jm = NULL;
                CCTK_REAL * flux_jp = NULL;
                CCTK_REAL * d_ptr = NULL;

                try {
                    u = new CCTK_REAL[claw::nequations*lsh[2]];
                    f = new CCTK_REAL[claw::nequations*lsh[2]];
                    smooth = new CCTK_REAL[lsh[2]];
                    eigenvalue = NULL;
                    if(system_splitter::requires_all_eigenvalues) {
                        eigenvalue = new CCTK_REAL[claw::nequations*
                            clawsolver::_M_cctkGH->cctk_lsh[dir]];
                    }
                    if(system_splitter::requires_all_eigenvectors) {
                        left_eigenvector = new CCTK_REAL[claw::nequations*
                            claw::nequations*clawsolver::
                                _M_cctkGH->cctk_lsh[dir]];
                        right_eigenvector = new CCTK_REAL[claw::nequations*
                            claw::nequations*clawsolver::
                                _M_cctkGH->cctk_lsh[dir]];
                    }
                    flux_jm = new CCTK_REAL[claw::nequations];
                    flux_jp = new CCTK_REAL[claw::nequations];
                    for(int eq = 0; eq < claw::nequations; ++eq) {
                        flux_jm[eq] = 0;
                        flux_jp[eq] = 0;
                    }
                }
                catch(std::bad_alloc & e) {
#pragma omp critical
                    CCTK_WARN(CCTK_WARN_ABORT, "Out of memory!");
                }

#pragma omp for collapse(2)
                for(__i = nghost[0]; __i < lsh[0] - nghost[0]; ++__i)
                for(__j = nghost[1]; __j < lsh[1] - nghost[1]; ++__j) {
                    // 1st pass: store the data in a contiguous area of memory
                    // Note that in this way we ensure that the data that we
                    // really need is in the same NUMA domain as the core on
                    // which we are running (for OpenMP on ccNUMA).
                    for(__k = 0; __k < lsh[2]; ++__k) {
                        index[0] = __i;
                        index[1] = __j;
                        index[2] = __k;
                        ijk = CCTK_GFINDEX3D(clawsolver::_M_cctkGH, i, j, k);
                        observer.jump_to_location(i, j, k);
                        clawsolver::_M_claw->
                            claw::template fluxes<dir>(observer);
                        for(int v = 0; v < claw::nequations; ++v) {
                            u[index[2] + v*lsh[2]] =
                                clawsolver::_M_conserved[v][ijk];
                            f[index[2] + v*lsh[2]] =
                                observer.flux[dir][v];
                        }

                        if(system_splitter::requires_all_eigenvectors) {
                            observer.jump_to_location(i, j, k);
                            clawsolver::_M_claw->
                                claw::template eig<dir>(observer);
#ifdef HRSCC_CHECK_EIGENVECTORS
                            eig_err = observer.check_eigenvectors();
                            if(eig_err > HRSCC_EPSILON) {
#pragma omp critical
                                CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__,
                                    CCTK_THORNSTRING, "The left and right "
                                    "eigenvectors are not orthonormal! The "
                                    "error is %e", eig_err);
                                //observer.reset_eigenvectors();
                            }
#endif
                            for(int vi = 0; vi < claw::nequations; ++vi) {
                                eigenvalue[index[2] + vi*lsh[2]] =
                                    observer.eigenvalue[vi];
                                for(int vj = 0; vj < claw::nequations;
                                        ++vj) {
                                    left_eigenvector[index[2] + vj*lsh[2] +
                                            vi*lsh[2]*claw::nequations] =
                                        observer.left_eigenvector[vi][vj];
                                    right_eigenvector[index[2] + vj*lsh[2] +
                                            vi*lsh[2]*claw::nequations] =
                                        observer.right_eigenvector[vi][vj];
                                }
                            }
                        }
                        else if(system_splitter::requires_all_eigenvalues) {
                            observer.jump_to_location(i, j, k);
                            clawsolver::_M_claw->
                                claw::template eigenvalues<dir>(observer);
                            for(int v = 0; v < claw::nequations; ++v) {
                                eigenvalue[index[2] + v*lsh[2]] =
                                    observer.eigenvalue[v];
                            }
                        }
                    }

                    // 2nd pass: reconstruct the fluxes
                    for(__k = nghost[2] - 1; __k < lsh[2] - nghost[2]; ++__k) {
                        index[0] = __i;
                        index[1] = __j;
                        index[2] = __k;
                        ijk = CCTK_GFINDEX3D(clawsolver::_M_cctkGH, i, j, k);

                        if(system_splitter::requires_interface_eigenvectors)
                        {
                            observer.linterp(i, j, k, shift[0], shift[1],
                                    shift[2]);
                            clawsolver::_M_claw->prim_to_all(observer);
                            clawsolver::_M_claw->
                                claw::template eig<dir>(observer);
#ifdef HRSCC_CHECK_EIGENVECTORS
                            eig_err = observer.check_eigenvectors();
                            if(eig_err > HRSCC_EPSILON) {
#pragma omp critical
                                CCTK_VWarn(CCTK_WARN_ALERT, __LINE__,
                                    __FILE__, CCTK_THORNSTRING, "The left "
                                    "and right eigenvectors are not "
                                    "orthonormal! The error is %e",eig_err);
                                //observer.reset_eigenvectors();
                            }
#endif
                        }
                        else if(system_splitter::
                                requires_interface_eigenvalues){
                            observer.linterp(i, j, k, shift[0], shift[1],
                                    shift[2]);
                            clawsolver::_M_claw->prim_to_all(observer);
                            clawsolver::_M_claw->
                                claw::template eigenvalues<dir>(observer);
                        }

                        // Maximum propagation speed
                        CCTK_REAL cbound = config::param::maxspeed;
                        if(!config::param::cartesian) {
                            cbound = cbound*clawsolver::_M_lapse[ijk] /
                                std::sqrt(clawsolver::_M_metric[dir][ijk]);
                        }

                        ssplitter.compute_fluxes(
                                delta,
                                lsh[2],
                                cbound,
                                &u[index[2]],
                                &f[index[2]],
                                &observer.eigenvalue[0],
                                &eigenvalue[index[2]],
                                &observer.left_eigenvector[0][0],
                                &observer.right_eigenvector[0][0],
                                &left_eigenvector[index[2]],
                                &right_eigenvector[index[2]],
                                flux_jp);

                        // Positivity preserving limiter
                        if(pplimiter) {
                            // Lax-Friedrichs fluxes
                            CCTK_REAL flux_lf_jp[claw::nequations];
                            for(int v = 0; v < claw::nequations; ++v) {
                                CCTK_REAL const * uv = &u[v*lsh[2]];
                                CCTK_REAL const * fv = &f[v*lsh[2]];

                                CCTK_REAL const u_j  = uv[index[2]];
                                CCTK_REAL const u_jp = uv[index[2]+1];
                                CCTK_REAL const f_j  = fv[index[2]];
                                CCTK_REAL const f_jp = fv[index[2]+1];

                                flux_lf_jp[v] = 0.5*(f_j  + f_jp +
                                        cbound*(u_j - u_jp));
                            }

                            // Hybridization coefficient
                            CCTK_REAL theta = 1.0;
                            for(int v = 0; v < claw::nequations; ++v) {
                                // NaN is used to flag conservative
                                // variables that do not need limiting
                                if(std::isnan(claw::conserved_lbound[v])) {
                                    continue;
                                }

                                CCTK_REAL const * uv = &u[v*lsh[2]];
                                CCTK_REAL const u_j  = uv[index[2]];
                                CCTK_REAL const u_jp = uv[index[2]+1];

                                CCTK_REAL theta_m = 1;
                                if(u_j - a2CFL*flux_jp[v] <
                                        claw::conserved_lbound[v]) {
                                    theta_m = (u_j - a2CFL*flux_lf_jp[v]
                                            - claw::conserved_lbound[v])/
                                            (a2CFL*(flux_jp[v] -
                                                    flux_lf_jp[v]));
                                    theta_m = std::max(0.0, theta_m);

#ifdef HRSCC_DEBUG
                                    assert(u_j>=claw::conserved_lbound[v]);
                                    assert(!std::isnan(theta_m));
                                    assert(!std::isinf(theta_m));
#endif
                                }
                                theta = std::min(theta, theta_m);

                                CCTK_REAL theta_p = 1;
                                if(u_jp + a2CFL*flux_jp[v] <
                                        claw::conserved_lbound[v]) {
                                    theta_p = (u_jp + a2CFL*flux_lf_jp[v] -
                                        claw::conserved_lbound[v])/
                                            (a2CFL*(flux_lf_jp[v] -
                                                    flux_jp[v]));
                                    theta_p = std::max(0.0, theta_p);

#ifdef HRSCC_DEBUG
                                    assert(u_jp>=claw::conserved_lbound[v]);
                                    assert(!std::isnan(theta_p));
                                    assert(!std::isinf(theta_p));
#endif
                                }
                                theta = std::min(theta, theta_p);
                            }

                            // Hybridize with Lax-Friedrichs
                            for(int v = 0; v < claw::nequations; ++v) {
                                flux_jp[v] = theta*flux_jp[v] +
                                    (1 - theta)*flux_lf_jp[v];
                            }
                        }

                        for(int v = 0; v < claw::nequations; ++v) {
                            clawsolver::_M_RHS[v][ijk] +=
                                idelta*(flux_jm[v] - flux_jp[v])*
                                static_cast<CCTK_REAL>(
                                    i >= clawsolver::_M_cctkGH->
                                            cctk_nghostzones[0]
                                 && i <  clawsolver::_M_cctkGH->cctk_lsh[0] -
                                         clawsolver::_M_cctkGH->
                                            cctk_nghostzones[0]
                                 && j >= clawsolver::_M_cctkGH->
                                            cctk_nghostzones[1]
                                 && j <  clawsolver::_M_cctkGH->cctk_lsh[1] -
                                         clawsolver::_M_cctkGH->
                                            cctk_nghostzones[1]
                                 && k >= clawsolver::_M_cctkGH->
                                            cctk_nghostzones[2]
                                 && k <  clawsolver::_M_cctkGH->cctk_lsh[2] -
                                         clawsolver::_M_cctkGH->
                                            cctk_nghostzones[2]);
#ifdef HRSCC_CHECK_FOR_NANS
                            {
                                using namespace std;
                                assert(!isnan(
                                        clawsolver::_M_RHS[v][ijk])
                                    && !isinf(
                                        clawsolver::_M_RHS[v][ijk]));
                            }
#endif
                        }
                        d_ptr = flux_jm;
                        flux_jm = flux_jp;
                        flux_jp = d_ptr;
                    }
                }

                delete[] u;
                delete[] f;
                delete[] smooth;
                delete[] eigenvalue;
                delete[] left_eigenvector;
                delete[] right_eigenvector;
                delete[] flux_jm;
                delete[] flux_jp;
            }
        }
};

} // namespace

#endif
