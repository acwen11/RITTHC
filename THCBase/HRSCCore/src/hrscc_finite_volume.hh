//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2013, David Radice <david.radice@aei.mpg.de>
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


#ifndef HRSCC_FINITE_VOLUME_HH
#define HRSCC_FINITE_VOLUME_HH

#include <cmath>

#include <cctk.h>
#include <cctk_WarnLevel.h>

#include <finite_difference.h>

#include <hrscc_claw_solver.hh>
#include <hrscc_config_par.hh>
#include <hrscc_macro.hh>
#include <hrscc_observer.hh>
#include <hrscc_traits.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

template<typename law_t, typename reconstructor_t,
    typename riemann_solver_t, bool pplimiter, bool refluxing>
class FiniteVolume;


//! traits of the FiniteVolume class
template<typename law_t, typename reconstructor_t,
    typename riemann_solver_t, bool pplimiter, bool refluxing_>
class traits<FiniteVolume<law_t, reconstructor_t, riemann_solver_t,
    pplimiter, refluxing_> > {
        public:
            //! number of ghost zones to use
            enum {nghostzones = reconstructor_t::width};
            //! is this method positivity preserving
            static bool const pospres = pplimiter;
            //! does this method support refluxing
            static bool const refluxing = refluxing_;
};


//! class implementing the 2nd order Finite Volume HRSC method
/*!
 *  \tparam law_t conservation law that we are solving
 *  \tparam reconstructor_t method for primitive variables reconstruction
 *  \tparam riemann_solver_t Riemann solver to use
 *  \tparam pplimiter use positivity-preserving limiter
 *             (see Hu et al. JCP 242 (2013) 169-180)
 *
 *  This class implements a very popular variant of the finite-volume scheme
 *  employing reconstruction in primitive variables and 2nd order accurate flux
 *  evaluations using an approximated Riemann solver
 */
template<typename law_t, typename reconstructor_t,
    typename riemann_solver_t, bool pplimiter = false,
    bool refluxing = false>
class FiniteVolume: public CLawSolver<law_t, FiniteVolume<law_t,
        reconstructor_t, riemann_solver_t, pplimiter, refluxing> > {
    public:
        typedef law_t law;
        typedef CLaw<law> claw;
        typedef reconstructor_t reconstructor;
        typedef riemann_solver_t riemann_solver;
        typedef FiniteVolume<law, reconstructor, riemann_solver,
                pplimiter, refluxing> method;
        typedef CLawSolver<law, method> clawsolver;

        //! number of ghost zones to use
        enum {nghostzones = reconstructor::width};
        //! is this method positivity preserving
        static bool const pospres = pplimiter;

        //! the constructor
        FiniteVolume(
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

                // Observer object
                Observer<claw> observer(
                        clawsolver::_M_cctkGH,
                        clawsolver::_M_coordinates,
                        clawsolver::_M_grid_variable,
                        clawsolver::_M_bitmask);

                // Reconstructor Operator
                reconstructor ROp;

                // Rimenn Solver
                riemann_solver RS;

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

                // Some scratch space
                CCTK_REAL * cons = NULL;
                CCTK_REAL * prim = NULL;
                CCTK_REAL * flux = NULL;

                CCTK_REAL Rprim[2][claw::nequations];
                CCTK_REAL Rcons[2][claw::nequations];
                CCTK_REAL Rflux[2][claw::nequations];

                CCTK_REAL eigenvalue[2][claw::nequations];
                CCTK_REAL left_eigenvector[2][
                    claw::nequations][claw::nequations];
                CCTK_REAL right_eigenvector[2][
                    claw::nequations][claw::nequations];

                // We need flux_jm and flux_jp to be pointers so that
                // we can rotate them
                CCTK_REAL * flux_jm = NULL;
                CCTK_REAL * flux_jp = NULL;
                CCTK_REAL * d_ptr   = NULL;
                try {
                    cons    = new CCTK_REAL[claw::nequations*lsh[2]];
                    prim    = new CCTK_REAL[claw::nequations*lsh[2]];
                    flux    = new CCTK_REAL[claw::nequations*lsh[2]];

                    flux_jm = new CCTK_REAL[claw::nequations];
                    flux_jp = new CCTK_REAL[claw::nequations];
                    for(int c = 0; c < claw::nequations; ++c) {
                        flux_jm[c] = 0;
                        flux_jp[c] = 0;
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
                            cons[index[2] + v*lsh[2]] = observer.conserved[v];
                            prim[index[2] + v*lsh[2]] = observer.primitive[v];
                            flux[index[2] + v*lsh[2]] = observer.flux[dir][v];
                        }
                    }

                    // 2nd pass: compute the fluxes
                    for(__k = nghost[2] - 1; __k < lsh[2] - nghost[2]; ++__k) {
                        index[0] = __i;
                        index[1] = __j;
                        index[2] = __k;
                        ijk = CCTK_GFINDEX3D(clawsolver::_M_cctkGH, i, j, k);

                        // 2.1 compute left state
                        observer.linterp(i, j, k, shift[0], shift[1],
                                shift[2]);
                        for(int c = 0; c < claw::nequations; ++c) {
                            Rprim[0][c] = ROp.reconstructor::template
                                reconstruct<policy::minus>(delta,
                                        &prim[index[2] + c*lsh[2]]);
                            observer.primitive[c] = Rprim[0][c];
                        }
                        clawsolver::_M_claw->prim_to_all(observer);
                        clawsolver::_M_claw->
                            claw::template fluxes<dir>(observer);
                        for(int c = 0; c < claw::nequations; ++c) {
                            Rcons[0][c] = observer.conserved[c];
                            Rflux[0][c] = observer.flux[dir][c];
                        }

                        if(riemann_solver::requires_eigenvectors) {
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
                            for(int a = 0; a < claw::nequations; ++a) {
                                eigenvalue[0][a] = observer.eigenvalue[a];
                                for(int b = 0; b < claw::nequations; ++b) {
                                    left_eigenvector[0][a][b] =
                                        observer.left_eigenvector[a][b];
                                    right_eigenvector[0][a][b] =
                                        observer.right_eigenvector[a][b];
                                }
                            }
                        }
                        else if(riemann_solver::requires_eigenvalues) {
                            clawsolver::_M_claw->
                                claw::template eigenvalues<dir>(observer);
                            for(int c = 0; c < claw::nequations; ++c) {
                                eigenvalue[0][c] = observer.eigenvalue[c];
                            }
                        }

                        // 2.2 compute right state
                        //observer.linterp(i, j, k, shift[0], shift[1],
                        //        shift[2]);
                        for(int c = 0; c < claw::nequations; ++c) {
                            Rprim[1][c] = ROp.reconstructor::template
                                reconstruct<policy::plus>(delta,
                                        &prim[index[2] + c*lsh[2]]);
                            observer.primitive[c] = Rprim[1][c];
                        }
                        clawsolver::_M_claw->prim_to_all(observer);
                        clawsolver::_M_claw->
                            claw::template fluxes<dir>(observer);
                        for(int c = 0; c < claw::nequations; ++c) {
                            Rcons[1][c] = observer.conserved[c];
                            Rflux[1][c] = observer.flux[dir][c];
                        }

                        if(riemann_solver::requires_eigenvectors) {
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

                            for(int a = 0; a < claw::nequations; ++a) {
                                eigenvalue[1][a] = observer.eigenvalue[a];
                                for(int b = 0; b < claw::nequations; ++b) {
                                    left_eigenvector[1][a][b] =
                                        observer.left_eigenvector[a][b];
                                    right_eigenvector[1][a][b] =
                                        observer.right_eigenvector[a][b];
                                }
                            }
                        }
                        else if(riemann_solver::requires_eigenvalues) {
                            clawsolver::_M_claw->
                                claw::template eigenvalues<dir>(observer);
                            for(int c = 0; c < claw::nequations; ++c) {
                                eigenvalue[1][c] = observer.eigenvalue[c];
                            }
                        }

                        // 2.3 Maximum propagation speed
                        CCTK_REAL cbound = config::param::maxspeed;
                        if(!config::param::cartesian) {
                            cbound = cbound*clawsolver::_M_lapse[ijk] /
                                std::sqrt(clawsolver::_M_metric[dir][ijk]);
                        }

                        // 2.4 Riemann solver
                        RS.compute_fluxes(
                                cbound,
                                const_cast<CCTK_REAL const *>(&Rcons[0][0]),
                                const_cast<CCTK_REAL const *>(&Rprim[0][0]),
                                const_cast<CCTK_REAL const *>(&Rflux[0][0]),
                                const_cast<CCTK_REAL const *>(
                                    &eigenvalue[0][0]),
                                const_cast<CCTK_REAL const *>(
                                    &left_eigenvector[0][0][0]),
                                const_cast<CCTK_REAL const *>(
                                    &right_eigenvector[0][0][0]),
                                flux_jp);

                        // 3rd step: positivity preserving limiter
                        if(pplimiter) {
                            // Lax-Friedrichs fluxes
                            CCTK_REAL flux_lf_jp[claw::nequations];
                            for(int c = 0; c < claw::nequations; ++c) {
                                // Some helpful alias
                                CCTK_REAL const * consc   = &cons[c*lsh[2]];
                                CCTK_REAL const * fluxc   = &flux[c*lsh[2]];
                                CCTK_REAL const consc_j   = consc[index[2]];
                                CCTK_REAL const consc_jp1 = consc[index[2]+1];
                                CCTK_REAL const fluxc_j   = fluxc[index[2]];
                                CCTK_REAL const fluxc_jp1 = fluxc[index[2]+1];

                                flux_lf_jp[c] = 0.5*(fluxc_j + fluxc_jp1 +
                                        cbound*(consc_j - consc_jp1));
                            }

                            // Hybridization coefficient
                            CCTK_REAL theta = 1.0;
                            for(int c = 0; c < claw::nequations; ++c) {
                                // NaN is used to flag conserved
                                // variables that do not need limiting
                                if(std::isnan(claw::conserved_lbound[c])) {
                                    continue;
                                }

                                // Some helpful alias
                                CCTK_REAL const * consc = &cons[c*lsh[2]];
                                CCTK_REAL const consc_j   = consc[index[2]];
                                CCTK_REAL const consc_jp1 = consc[index[2]+1];

                                // Ensure positivity of the left state
                                CCTK_REAL theta_m = 1;
                                if(consc_j - a2CFL*flux_jp[c] <
                                        claw::conserved_lbound[c]) {
                                    theta_m = (consc_j - a2CFL*flux_lf_jp[c]
                                            - claw::conserved_lbound[c])/
                                            (a2CFL*(flux_jp[c] -
                                                    flux_lf_jp[c]));
                                    theta_m = std::max(0.0, theta_m);
#ifdef HRSCC_DEBUG
                                    assert(consc_j >=
                                            claw::conserved_lbound[c]);
                                    assert(!std::isnan(theta_m));
                                    assert(!std::isinf(theta_m));
#endif
                                }
                                theta = std::min(theta, theta_m);

                                // Ensure positivity of the right state
                                CCTK_REAL theta_p = 1;
                                if(consc_jp1 + a2CFL*flux_jp[c] <
                                        claw::conserved_lbound[c]) {
                                    theta_p = (consc_jp1 + a2CFL*flux_lf_jp[c]
                                            - claw::conserved_lbound[c])/
                                            (a2CFL*(flux_lf_jp[c] -
                                                    flux_jp[c]));
                                    theta_p = std::max(0.0, theta_p);
#ifdef HRSCC_DEBUG
                                    assert(consc_jp1 >=
                                            claw::conserved_lbound[c]);
                                    assert(!std::isnan(theta_p));
                                    assert(!std::isinf(theta_p));
#endif
                                }
                                theta = std::min(theta, theta_p);
                            }

                            // Hybridize with Lax-Friedrichs
                            for(int c = 0; c < claw::nequations; ++c) {
                                flux_jp[c] = theta*flux_jp[c] +
                                    (1 - theta)*flux_lf_jp[c];
                            }
                        }

                        // 4th step: add everything to the RHS
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
                        if(refluxing) {
                            for(int v = 0; v < claw::nequations; ++v) {
                                clawsolver::_M_num_flux[dir][v][ijk] =
                                    flux_jm[v];
                            }
                            if(__k == lsh[2] - nghost[2] - 1) {
                                index[0] = __i;
                                index[1] = __j;
                                index[2] = __k + 1;
                                ijk = CCTK_GFINDEX3D(clawsolver::_M_cctkGH,
                                        i, j, k);
                                for(int v = 0; v < claw::nequations; ++v) {
                                    clawsolver::_M_num_flux[dir][v][ijk] =
                                        flux_jp[v];
                                }
                            }
                        }
                        d_ptr = flux_jm;
                        flux_jm = flux_jp;
                        flux_jp = d_ptr;
                    }
                }

                delete[] flux_jp;
                delete[] flux_jm;
                delete[] flux;
                delete[] prim;
                delete[] cons;
            }
        }
};

} // namespace

#endif
