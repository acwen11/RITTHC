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


#ifndef HRSCC_SDG_METHOD_HH
#define HRSCC_SDG_METHOD_HH

#include <algorithm>

#include <cctk.h>
#include <cctk_WarnLevel.h>

#include <utils.hh>

#include <hrscc_claw.hh>
#include <hrscc_claw_solver.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

template<typename law_t, typename grid_t, typename riemann_solver_t>
class SDGMethod;

//! traits of the SDGMethod class
template<typename law_t, typename grid_t, typename riemann_solver_t>
class traits<SDGMethod<law_t, grid_t, riemann_solver_t> > {
    public:
        //! number of ghost zones to use
        enum {nghostzones = 1};
            //! this method is not positivity preserving
            static bool const pospres = false;
            //! this method does not support refluxing
            static bool const refluxing = false;
};

//! class implementing the Spectral Discontinuous Galerkin Method
/*!
 *  \tparam law_t conservation law that we are solving
 *  \tparam grid_t numerical grid to use
 *  \tparam riemann_solver_t Riemann solver used to compute the fluxes
 *
 *  Currently only the GNIGrid is supported by this class.
 */
template<typename law_t, typename grid_t, typename riemann_solver_t>
class SDGMethod: public CLawSolver<law_t, SDGMethod<law_t, grid_t,
        riemann_solver_t> >  {
    public:
        typedef law_t law;
        typedef CLaw<law> claw;
        typedef grid_t grid;
        typedef riemann_solver_t riemann_solver;
        typedef SDGMethod<law, grid, riemann_solver> method;
        typedef CLawSolver<law, method> clawsolver;

        //! number of ghost zones to use
        enum {nghostzones = 1};
        //! pplimiter for DG not implemented yet
        static bool const pospres = false;

        //! the constructor
        SDGMethod(
                cGH const * const cctkGH
                ): clawsolver(cctkGH), mesh(cctkGH) {}

        //! computes the pseudo-spectral derivative of a given function
        /*!
         *  \tparam dir direction of the derivative
         *  \see GNIGrid::diff
         */
        template<policy::direction_t dir>
        void diff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [out] derivative of the given grid function
                CCTK_REAL * diff_grid_function
                ) const {
            mesh.template diff<dir>(grid_function,
                    diff_grid_function);
        }
        //! \brief computes the pseudo-spectral derivatives of a given function
        //! in a point
        /*!
         *  \tparam dir direction of the derivative
         *  \see GNIGrid::diff
         */
        template<policy::direction_t dir>
        CCTK_REAL diff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int const i, int const j, int const k
                ) const {
            return mesh.template diff<dir>(grid_function, i, j, k);
        }

        //! \brief computes the weak derivative of a given grid function
        template<policy::direction_t dir>
        CCTK_REAL wdiff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int const i, int const j, int const k
                ) const {
            return mesh.template wdiff<dir>(grid_function, i, j, k);
        }
        //! \brief computes the weighted value of the jump of a given
        //! grid function in a point
        /*!
         *  The jump in a given direction, \f$d\f$, is usually defined as
         *  \f[
                (u_+ \vec{n}_+ + u_- \vec{n}_-) \cdot \vec{d} = u_L - u_R,
            \f]
         *  where \f$u_L\f$ and \f$u_R\f$ are the left and right states
         *  along \f$\vec{d}\f$.
         *
         *  The weighted jump operator returns the integrated jumps divided by
         *  the mass matrix, i.e.
         *  \f[
                \frac{u_L - u_R}{\Delta x\, w},
            \f]
         *  where \f$w\f$ is the quadrature weight of the given grid point
         */
        template<policy::direction_t dir>
        CCTK_REAL jump(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int const i, int const j, int const k
                ) const {
            return mesh.template jump<dir>(grid_function, i, j, k);
        }

        //! actually computes the right-hand-side
        /*!
         *  This function computes the flux part of the RHS, taking care of
         *  multiplying it by the inverse the mass-matrix, and adds it to the
         *  RHS vector.
         *
         *  The user should take care of initializing the RHS with the source
         *  terms. Note that, in the case of a collocation method, no quadrature
         *  is required for the calculation of the sources.
         *
         *  NOTE: currently this method is not very efficient because we have
         *  to compute the flux function at the boundary of the elements and
         *  solve the Riemann problems twice.
         *
         *  TODO: I must cleanup this routine and split it into several parts
         *  as it is now it's a mess
         */
        void compute_rhs() {
            // total number of points in each element
            int const elemsize = grid::element::ndof;
            // total number of points in each element including boundaries
            int const xelemsize = utils::pow<3>(grid::element::npoints + 2);
            // stride for local arrays
            int const stride[3] = {
                grid::element::npoints*grid::element::npoints,
                grid::element::npoints, 1};
            // stride for extended local arrays
            int const xstride[3] = {
                utils::pow<2>(grid::element::npoints + 2),
                grid::element::npoints + 2, 1};

            // shape of the elements
            CCTK_REAL const shape[3] = {mesh.delta_space[0],
                mesh.delta_space[1], mesh.delta_space[2]};

            // inverse mass matrix (assuming that the mass matrix is diagonal)
            CCTK_REAL imass[elemsize];
            for(int lk = 0; lk < grid::element::npoints; ++lk)
            for(int lj = 0; lj < grid::element::npoints; ++lj)
            for(int li = 0; li < grid::element::npoints; ++li) {
                int const idx = li*stride[0] + lj*stride[1] + lk*stride[2];
                imass[idx] = 1.0 / (
                        grid::element::weight[li]*0.5*shape[0]*
                        grid::element::weight[lj]*0.5*shape[1]*
                        grid::element::weight[lk]*0.5*shape[2]);
            }

#pragma omp parallel
            {
                Observer<claw> observer(
                        clawsolver::_M_cctkGH,
                        clawsolver::_M_coordinates,
                        clawsolver::_M_grid_variable,
                        clawsolver::_M_bitmask);

                // Riemann solver
                riemann_solver RS;

                // some scrach space
                CCTK_REAL __scratch[6][claw::nequations][xelemsize];

                int const offset = xstride[0] + xstride[1] + xstride[2];
                CCTK_REAL * scratch[6][claw::nequations];
                for(int i = 0; i < 6; ++i) {
                    for(int c = 0; c < claw::nequations; ++c) {
                        scratch[i][c] = &__scratch[i][c][offset];
                    }
                }

                // storage for the fluxes from the Riemann solver
                CCTK_REAL locflux[claw::nequations];

                // aliases (used to define the Riemann problems)
                CCTK_REAL u[2][claw::nequations];
                CCTK_REAL p[2][claw::nequations];
                CCTK_REAL f[2][claw::nequations];
                CCTK_REAL eigenvalue[2][claw::nequations];
                CCTK_REAL left_eigenvector[2][
                        claw::nequations][claw::nequations];
                CCTK_REAL right_eigenvector[2][
                        claw::nequations][claw::nequations];

#pragma omp for collapse(3)
                for(int ek = 0; ek < int(mesh.nelem[2]); ++ek)
                for(int ej = 0; ej < int(mesh.nelem[1]); ++ej)
                for(int ei = 0; ei < int(mesh.nelem[0]); ++ei) {
                    int const element[3] = {ei, ej, ek};
                    int const iorigin[3] = {
                        mesh.grid::template element_iorigin<policy::x>(element),
                        mesh.grid::template element_iorigin<policy::y>(element),
                        mesh.grid::template element_iorigin<policy::z>(element)
                    };

                    // compute the flux function in the element (and ghost
                    // regions)
                    for(int lk = -1; lk < grid::element::npoints+1; ++lk)
                    for(int lj = -1; lj < grid::element::npoints+1; ++lj)
                    for(int li = -1; li < grid::element::npoints+1; ++li) {
                        int const i = li + iorigin[0];
                        int const j = lj + iorigin[1];
                        int const k = lk + iorigin[2];
                        observer.jump_to_location(i, j, k);

                        clawsolver::_M_claw->claw::template fluxes<policy::x>(
                                observer);
                        clawsolver::_M_claw->claw::template fluxes<policy::y>(
                                observer);
                        clawsolver::_M_claw->claw::template fluxes<policy::z>(
                                observer);
                        for(int dir = 0; dir < 3; ++dir) {
                            for(int c = 0; c < claw::nequations; ++c) {
                                scratch[dir][c][li*xstride[0] + lj*xstride[1] +
                                    lk*xstride[2]] = observer.flux[dir][c];
                            }
                        }
                    }

                    // compute the codifferential of the fluxes
                    for(int c = 0; c < claw::nequations; ++c) {
                        grid::element::template codiff<policy::x>(shape,
                                &scratch[0][c][0], xstride,
                                &scratch[3][c][0], xstride);
                        grid::element::template codiff<policy::y>(shape,
                                &scratch[1][c][0], xstride,
                                &scratch[4][c][0], xstride);
                        grid::element::template codiff<policy::z>(shape,
                                &scratch[2][c][0], xstride,
                                &scratch[5][c][0], xstride);
                    }

                    // TODO: this has to be split out as a function
                    // solve the Riemann problems
#define HRSCC_SDGMETHOD_SOLVERIEMANN(DIRECTION,SIDE)                           \
                    for(int lb = 0; lb < grid::element::npoints; ++lb)         \
                    for(int la = 0; la < grid::element::npoints; ++la) {       \
                        int const li0 = DIRECTION == 0 ?                       \
                            (SIDE)*(grid::element::npoints) - 1 : la;          \
                        int const lj0 = DIRECTION == 1 ?                       \
                            (SIDE)*(grid::element::npoints) - 1 :              \
                            (DIRECTION == 0 ? la : lb);                        \
                        int const lk0 = DIRECTION == 2 ?                       \
                            (SIDE)*(grid::element::npoints) - 1 : lb;          \
                        int const lijk0 = li0*xstride[0] + lj0*xstride[1] +    \
                            lk0*xstride[2];                                    \
                        int const li1 = li0 + static_cast<int>(DIRECTION==0);  \
                        int const lj1 = lj0 + static_cast<int>(DIRECTION==1);  \
                        int const lk1 = lk0 + static_cast<int>(DIRECTION==2);  \
                        int const lijk1 = li1*xstride[0] + lj1*xstride[1] +    \
                            lk1*xstride[2];                                    \
                                                                               \
                        int const i0 = iorigin[0] + li0;                       \
                        int const j0 = iorigin[1] + lj0;                       \
                        int const k0 = iorigin[2] + lk0;                       \
                        int const ijk0 = CCTK_GFINDEX3D(clawsolver::_M_cctkGH, \
                                i0, j0, k0);                                   \
                        int const i1 = iorigin[0] + li1;                       \
                        int const j1 = iorigin[1] + lj1;                       \
                        int const k1 = iorigin[2] + lk1;                       \
                        int const ijk1 = CCTK_GFINDEX3D(clawsolver::_M_cctkGH, \
                                i1, j1, k1);                                   \
                                                                               \
                        CCTK_REAL const shapefact = 0.5*0.5*(DIRECTION == 0 ?  \
                            shape[1]*shape[2] : (DIRECTION == 1 ?              \
                                    shape[0]*shape[2] : shape[0]*shape[1]));   \
                                                                               \
                        for(int c = 0; c < claw::nequations; ++c) {            \
                            u[0][c] = clawsolver::_M_conserved[c][ijk0];       \
                            u[1][c] = clawsolver::_M_conserved[c][ijk1];       \
                            p[0][c] = clawsolver::_M_primitive[c][ijk0];       \
                            p[1][c] = clawsolver::_M_primitive[c][ijk1];       \
                            f[0][c] = scratch[DIRECTION][c][lijk0];            \
                            f[1][c] = scratch[DIRECTION][c][lijk1];            \
                        }                                                      \
                                                                               \
                        if(riemann_solver::requires_eigenvectors) {            \
                            observer.jump_to_location(i0, j0, k0);             \
                            clawsolver::_M_claw->claw::template                \
                                eig<DIRECTION>(observer);                      \
                            for(int i = 0; i < claw::nequations; ++i) {        \
                                eigenvalue[0][i] =                             \
                                        observer.eigenvalue[i];                \
                                for(int j = 0; j < claw::nequations; ++j) {    \
                                    left_eigenvector[0][i][j] =                \
                                            observer.left_eigenvector[i][j];   \
                                    right_eigenvector[0][i][j] =               \
                                            observer.right_eigenvector[i][j];  \
                                }                                              \
                            }                                                  \
                                                                               \
                            observer.jump_to_location(i1, j1, k1);             \
                            clawsolver::_M_claw->claw::template                \
                                eig<DIRECTION>(observer);                      \
                            for(int i = 0; i < claw::nequations; ++i) {        \
                                eigenvalue[1][i] =                             \
                                        observer.eigenvalue[i];                \
                                for(int j = 0; j < claw::nequations; ++j) {    \
                                    left_eigenvector[1][i][j] =                \
                                            observer.left_eigenvector[i][j];   \
                                    right_eigenvector[1][i][j] =               \
                                            observer.right_eigenvector[i][j];  \
                                }                                              \
                            }                                                  \
                        }                                                      \
                        else if(riemann_solver::requires_eigenvalues) {        \
                            observer.jump_to_location(i0, j0, k0);             \
                            clawsolver::_M_claw->claw::template                \
                                eigenvalues<DIRECTION>(observer);              \
                            for(int i = 0; i < claw::nequations; ++i) {        \
                                eigenvalue[0][i] =                             \
                                        observer.eigenvalue[i];                \
                            }                                                  \
                                                                               \
                            observer.jump_to_location(i1, j1, k1);             \
                            clawsolver::_M_claw->claw::template                \
                                eigenvalues<DIRECTION>(observer);              \
                            for(int i = 0; i < claw::nequations; ++i) {        \
                                eigenvalue[1][i] =                             \
                                        observer.eigenvalue[i];                \
                            }                                                  \
                        }                                                      \
                                                                               \
                        RS.compute_fluxes(                                     \
                                config::param::maxspeed,                       \
                                const_cast<CCTK_REAL const *>(&u[0][0]),       \
                                const_cast<CCTK_REAL const *>(&p[0][0]),       \
                                const_cast<CCTK_REAL const *>(&f[0][0]),       \
                                const_cast<CCTK_REAL const *>(                 \
                                    &eigenvalue[0][0]),                        \
                                const_cast<CCTK_REAL const *>(                 \
                                    &left_eigenvector[0][0][0]),               \
                                const_cast<CCTK_REAL const *>(                 \
                                    &right_eigenvector[0][0][0]),              \
                                &locflux[0]);                                  \
                                                                               \
                        /* Notice the non-standard sign convention! */         \
                        for(int c = 0; c < claw::nequations; ++c) {            \
                            if(SIDE == 0) {                                    \
                                scratch[DIRECTION+3][c][lijk1] -=              \
                                        locflux[c] *                           \
                                        shapefact *                            \
                                        grid::element::weight[la] *            \
                                        grid::element::weight[lb];             \
                            }                                                  \
                            else {                                             \
                                scratch[DIRECTION+3][c][lijk0] +=              \
                                        locflux[c] *                           \
                                        shapefact *                            \
                                        grid::element::weight[la] *            \
                                        grid::element::weight[lb];             \
                            }                                                  \
                        }                                                      \
                    }

                    // TODO: solve only one Riemann problem for each interface
                    // quadrature point
                    HRSCC_SDGMETHOD_SOLVERIEMANN(policy::x, 0)
                    HRSCC_SDGMETHOD_SOLVERIEMANN(policy::x, 1)
                    HRSCC_SDGMETHOD_SOLVERIEMANN(policy::y, 0)
                    HRSCC_SDGMETHOD_SOLVERIEMANN(policy::y, 1)
                    HRSCC_SDGMETHOD_SOLVERIEMANN(policy::z, 0)
                    HRSCC_SDGMETHOD_SOLVERIEMANN(policy::z, 1)
#undef HRSCC_SDGMETHOD_SOLVERIEMANN

                    // Add everything to the RHS vector
                    for(int lk = 0; lk < grid::element::npoints; ++lk)
                    for(int lj = 0; lj < grid::element::npoints; ++lj)
                    for(int li = 0; li < grid::element::npoints; ++li) {
                        int const i = li + iorigin[0];
                        int const j = lj + iorigin[1];
                        int const k = lk + iorigin[2];
                        int const ijk = CCTK_GFINDEX3D(clawsolver::_M_cctkGH,
                                i, j, k);

                        int const idx = li*stride[0] + lj*stride[1] +
                            lk*stride[2];
                        int const xidx = li*xstride[0] + lj*xstride[1] +
                            lk*xstride[2];

                        for(int c = 0; c < claw::nequations; ++c) {
                            clawsolver::_M_RHS[c][ijk] -=
                                (scratch[3][c][xidx] + scratch[4][c][xidx] +
                                 scratch[5][c][xidx])*imass[idx];
                        }
                    }
                }
            }
        }
    protected:
        grid const mesh;
};

} // namespace

#endif
