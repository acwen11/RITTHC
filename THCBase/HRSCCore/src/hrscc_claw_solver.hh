//  HRSCCore: HRSC methods for Cactus
//  Copyright (C) 2012, David Radice <david.radice@aei.mpg.de>
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


#ifndef HRSCC_CLAW_SOLVER_HH
#define HRSCC_CLAW_SOLVER_HH

#include <cctk.h>

#include <hrscc_gridinfo.hh>
#include <hrscc_metric_info.hh>
#include <hrscc_macro.hh>
#include <hrscc_observer.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! Conservation Laws solver
/*!
 *  This is an "abstract" class representing an abstract HRSC numerical scheme
 *
 *  \tparam law_t a class describing a conservation law, it should inherit from
 *          CLaw
 *  \tparam method_t should be a class inheriting from CLawSolver defining
 *          a numerical method
 */
template<typename law_t, typename method_t>
class CLawSolver {
    public:
        typedef law_t law;
        typedef CLaw<law> claw;
        typedef method_t method;

        //! number of ghost zones to use
        enum {nghostzones = traits<method>::nghostzones};
        //! is this method positivity preserving
        static bool const pospres = traits<method>::pospres;
        //! does this method support refluxing
        static bool const refluxing = traits<method>::refluxing;

        //! constructor
        CLawSolver(
                //! [in] Cactus grid hierarchy
                cGH const * const cctkGH
                ): _M_cctkGH(cctkGH) {
            _M_claw = new law;

            _M_coordinates[0] = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0,
                        gridinfo::idx_coordinates[0]));
            _M_coordinates[1] = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0,
                        gridinfo::idx_coordinates[1]));
            _M_coordinates[2] = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0,
                        gridinfo::idx_coordinates[2]));

            _M_lapse = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, metric_info::idx_lapse));
            _M_metric[0] = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0,
                        metric_info::idx_metric[0]));
            _M_metric[1] = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0,
                        metric_info::idx_metric[1]));
            _M_metric[2] = static_cast<CCTK_REAL const *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0,
                        metric_info::idx_metric[2]));

            for(int v = 0; v < claw::nequations; ++v) {
                _M_conserved[v] = static_cast<CCTK_REAL *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, claw::conserved_idx[v]));
                _M_primitive[v] = static_cast<CCTK_REAL *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, claw::primitive_idx[v]));
                _M_RHS[v] = static_cast<CCTK_REAL *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, claw::rhs_idx[v]));
            }
            for(int v = 0; v < claw::nexternal; ++v) {
                _M_field[v] = static_cast<CCTK_REAL *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, claw::field_idx[v]));
            }

            CCTK_REAL ** p = &_M_grid_variable[0];
            for(int v = 0; v < claw::nequations; ++v) {
                *(p++) = _M_conserved[v];
            }
            for(int v = 0; v < claw::nequations; ++v) {
                *(p++) = _M_primitive[v];
            }
            for(int v = 0; v < claw::nequations; ++v) {
                *(p++) = _M_RHS[v];
            }
            for(int v = 0; v < claw::nexternal; ++v) {
                *(p++) = _M_field[v];
            }

            for(int b = 0; b < claw::nbitmasks; ++b) {
                _M_bitmask[b] = static_cast<CCTK_INT *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, claw::bitmask_idx[b]));
            }

            if(refluxing) {
                for(int d = 0; d < 3; ++d)
                for(int v = 0; v < claw::nequations; ++v) {
                    _M_num_flux[d][v] = static_cast<CCTK_REAL *>(
                            CCTK_VarDataPtrI(_M_cctkGH, 0,
                                claw::num_flux_idx[3*v + d]));
                }
            }
        }
        //! destructor
        ~CLawSolver() {
            delete _M_claw;
        }

        //! get a read-only reference to the conservation law
        claw const * get_claw() const {
            return const_cast<claw const *>(_M_claw);
        }

        //! make an observer
        Observer<claw> * observer_alloc() const {
            return new Observer<claw>(_M_cctkGH, _M_coordinates,
                    _M_grid_variable, _M_bitmask);
        }
        //! de-allocate an observer
        void observer_free(Observer<claw> * observer) const {
            delete observer;
        }

        //! compute all the variables from the primitives on the grid
        void prim_to_all() {
#pragma omp parallel
            {
                Observer<claw> observer(_M_cctkGH, _M_coordinates,
                        _M_grid_variable, _M_bitmask);

                UTILS_LOOP3(prim_to_all,
                        k, 0, _M_cctkGH->cctk_lsh[2],
                        j, 0, _M_cctkGH->cctk_lsh[1],
                        i, 0, _M_cctkGH->cctk_lsh[0]) {
                    observer.jump_to_location(i, j, k);
                    _M_claw->prim_to_all(observer);
                    observer.record();
                } UTILS_ENDLOOP3(prim_to_all);
            }
        }

        //! compute all the variables from the conservatives on the grid
        void cons_to_all() {
#pragma omp parallel
            {
                Observer<claw> observer(_M_cctkGH, _M_coordinates,
                        _M_grid_variable, _M_bitmask);

                UTILS_LOOP3(cons_to_all,
                        k, 0, _M_cctkGH->cctk_lsh[2],
                        j, 0, _M_cctkGH->cctk_lsh[1],
                        i, 0, _M_cctkGH->cctk_lsh[0]) {
                    observer.jump_to_location(i, j, k);
                    _M_claw->cons_to_all(observer);
                    observer.record();
                } UTILS_ENDLOOP3(cons_to_all);
            }
        }

        //! \brief "virtual" method to compute the derivative of a grid
        //! function in a point
        /*!
         *  \tparam dir direction of the derivative
         */
        template<policy::direction_t dir>
        CCTK_REAL diff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int i, int j, int k
                ) const {
            return static_cast<method const *>(this)->method::template
                diff<dir>(grid_function, i, j, k);
        }

        //! \brief "virtual" method to compute the weak derivative of a grid
        //! function in a point
        /*!
         *  \tparam dir direction of the derivative
         */
        template<policy::direction_t dir>
        CCTK_REAL wdiff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int i, int j, int k
                ) const {
            return static_cast<method const *>(this)->method::template
                wdiff<dir>(grid_function, i, j, k);
        }
        //! \brief "virtual" method to compute the weighted value of the
        //! jump of a given grid function in a point
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
                int i, int j, int k
                ) const {
            return static_cast<method const *>(this)->method::template
                jump<dir>(grid_function, i, j, k);
        }

        //! set the RHS to zero
        /*!
         *  To be called before compute_rhs if the conservation law has no
         *  souce terms
         */
        void reset_rhs() {
#pragma omp parallel
            {
                int eq;
                int ijk;
                UTILS_LOOP3(reset_rhs,
                        k, 0, _M_cctkGH->cctk_lsh[2],
                        j, 0, _M_cctkGH->cctk_lsh[1],
                        i, 0, _M_cctkGH->cctk_lsh[0]) {
                    ijk = CCTK_GFINDEX3D(_M_cctkGH, i, j, k);
                    for(eq = 0; eq < claw::nequations; ++eq) {
                        _M_RHS[eq][ijk] = 0;
                    }
                } UTILS_ENDLOOP3(reset_rhs);
            }
        }
        //! "virtual" method to compute the RHS of the equations
        /*!
         *  This method should *ADD* to the RHS variables the values determined
         *  from the spatial discretization. The user should take care to either
         *  initialize those grid variables to zero before this is scheduled or
         *  to initialize them with the values of the sources of the balance
         *  law.
         *
         *  If CLaw::pureclaw is set to true then the initialization of the RHS
         *  is not necessary as CLawSolver will take care of setting it to zero
         */
        inline void compute_rhs() {
            if(claw::pure) {
                this->reset_rhs();
            }
            static_cast<method *>(this)->compute_rhs();
        }
    protected:
        //! cactus grid hierarchy
        cGH const * const _M_cctkGH;
        //! physical driver class
        claw * _M_claw;

        //! grid coordinates
        CCTK_REAL const * _M_coordinates[3];

        //! Lapse function
        CCTK_REAL const * _M_lapse;
        //! Diagonal components of the metric
        CCTK_REAL const * _M_metric[3];

        //! global index for all the grid variables registered by claw
        /*!
         *  this index does not include the numerical fluxes
         */
        CCTK_REAL * _M_grid_variable[claw::nvariables];
        //! alias for the conserved variables
        CCTK_REAL * _M_conserved[claw::nequations];
        //! alias for the primitive variables
        CCTK_REAL * _M_primitive[claw::nequations];
        //! alias for the grid variables containing the RHS
        CCTK_REAL * _M_RHS[claw::nequations];
        //! alias for the external fields registered by claw
        CCTK_REAL * _M_field[claw::nexternal];

        //! bitmasks requested by claw
        CCTK_INT  * _M_bitmask[claw::nbitmasks];

        //! numerical fluxes
        CCTK_REAL * _M_num_flux[3][claw::nequations];
};

} // namespace

#endif
