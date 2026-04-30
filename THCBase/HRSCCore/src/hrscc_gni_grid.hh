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


#ifndef HRSCC_GNI_GRID_HH
#define HRSCC_GNI_GRID_HH

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>

#include <cctk.h>

#include <hrscc_gridinfo.hh>
#include <hrscc_typedefs.hh>

#include <utils.hh>

namespace hrscc {

//! A rather basic grid class for Galerkin-Numerical Integration schemes
template<typename element_t>
class GNIGrid {
    public:
        typedef element_t element;

        //! the constructor
        GNIGrid(cGH const * const cctkGH): _M_cctkGH(cctkGH) {
            _M_coordinates[0] = static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                        _M_cctkGH, 0, gridinfo::idx_coordinates[0]));
            _M_coordinates[1] = static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                        _M_cctkGH, 0, gridinfo::idx_coordinates[1]));
            _M_coordinates[2] = static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                        _M_cctkGH, 0, gridinfo::idx_coordinates[2]));
            _M_coordinates[3] = static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                        _M_cctkGH, 0, gridinfo::idx_coordinates[3]));

            // size of an element on this refinement level
            delta_space[0] =
              (_M_cctkGH->cctk_delta_space[0] * element::npoints /
               _M_cctkGH->cctk_levfac[0]);
            delta_space[1] =
              (_M_cctkGH->cctk_delta_space[1] * element::npoints /
               _M_cctkGH->cctk_levfac[1]);
            delta_space[2] =
              (_M_cctkGH->cctk_delta_space[2] * element::npoints /
               _M_cctkGH->cctk_levfac[2]);

            // ensure the grid has the correct size, and that Cactus
            // parallelised the grid at an element boundary
            assert ((_M_cctkGH->cctk_lsh[0] - 2) % element::npoints == 0);
            assert ((_M_cctkGH->cctk_lsh[1] - 2) % element::npoints == 0);
            assert ((_M_cctkGH->cctk_lsh[2] - 2) % element::npoints == 0);

            // number of elements in this component
            nelem[0] = (_M_cctkGH->cctk_lsh[0] - 2) / element::npoints;
            nelem[1] = (_M_cctkGH->cctk_lsh[1] - 2) / element::npoints;
            nelem[2] = (_M_cctkGH->cctk_lsh[2] - 2) / element::npoints;

            // assert that the Cactus index space is aligned at an
            // element boundary
            assert (_M_cctkGH->cctk_lbnd[0] % element::npoints == 0);
            assert (_M_cctkGH->cctk_lbnd[1] % element::npoints == 0);
            assert (_M_cctkGH->cctk_lbnd[2] % element::npoints == 0);

            // coordinate location of origin of this component
            for (int d=0; d<3; ++d) {
                CCTK_REAL const cactus_delta_space =
                    _M_cctkGH->cctk_delta_space[d] / _M_cctkGH->cctk_levfac[d];
                CCTK_REAL const origin_space =
                    _M_cctkGH->cctk_origin_space[d] + cactus_delta_space *
                    _M_cctkGH->cctk_levoff[d] / _M_cctkGH->cctk_levoffdenom[d];
                int const lbnd = _M_cctkGH->cctk_lbnd[d];
                // we add 1/2 to lbnd, because this component contains a
                // layer of ghost or boundary points, and (from Cactus's
                // point of view) the interior begins halfway between
                // this layer and the next.
                origin[d] = origin_space + (lbnd + 0.5) * cactus_delta_space;
            }

            // coordinate location of last point of this component
            end[0] = origin[0] + nelem[0]*delta_space[0];
            end[1] = origin[1] + nelem[1]*delta_space[1];
            end[2] = origin[2] + nelem[2]*delta_space[2];

            // index strides for array accesses
            stride[0] = CCTK_GFINDEX3D(_M_cctkGH, 1, 0, 0);
            stride[1] = CCTK_GFINDEX3D(_M_cctkGH, 0, 1, 0);
            stride[2] = CCTK_GFINDEX3D(_M_cctkGH, 0, 0, 1);

            // poison field
            _M_evaluator_point[0] = std::numeric_limits<CCTK_REAL>::quiet_NaN();
            _M_evaluator_point[1] = std::numeric_limits<CCTK_REAL>::quiet_NaN();
            _M_evaluator_point[2] = std::numeric_limits<CCTK_REAL>::quiet_NaN();
        }

        //! setups the coordinate vectors for Cactus
        void setup_coordinates() {
            int const gridsize = _M_cctkGH->cctk_lsh[0]*_M_cctkGH->cctk_lsh[1]*
                _M_cctkGH->cctk_lsh[2];
#pragma omp parallel
            {
                // Setup coordinates in the bulk
                UTILS_LOOP3(setup_coord,
                        ek, 0, nelem[2],
                        ej, 0, nelem[1],
                        ei, 0, nelem[0]) {
                    int const elem[3] = {ei, ej, ek};
                    CCTK_REAL * x = this->restriction(elem,
                            _M_coordinates[0]);
                    CCTK_REAL * y = this->restriction(elem,
                            _M_coordinates[1]);
                    CCTK_REAL * z = this->restriction(elem,
                            _M_coordinates[2]);

                    for(int lk = 0; lk < element::npoints; ++lk)
                    for(int lj = 0; lj < element::npoints; ++lj)
                    for(int li = 0; li < element::npoints; ++li) {
                        int const idx = li*stride[0]+lj*stride[1]+lk*stride[2];
                        assert(idx + this->element_offset(elem) < gridsize);

                        x[idx] = origin[0] + delta_space[0]*(elem[0] +
                                0.5*(1.0 + element::node[li]));
                        y[idx] = origin[1] + delta_space[1]*(elem[1] +
                                0.5*(1.0 + element::node[lj]));
                        z[idx] = origin[2] + delta_space[2]*(elem[2] +
                                0.5*(1.0 + element::node[lk]));
                    }
                } UTILS_ENDLOOP3(setup_coord);

                // Ghost-regions at x = 0, x = xmax
                for(int k = 0; k < _M_cctkGH->cctk_lsh[2]; ++k)
                for(int j = 0; j < _M_cctkGH->cctk_lsh[1]; ++j) {
                    int const idx0 = CCTK_GFINDEX3D(_M_cctkGH, 0, j, k);
                    int const idx1 = CCTK_GFINDEX3D(_M_cctkGH, 1, j, k);

                    _M_coordinates[0][idx0] = origin[0];
                    _M_coordinates[1][idx0] = _M_coordinates[1][idx1];
                    _M_coordinates[2][idx0] = _M_coordinates[2][idx1];

                    int const idxN1 = CCTK_GFINDEX3D(_M_cctkGH,
                            _M_cctkGH->cctk_lsh[0]-1, j, k);
                    int const idxN0 = CCTK_GFINDEX3D(_M_cctkGH,
                            _M_cctkGH->cctk_lsh[0]-2, j, k);

                    _M_coordinates[0][idxN1] = end[0];
                    _M_coordinates[1][idxN1] = _M_coordinates[1][idxN0];
                    _M_coordinates[2][idxN1] = _M_coordinates[2][idxN0];
                }

                // Ghost-regions at y = 0, y = ymax
                for(int k = 0; k < _M_cctkGH->cctk_lsh[2]; ++k)
                for(int i = 0; i < _M_cctkGH->cctk_lsh[0]; ++i) {
                    int const idx0 = CCTK_GFINDEX3D(_M_cctkGH, i, 0, k);
                    int const idx1 = CCTK_GFINDEX3D(_M_cctkGH, i, 1, k);

                    _M_coordinates[0][idx0] = _M_coordinates[0][idx1];
                    _M_coordinates[1][idx0] = origin[1];
                    _M_coordinates[2][idx0] = _M_coordinates[2][idx1];

                    int const idxN1 = CCTK_GFINDEX3D(_M_cctkGH,
                            i, _M_cctkGH->cctk_lsh[1]-1, k);
                    int const idxN0 = CCTK_GFINDEX3D(_M_cctkGH,
                            i, _M_cctkGH->cctk_lsh[1]-2, k);

                    _M_coordinates[0][idxN1] = _M_coordinates[0][idxN0];
                    _M_coordinates[1][idxN1] = end[1];
                    _M_coordinates[2][idxN1] = _M_coordinates[2][idxN0];
                }

                // Ghost-regions at z = 0, z = zmax
                for(int j = 0; j < _M_cctkGH->cctk_lsh[1]; ++j)
                for(int i = 0; i < _M_cctkGH->cctk_lsh[0]; ++i) {
                    int const idx0 = CCTK_GFINDEX3D(_M_cctkGH, i, j, 0);
                    int const idx1 = CCTK_GFINDEX3D(_M_cctkGH, i, j, 1);

                    _M_coordinates[0][idx0] = _M_coordinates[0][idx1];
                    _M_coordinates[1][idx0] = _M_coordinates[1][idx1];
                    _M_coordinates[2][idx0] = origin[2];

                    int const idxN1 = CCTK_GFINDEX3D(_M_cctkGH,
                            i, j, _M_cctkGH->cctk_lsh[2]-1);
                    int const idxN0 = CCTK_GFINDEX3D(_M_cctkGH,
                            i, j, _M_cctkGH->cctk_lsh[2]-2);

                    _M_coordinates[0][idxN1] = _M_coordinates[0][idxN0];
                    _M_coordinates[1][idxN1] = _M_coordinates[1][idxN0];
                    _M_coordinates[2][idxN1] = end[2];
                }

                // Setup the radial coordinate
#pragma omp for
                for(int idx = 0; idx < _M_cctkGH->cctk_lsh[0]*
                        _M_cctkGH->cctk_lsh[1]*_M_cctkGH->cctk_lsh[2]; ++idx) {
                    _M_coordinates[3][idx] = std::sqrt(
                            utils::pow<2>(_M_coordinates[0][idx]) +
                            utils::pow<2>(_M_coordinates[1][idx]) +
                            utils::pow<2>(_M_coordinates[2][idx]));
                }
            }
        }

        //! setups the weights for the integration
        void setup_weights() const {
            CCTK_REAL * weight = static_cast<CCTK_REAL *>(
                    CCTK_VarDataPtrI(_M_cctkGH, 0, gridinfo::idx_weight));
            int const gridsize = _M_cctkGH->cctk_lsh[0]*_M_cctkGH->cctk_lsh[1]*
                _M_cctkGH->cctk_lsh[2];
#pragma omp parallel
            {
                // Setup weights in the bulk
                UTILS_LOOP3(setup_weights,
                        ek, 0, nelem[2],
                        ej, 0, nelem[1],
                        ei, 0, nelem[0]) {
                    int const elem[3] = {ei, ej, ek};
                    CCTK_REAL * w = this->restriction(elem, weight);
                    CCTK_REAL const coarse_space[3] = {
                        delta_space[0]*_M_cctkGH->cctk_levfac[0],
                        delta_space[1]*_M_cctkGH->cctk_levfac[1],
                        delta_space[2]*_M_cctkGH->cctk_levfac[2]};
                    for(int lk = 0; lk < element::npoints; ++lk)
                    for(int lj = 0; lj < element::npoints; ++lj)
                    for(int li = 0; li < element::npoints; ++li) {
                        int const idx = li*stride[0]+lj*stride[1]+lk*stride[2];
                        assert(idx + this->element_offset(elem) < gridsize);

                        w[idx] *= 0.5*coarse_space[0]*element::weight[li] *
                                  0.5*coarse_space[1]*element::weight[lj] *
                                  0.5*coarse_space[2]*element::weight[lk];
                    }
                } UTILS_ENDLOOP3(setup_weights);

                // Ghost-regions at x = 0, x = xmax
                for(int k = 0; k < _M_cctkGH->cctk_lsh[2]; ++k)
                for(int j = 0; j < _M_cctkGH->cctk_lsh[1]; ++j) {
                    int const idx = CCTK_GFINDEX3D(_M_cctkGH, 0, j, k);
                    weight[idx] = 0.0;

                    int const idxN = CCTK_GFINDEX3D(_M_cctkGH,
                            _M_cctkGH->cctk_lsh[0]-1, j, k);
                    weight[idxN] = 0.0;
                }

                // Ghost-regions at y = 0, y = ymax
                for(int k = 0; k < _M_cctkGH->cctk_lsh[2]; ++k)
                for(int i = 0; i < _M_cctkGH->cctk_lsh[0]; ++i) {
                    int const idx0 = CCTK_GFINDEX3D(_M_cctkGH, i, 0, k);
                    weight[idx0] = 0.0;

                    int const idxN = CCTK_GFINDEX3D(_M_cctkGH,
                            i, _M_cctkGH->cctk_lsh[1]-1, k);
                    weight[idxN] = 0.0;
                }

                // Ghost-regions at z = 0, z = zmax
                for(int j = 0; j < _M_cctkGH->cctk_lsh[1]; ++j)
                for(int i = 0; i < _M_cctkGH->cctk_lsh[0]; ++i) {
                    int const idx = CCTK_GFINDEX3D(_M_cctkGH, i, j, 0);
                    weight[idx] = 0.0;

                    int const idxN = CCTK_GFINDEX3D(_M_cctkGH,
                            i, j, _M_cctkGH->cctk_lsh[2]-1);
                    weight[idxN] = 0.0;
                }
            }
        }

        //! \brief finds the absolute index of the given point, i.e. in an
        //! indexing where duplicated dof are not taken into account
        /*!
         *  NOTE: this does not work for the P0 case!
         */
        inline int absolute_index(
                //! [in] 1D index
                int idx
                ) const {
            assert(element::npoints > 0);
            int li = idx % element::npoints;
            int el = idx / element::npoints;
            return el*(element::npoints - 1) + std::max(li - 1, 0);
        }

        //! \brief finds the index of the first point (in the lexicographical
        //! order) of the given element
        template<policy::direction_t dir>
        inline int element_iorigin(
                //! [in] index of the element
                int const elem[3]
                ) const {
            return elem[dir]*element::npoints + 1;
        }

        //! \brief returns the coordinate position of the first point (in the
        //! lexicographical order) of the given element
        template<policy::direction_t dir>
        inline CCTK_REAL element_origin(
                //! [in] index of the element
                int const elem[3]
                ) const {
            return origin[dir] + delta_space[dir]*elem[dir];
        }

        //! \brief finds the 3d index of the first point (in the lexicographical
        //! order) of a given element
        inline int element_offset(
                //! [in] index of the element
                int const elem[3]
                ) const {
            int const i = this->element_iorigin<policy::x>(elem);
            int const j = this->element_iorigin<policy::y>(elem);
            int const k = this->element_iorigin<policy::z>(elem);
            return CCTK_GFINDEX3D(_M_cctkGH, i, j, k);
        }

        //! mapping from the reference element into the physical space
        template<policy::direction_t dir>
        inline CCTK_REAL element_map(
                //! [in] index of the element
                int const elem[3],
                //! [in] point in the reference element
                CCTK_REAL const refpos[3]
                ) const {
            return this->element_origin<dir>(elem) + 0.5*delta_space[dir] *
                (1.0 + refpos[dir]);
        }

        //! mapping from the physical space into the reference element
        template<policy::direction_t dir>
        inline CCTK_REAL element_imap(
                //! [in] index of the element
                int const elem[3],
                //! [in] point in the physical space
                CCTK_REAL const physpos[3]
                ) const {
            return 2.0*(physpos[dir] - this->element_origin<dir>(elem)) /
                delta_space[dir];
        }

        //! restrict a given grid-function to the selected element
        inline CCTK_REAL * restriction(
                //! [in] index of the element
                int const elem[3],
                //! [in] target function
                CCTK_REAL * grid_function
                ) const {
            return &grid_function[this->element_offset(elem)];
        }
        //! restrict a given grid-function to the selected element
        inline CCTK_REAL const * restriction(
                //! [in] index of the element
                int const elem[3],
                //! [in] target function
                CCTK_REAL const * grid_function
                ) const {
            return const_cast<CCTK_REAL const *>(this->restriction(elem,
                        const_cast<CCTK_REAL *>(grid_function)));
        }

        //! computes the pseudo-spectral derivative of a grid function
        /*!
         *  \tparam dir direction of the derivative
         *
         *  In the reference element we have
         *  \f[
                [\partial_x u]_{ijk} = \sum_{q} u_{qjk}\ l'_q(x_i), \quad
                [\partial_y u]_{ijk} = \sum_{q} u_{iqk}\ l'_q(y_j), \quad
                [\partial_z u]_{ijk} = \sum_{q} u_{ijq}\ l'_q(z_k).
            \f]
         *
         *  This method is supposed to be used to compute derivative of
         *  quantities appearing in the source terms or for data-analysis.
         *  This is *NOT* the proper Galerkin derivative and should not
         *  be used for the flux terms.
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically and divides the result by the mass matrix.
         */
        template<policy::direction_t dir>
        CCTK_REAL diff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int const i, int const j, int const k
                ) const {
            int elem[3];
            elem[0] = (i - 1) / element::npoints;
            elem[1] = (j - 1) / element::npoints;
            elem[2] = (k - 1) / element::npoints;

            int const lindex[3] = {
                i - this->element_iorigin<policy::x>(elem),
                j - this->element_iorigin<policy::y>(elem),
                k - this->element_iorigin<policy::z>(elem)};

            return element::template diff<dir>(delta_space,
                    this->restriction(elem, grid_function),
                    stride, lindex);
        }

        //! computes the weak derivative of a grid function in a point
        /*!
         *  \tparam dir direction of the derivative
         *
         *  This is a weak derivative computed using a centered flux
         *
         *  In the reference element the pseudo-spectral derivative is
         *  \f[
                [\partial_x u]_{ijk} = \sum_{q} u_{qjk}\ l'_q(x_i), \quad
                [\partial_y u]_{ijk} = \sum_{q} u_{iqk}\ l'_q(y_j), \quad
                [\partial_z u]_{ijk} = \sum_{q} u_{ijq}\ l'_q(z_k).
            \f]
         *
         *  The weak derivative is computed as
         *  \f[
                D u = \partial u + \llbracket u \rrbracket \cdot \vec{d}
                \mathbin{\vrule height 1.6ex depth 0pt width
                 0.13ex\vrule height 0.13ex depth 0pt width 1.3ex} j_u,
            \f]
         *  where \f$j_u\f$ is the jump set of \f$u\f$, i.e. the inter-element
         *  boundary and \f$\llbracket u \rrbracket\f$ is the jump-operator
         *  \f[
         *      \llbracket u \rrbracket = u_1 \vec{n}_1 + u_2 \vec{n_2},
         *  \f]
         *  such that the jump in a given direction, \f$\vec{d}\f$, is defined
         *  as
         *  \f[
                \llbracket u \rrbracket \cdot \vec{d} =
                (u_1 \vec{n}_1 + u_2 \vec{n}_2) \cdot \vec{d} = u_L - u_R,
            \f]
         *  where \f$u_L\f$ and \f$u_R\f$ are the left and right states
         *
         *  This method is supposed to be used to compute derivative of
         *  quantities appearing in the source terms or for data-analysis.
         *  This is *NOT* the proper Galerkin derivative and should not
         *  be used for the flux terms.
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically and divides the result by the mass matrix.
         */
        template<policy::direction_t dir>
        CCTK_REAL wdiff(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point index
                int const i, int const j, int const k
                ) const {
            int elem[3];
            elem[0] = (i - 1) / element::npoints;
            elem[1] = (j - 1) / element::npoints;
            elem[2] = (k - 1) / element::npoints;

            int const lindex[3] = {
                i - this->element_iorigin<policy::x>(elem),
                j - this->element_iorigin<policy::y>(elem),
                k - this->element_iorigin<policy::z>(elem)};

            return element::template wdiff<dir>(delta_space,
                    this->restriction(elem, grid_function),
                    stride, lindex);
        }
        //! \brief computes the weighted value of the jump of a given
        //! grid function in a point
        /*!
         *  The jump in a given direction, \f$\vec{d}\f$, is defined as
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
            int elem[3];
            elem[0] = (i - 1) / element::npoints;
            elem[1] = (j - 1) / element::npoints;
            elem[2] = (k - 1) / element::npoints;

            int const lindex[3] = {
                i - this->element_iorigin<policy::x>(elem),
                j - this->element_iorigin<policy::y>(elem),
                k - this->element_iorigin<policy::z>(elem)};

            return element::template jump<dir>(delta_space,
                    this->restriction(elem, grid_function),
                    stride, lindex);
        }

        //! integrates a grid function using a Gaussian quadrature formula
        CCTK_REAL integrate(
                //! [in] target grid function
                CCTK_REAL const * grid_function
                ) const {
            CCTK_REAL s = 0;
#pragma omp parallel for collapse(3) reduction(+:s)
            for(int k = 0; k < nelem[2]; ++k)
            for(int j = 0; j < nelem[1]; ++j)
            for(int i = 0; i < nelem[0]; ++i) {
                int const elem[3] = {i, j, k};
                s += element::integrate(delta_space,
                        this->restriction(elem, grid_function),
                        stride);
            }
            return s;
        }

        //! performs a P0 projection of a given grid function
        void flatten(
                //! [in,out] target grid function
                CCTK_REAL * grid_function
                ) const {
            CCTK_REAL const ones[3] = {1.0, 1.0, 1.0};
#pragma omp parallel
            {
                UTILS_LOOP3(flatten,
                        k, 0, nelem[2],
                        j, 0, nelem[1],
                        i, 0, nelem[0]) {
                    int const elem[3] = {i, j, k};
                    CCTK_REAL * gf = this->restriction(elem, grid_function);
                    CCTK_REAL const avg = element::integrate(ones, gf, stride);

                    for(int lk = 0; lk < element::npoints; ++lk)
                    for(int lj = 0; lj < element::npoints; ++lj)
                    for(int li = 0; li < element::npoints; ++li) {
                        gf[li*stride[0] + lj*stride[1] + lk*stride[2]] = avg;
                    }
                } UTILS_ENDLOOP3(flatten);
            }
        }

        //! locate the element containing the given point
        inline void locate(
                //! [in] target point
                CCTK_REAL const point[3],
                //! [out] element
                int elem[3]
                ) const {
            elem[0] = static_cast<int>((point[0] - origin[0])/delta_space[0]);
            elem[1] = static_cast<int>((point[1] - origin[1])/delta_space[1]);
            elem[2] = static_cast<int>((point[2] - origin[2])/delta_space[2]);
        }

        //! finds the closest *interior* point to the given point
        void project_point(
                CCTK_REAL const point[3],
                CCTK_REAL ppoint[3]
                ) const {
            CCTK_REAL const delta[3] = {
                std::numeric_limits<CCTK_REAL>::epsilon()*(end[0] - origin[0]),
                std::numeric_limits<CCTK_REAL>::epsilon()*(end[1] - origin[1]),
                std::numeric_limits<CCTK_REAL>::epsilon()*(end[2] - origin[2])
            };

            for(int d = 0; d < 3; ++d) {
                if(point[d] > origin[d] && point[d] < end[d]) {
                    ppoint[d] = point[d];
                }
                else if(point[d] <= origin[d]) {
                    ppoint[d] = origin[d] + delta[d];
                }
                else {
                    ppoint[d] = end[d] - delta[d];
                }
            }
        }
        //! evaluates a grid function in a point
        /*!
         *  This method evaluates the grid function using its expansion on the
         *  collocation basis:
         *  \f[
                u(x) = \sum_{ijk} u_{ijk} l_i(x) l_j(y) l_k(z).
            \f]
         *  The values of the basis elements in the previous expression,
         *  \f$l_i(x)\f$, are computed using the definition:
         *  \f[
                l_i(x) = \prod_{k\neq i} \frac{x-x_k}{x_i-x_k} =
                    \frac{1}{q_i} \prod_{k\neq i} (x - x_k),
            \f]
         *  where the coefficients \f$1/q_i\f$ are prestored.
         *
         *  This method takes care of handling the mappings from/to the
         *  reference element.
         *
         *  The values of the basis functions at the interpolation point are
         *  cached and re-used for successive evaluations of different functions
         *  at the same grid location.
         *
         *  If the given point is not part of the grid this returns a NaN.
         *
         *  This method is *NOT* thread safe. In order to use in a multi-
         *  threaded environment the user should take care of creating a
         *  separate instance of this class for each thread planning to
         *  call this method.
         *
         *  \see GLLElement::icoeff
         */
        CCTK_REAL eval(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point
                CCTK_REAL const point[3]
                ) const {
            if(    (_M_evaluator_point[0] != point[0])
                || (_M_evaluator_point[1] != point[1])
                || (_M_evaluator_point[2] != point[2])) {
                this->locate(point, _M_evaluator_elem);
                for(int d = 0; d < 3; ++d) {
                    if(_M_evaluator_elem[d] < 0 ||
                            _M_evaluator_elem[d] >= nelem[d]) {
                        return std::numeric_limits<CCTK_REAL>::signaling_NaN();
                    }
                }

                _M_evaluator_point[0] = point[0];
                _M_evaluator_point[1] = point[1];
                _M_evaluator_point[2] = point[2];

                // Point in the reference element
                CCTK_REAL const x[3] = {
                    this->element_imap<policy::x>(_M_evaluator_elem,
                            _M_evaluator_point),
                    this->element_imap<policy::y>(_M_evaluator_elem,
                            _M_evaluator_point),
                    this->element_imap<policy::z>(_M_evaluator_elem,
                            _M_evaluator_point)
                };

                for(int d = 0; d < 3; ++d) {
                    for(int p = 0; p < element::npoints; ++p) {
                        _M_evaluator_operator[d][p] = element::icoeff[p];
                        for(int q = 0; q < p; ++q) {
                            _M_evaluator_operator[d][p] =
                                _M_evaluator_operator[d][p] *
                                    (x[d] - element::x[q]);
                        }
                        for(int q = p + 1; q < element::npoints; ++q) {
                            _M_evaluator_operator[d][p] =
                                _M_evaluator_operator[d][p] *
                                    (x[d] - element::x[q]);
                        }
                    }
                }
            }
            CCTK_REAL const * func = this->restriction(_M_evaluator_elem,
                    grid_function);
            CCTK_REAL sum = 0;
            for(int k = 0; k < element::npoints; ++k)
            for(int j = 0; j < element::npoints; ++j)
            for(int i = 0; i < element::npoints; ++i) {
                sum += func[i*stride[0] + j*stride[1] + k*stride[2]] *
                    _M_evaluator_operator[0][i] * _M_evaluator_operator[1][j] *
                    _M_evaluator_operator[2][k];
            }
            return sum;
        }
        //! 0th order interpolation of a grid function in a point
        /*
         *  If the given point is not part of the grid this returns a NaN.
         *
         *  This method is *NOT* thread safe. In order to use in a multi-
         *  threaded environment the user should take care of creating a
         *  separate instance of this class for each thread planning to
         *  call this method.
         */
        CCTK_REAL eval0(
                //! [in] target grid function
                CCTK_REAL const * grid_function,
                //! [in] grid point
                CCTK_REAL const point[3]
                ) const {
            if(    (_M_evaluator_point[0] != point[0])
                || (_M_evaluator_point[1] != point[1])
                || (_M_evaluator_point[2] != point[2])) {
                this->locate(point, _M_evaluator_elem);
                for(int d = 0; d < 3; ++d) {
                    if(_M_evaluator_elem[d] < 0 ||
                            _M_evaluator_elem[d] >= nelem[d]) {
                        return std::numeric_limits<CCTK_REAL>::signaling_NaN();
                    }
                }

                _M_evaluator_point[0] = point[0];
                _M_evaluator_point[1] = point[1];
                _M_evaluator_point[2] = point[2];

                int const iorigin[3] = {
                    element_iorigin<0>(_M_evaluator_elem),
                    element_iorigin<1>(_M_evaluator_elem),
                    element_iorigin<2>(_M_evaluator_elem)
                };
                int iindex[3];
                for(int dir = 0; dir < 3; ++dir) {
                    iindex[0] = iorigin[0];
                    iindex[1] = iorigin[1];
                    iindex[2] = iorigin[2];
                    int ijk = CCTK_GFINDEX3D(_M_cctkGH, iindex[0], iindex[1],
                            iindex[2]);

                    int & i = iindex[dir];
                    while(_M_coordinates[dir][ijk] < point[dir]) {
                        ++i;
                        ijk = CCTK_GFINDEX3D(_M_cctkGH, iindex[0], iindex[1],
                            iindex[2]);
                    }
                    if(std::abs(_M_coordinates[dir][ijk] - point[dir]) <
                        std::abs(_M_coordinates[dir][ijk - stride[dir]] -
                            point[dir])) {
                        _M_evaluator_index[dir] = i;
                    }
                    else {
                        _M_evaluator_index[dir] = i - 1;
                    }
                }
            }
            int idx = CCTK_GFINDEX3D(_M_cctkGH, _M_evaluator_index[0],
                    _M_evaluator_index[1], _M_evaluator_index[2]);

            return grid_function[idx];
        }

        //! re-projects a given grid function to make it continuous
        void smooth(
                //! [in,out] target grid function
                CCTK_REAL * grid_function
                ) const {
            int const siz = _M_cctkGH->cctk_lsh[0]*_M_cctkGH->cctk_lsh[1]*
            _M_cctkGH->cctk_lsh[2];
            for(int k = 1; k < _M_cctkGH->cctk_lsh[2]-1; ++k)
            for(int j = 1; j < _M_cctkGH->cctk_lsh[1]-1; ++j)
            for(int i = 1; i < _M_cctkGH->cctk_lsh[0]-1; ++i) {
                int const ijk = CCTK_GFINDEX3D(_M_cctkGH, i, j, k);

                bool stencil[3][3][3];

                CCTK_REAL avg = 0;
                CCTK_REAL avg_norm = 0;
                for(int lk = -1; lk <= 1; ++lk)
                for(int lj = -1; lj <= 1; ++lj)
                for(int li = -1; li <= 1; ++li) {
                    int lijk = li*stride[0] + lj*stride[1] + lk*stride[2];
                    assert(lijk + ijk < siz && lijk + ijk >= 0);
                    if(
                           absolute_index(i) == absolute_index(i + li)
                        && absolute_index(j) == absolute_index(j + lj)
                        && absolute_index(k) == absolute_index(k + lk)
                        && !(i + li == 0 && _M_cctkGH->cctk_bbox[0])
                        && !(i + li == _M_cctkGH->cctk_lsh[0]-1 &&
                            _M_cctkGH->cctk_bbox[1])
                        && !(j + lj == 0 && _M_cctkGH->cctk_bbox[2])
                        && !(j + lj == _M_cctkGH->cctk_lsh[1]-1 &&
                            _M_cctkGH->cctk_bbox[3])
                        && !(k + lk == 0 && _M_cctkGH->cctk_bbox[4])
                        && !(k + lk == _M_cctkGH->cctk_lsh[2]-1 &&
                            _M_cctkGH->cctk_bbox[5])
                      ) {
                        avg += grid_function[ijk + lijk];
                        avg_norm += 1;
                        stencil[li+1][lj+1][lk+1] = true;
                    }
                    else {
                        stencil[li+1][lj+1][lk+1] = false;
                    }
                }
                assert(avg == avg);
                avg = avg / avg_norm;

                for(int lk = -1; lk <= 1; ++lk)
                for(int lj = -1; lj <= 1; ++lj)
                for(int li = -1; li <= 1; ++li) {
                    if(stencil[li+1][lj+1][lk+1]) {
                        int lijk = li*stride[0] + lj*stride[1] +
                            lk*stride[2];
                        grid_function[ijk + lijk] = avg;
                    }
                }
            }
        }

        //! filter a given grid function
        void filter(
                //! [in] filter values on the whole grid
                CCTK_REAL const filter[element::npoints],
                //! [in] filter strength (should be constant over each element)
                CCTK_REAL const * strength,
                //! [in, out] grid function to filter
                CCTK_REAL * grid_function
                ) const {
#pragma omp parallel
            {
                CCTK_REAL phi[element::npoints];
                UTILS_LOOP3(filter,
                        k, 0, nelem[2],
                        j, 0, nelem[1],
                        i, 0, nelem[0]) {
                    int const elem[3] = {i, j, k};

                    CCTK_REAL const w = *this->restriction(elem, strength);
                    for(int a = 0; a < element::npoints; ++a) {
                        phi[a] = std::pow(filter[a], w);
                    }

                    CCTK_REAL * gf = this->restriction(elem, grid_function);
                    element::filter(phi, gf, stride);
                } UTILS_ENDLOOP3(filter);
            }
        }

        //! position of the first grid point (in the lexicographical order)
        CCTK_REAL origin[3];
        //! position of the last grid point (in the lexicographical order)
        CCTK_REAL end[3];
        //! grid spacing in the three directions
        CCTK_REAL delta_space[3];
        //! number of elements in the three directions
        int nelem[3];
        //! stride used to access the data in the single elements
        int stride[3];
    private:
        //! Cactus grid hierarchy
        cGH const * const _M_cctkGH;
        //! Cactus coordinates grid functions: {x, y, z, r}
        CCTK_REAL * _M_coordinates[4];
        //! current evaluator point
        mutable CCTK_REAL _M_evaluator_point[3];
        //! index of the current evaluator element
        mutable int _M_evaluator_elem[3];
        //! index of the grid point used by the 0th order evaluator
        mutable int _M_evaluator_index[3];
        //! cached values of the collocation basis at the evaluator point
        mutable CCTK_REAL _M_evaluator_operator[3][element::npoints];
};

} // namespace

#endif
