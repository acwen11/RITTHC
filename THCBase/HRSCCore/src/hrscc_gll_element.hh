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


#ifndef HRSCC_GLL_ELEMENT_HH
#define HRSCC_GLL_ELEMENT_HH

#include <algorithm>
#include <cmath>

#include <cctk.h>

#include <hrscc_config_par.hh>
#include <hrscc_typedefs.hh>

#include <utils_pow.hh>
#include <utils_sequence.hh>

//! Maximum order of the implemented elements
#define HRSCC_GLL_ELEMENT_MAX_ORDER 25

namespace hrscc {

//! Gauss-Legendre-Lobatto finite-element
/*!
 *  This class represents a one-dimensional reference element in \f$[-1,1]\f$.
 *
 *  \tparam element_order order of the polynomials used to represent the data
 */
template<int element_order>
class GLLElement {
    public:
        //! order of the element
        enum {order = element_order};
        //! number of collocation points in each direction
        enum {npoints = order + 1};
        //! total number of degrees of freedom within the element
        enum {ndof = npoints*npoints*npoints };

        //! evaluates the co-differential of a grid function
        /*!
         *  \tparam dir direction of the derivative
         *
         *  This method evaluates the co-differential of a function over
         *  a single element, i.e. it computes
         *  \f[
                [\delta u]_{ijk} = - \int u \partial l_{ijk} dx,
            \f]
         *  where \f$l_{ijk}(x,y,z) = l_i(x) l_j(y) l_k(z)\f$ is the collocation
         *  basis and \f$l_i\f$ is the one-dimensional collocation basis.
         *
         *  In particular, in the reference element, we have
            \f[
            [\delta_x u]_{ijk} = - \sum_q u_{qjk} w_q w_j w_k l'_i(x_q),\qquad
            [\delta_y u]_{ijk} = - \sum_q u_{iqk} w_i w_q w_k l'_j(y_q),\qquad
            [\delta_z u]_{ijk} = - \sum_q u_{ijq} w_i w_j w_q l'_k(z_q).
            \f]
         *
         *  Input and outputs grid-arrays should be offsetted so that the 0th
         *  element represents the value of the grid function in the first
         *  point of the element (in the lexicographical ordering).
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically.
         *
         *  \see GNIGrid::restriction
         */
        template<policy::direction_t dir>
        static void codiff(
                //! [in] size of the element in each direction
                CCTK_REAL const shape[3],
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access grid_function array
                int const stride_grid_function[3],
                //! [out] *restricted* codifferential of the grid function
                CCTK_REAL * codiff_grid_function,
                //! [in] stride used to access the codiff_grid_function array
                int const stride_codiff_grid_function[3]) {
            int index[3];
            for(index[2] = 0; index[2] < npoints; ++index[2])
            for(index[1] = 0; index[1] < npoints; ++index[1])
            for(index[0] = 0; index[0] < npoints; ++index[0]) {
                CCTK_REAL & du = codiff_grid_function[
                    index[0]*stride_codiff_grid_function[0] +
                    index[1]*stride_codiff_grid_function[1] +
                    index[2]*stride_codiff_grid_function[2]];
                du = 0;

                // moving index
                int mindex[3] = {index[0], index[1], index[2]};
                int & i = mindex[dir];
                for(i = 0; i < npoints; ++i) {
                    du -= grid_function[
                            mindex[0]*stride_grid_function[0] +
                            mindex[1]*stride_grid_function[1] +
                            mindex[2]*stride_grid_function[2]] *
                        diffop[i][index[dir]] *
                        weight[mindex[0]]*weight[mindex[1]]*weight[mindex[2]];
                }
                for(int d = 0; d < 3; ++d) {
                    du = du * (d == dir ? 1 : 0.5*shape[d]);
                }
            }
        }

        //! evaluates the pseudo-spectral derivative of a grid function
        /*!
         *  \tparam dir direction of the derivative
         *
         *  This method evaluates the differential of a function over
         *  a single element.
         *
         * In the reference element we have
         *  \f[
                [\partial_x u]_{ijk} = \sum_{q} u_{qjk}\ l'_q(x_i), \qquad
                [\partial_y u]_{ijk} = \sum_{q} u_{iqk}\ l'_q(y_j), \qquad
                [\partial_z u]_{ijk} = \sum_{q} u_{ijq}\ l'_q(z_k);
            \f]
         *  where \f$l_i\f$ is the one-dimensional collocation basis.
         *
         *  Input and outputs grid-arrays should be offsetted so that the 0th
         *  element represents the value of the grid function in the first
         *  point of the element (in the lexicographical ordering).
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically.
         *
         *  \see GNIGrid::restriction
         */
        template<policy::direction_t dir>
        static void diff(
                //! [in] size of the element in each direction
                CCTK_REAL const shape[3],
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access grid_function array
                int const stride_grid_function[3],
                //! [out] *restricted* differential of the grid function
                CCTK_REAL * diff_grid_function,
                //! [in] stride used to access the diff_grid_function array
                int const stride_diff_grid_function[3]) {
            int index[3];
            for(index[2] = 0; index[2] < npoints; ++index[2])
            for(index[1] = 0; index[1] < npoints; ++index[1])
            for(index[0] = 0; index[0] < npoints; ++index[0]) {
                CCTK_REAL & du = diff_grid_function[
                    index[0]*stride_diff_grid_function[0] +
                    index[1]*stride_diff_grid_function[1] +
                    index[2]*stride_diff_grid_function[2]];
                du = 0;

                // moving index
                int mindex[3] = {index[0], index[1], index[2]};
                int & i = mindex[dir];
                for(i = 0; i < npoints; ++i) {
                    du += grid_function[
                            mindex[0]*stride_grid_function[0] +
                            mindex[1]*stride_grid_function[1] +
                            mindex[2]*stride_grid_function[2]] *
                        diffop[index[dir]][i];
                }
                du = du / (0.5*shape[dir]);
            }
        }

        //! \brief evaluates the pseudo-spectral derivative of a grid
        //! function in a point
        /*!
         *  \tparam dir direction of the derivative
         *
         *  This method evaluates the differential of a function over
         *  a single element.
         *
         * In the reference element we have
         *  \f[
                [\partial_x u]_{ijk} = \sum_{q} u_{qjk}\ l'_q(x_i), \qquad
                [\partial_y u]_{ijk} = \sum_{q} u_{iqk}\ l'_q(y_j), \qquad
                [\partial_z u]_{ijk} = \sum_{q} u_{ijq}\ l'_q(z_k);
            \f]
         *  where \f$l_i\f$ is the one-dimensional collocation basis.
         *
         *  Input and outputs grid-arrays should be offsetted so that the 0th
         *  element represents the value of the grid function in the first
         *  point of the element (in the lexicographical ordering).
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically and divides the result by the mass matrix.
         *
         *  \see GNIGrid::restriction
         */
        template<policy::direction_t dir>
        static CCTK_REAL diff(
                //! [in] size of the element in each direction
                CCTK_REAL const shape[3],
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access grid_function array
                int const stride_grid_function[3],
                //! [in] local index where to compute the derivative
                int const index[3]) {
            CCTK_REAL out = 0;
            int mindex[3] = {index[0], index[1], index[2]};
            int & i = mindex[dir];

            for(i = 0; i < npoints; ++i) {
                out += grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]] *
                    diffop[index[dir]][i];
            }
            return out / (0.5*shape[dir]);
        }

        //! computes the weak derivative in a point
        /*!
         *  The weak derivative is computed as
         *  \f[
                D u = \partial_x u + \hat{u} j_u,
            \f]
         *  where \f$j_u\f$ is the jump set of \f$u\f$, i.e. the inter-element
         *  boundary, and \f$\hat{u}\f$ is the numerical flux
         *  \f[
               \hat{u} = 1/2 (u_L + u_R).
            \f]
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically and divides the result by the mass matrix.
         */
        template<policy::direction_t dir>
        static CCTK_REAL wdiff(
                //! [in] size of the element in each direction
                CCTK_REAL const shape[3],
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access grid_function array
                int const stride_grid_function[3],
                //! [in] local index where to compute the derivative
                int const index[3]) {
            CCTK_REAL out = 0;
            int mindex[3] = {index[0], index[1], index[2]};
            int & i = mindex[dir];

            if(i == 0) {
                // Left and right states
                i = -1;
                CCTK_REAL const uL = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];
                i = 0;
                CCTK_REAL const uR = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];

                // Centered flux
                CCTK_REAL const flux = 0.5*(uL + uR);
                // Lax-Friedrichs flux
                //CCTK_REAL const flux = 0.5*(uL + uR) - 0.5*s*(uR - uL);

                // Flux remainder (the interior value is already taken
                // care of by the pseudo-spectral derivative)
                out = - flux + uR;
            }
            if(i == npoints - 1) {
                // Left and right states
                CCTK_REAL const uL = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];
                i = npoints;
                CCTK_REAL const uR = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];

                // Centered flux
                CCTK_REAL const flux = 0.5*(uL + uR);
                // Lax-Friedrichs flux
                //CCTK_REAL const flux = 0.5*(uL + uR) - 0.5*s*(uR - uL);

                // Flux remainder (the interior value is already taken
                // care of by the pseudo-spectral derivative)
                out = flux - uL;
            }

            out = out/(0.5*shape[dir]*weight[index[dir]]) +
                diff<dir>(shape, grid_function, stride_grid_function, index);

            return out;
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
        static CCTK_REAL jump(
                //! [in] size of the element in each direction
                CCTK_REAL const shape[3],
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access grid_function array
                int const stride_grid_function[3],
                //! [in] local index where to compute the derivative
                int const index[3]) {
            CCTK_REAL out = 0;
            int mindex[3] = {index[0], index[1], index[2]};
            int & i = mindex[dir];

            if(i == 0) {
                // Left and right states
                i = -1;
                CCTK_REAL const uL = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];
                i = 0;
                CCTK_REAL const uR = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];

                out = uL - uR;
            }
            if(i == npoints - 1) {
                // Left and right states
                CCTK_REAL const uL = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];
                i = npoints;
                CCTK_REAL const uR = grid_function[
                        mindex[0]*stride_grid_function[0] +
                        mindex[1]*stride_grid_function[1] +
                        mindex[2]*stride_grid_function[2]];

                out = uL - uR;
            }

            out = out/(0.5*shape[dir]*weight[index[dir]]);

            return out;
        }

        //! integrates a grid function using a Gaussian quadrature formula
        /*!
         *  Input and outputs grid-arrays should be offsetted so that the 0th
         *  element represents the value of the grid function in the first
         *  point of the element (in the lexicographical ordering).
         *
         *  This function takes care of the mapping from-to the reference
         *  element automatically.
         *
         *  \see GNIGrid::restriction
         */
        static CCTK_REAL integrate(
                //! [in] size of the element in each direction
                CCTK_REAL const shape[3],
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access grid_function array
                int const stride_grid_function[3]) {
            CCTK_REAL s = 0;
            for(int k = 0; k < npoints; ++k)
            for(int j = 0; j < npoints; ++j)
            for(int i = 0; i < npoints; ++i) {
                s += grid_function[i*stride_grid_function[0] +
                    j*stride_grid_function[1] + k*stride_grid_function[2]] *
                    weight[i]*weight[j]*weight[k];
            }
            return s*0.5*shape[0]*0.5*shape[1]*0.5*shape[2];
        }

        //! filter a given grid function
        /*!
         *  Since the FE space is a tensor-product space we can perform the
         *  filtering in a dimension-by-dimension way:
         *  \f[
                [Fu]_{ijk} = \mathrm{IDLT}_{ia} \mathrm{IDLT}_{jb}
                    \mathrm{IDLT}_{kc} \sigma_a \sigma_b \sigma_c \hat{u}_{abc},
            \f]
         *  that is:
         *  \f[
                [Fu]_{ijk} =
                    \mathrm{IDLT}_{ia} \sigma_a \mathrm{DLT}_{al}
                    \mathrm{IDLT}_{jb} \sigma_b \mathrm{DLT}_{bm}
                    \mathrm{IDLT}_{kc} \sigma_c \mathrm{DLT}_{cn} u_{lmn}
            \f]
         */
        static void filter(
                //! [in] filter, must be an npoints vector
                CCTK_REAL const filter[npoints],
                //! [in, out] target *restricted* grid function
                CCTK_REAL * grid_function,
                //! [in] stride used to access the grid_function
                int const stride_grid_function[3]) {
            CCTK_REAL scratch[npoints];
            for(int dir = 0; dir < 3; ++dir) {
                // Direction map
                int diridx[3];
                switch(dir) {
                    case 0:
                        diridx[0] = 0;
                        diridx[1] = 1;
                        diridx[2] = 2;
                        break;
                    case 1:
                        diridx[0] = 1;
                        diridx[1] = 0;
                        diridx[2] = 2;
                        break;
                    case 2:
                        diridx[0] = 2;
                        diridx[1] = 0;
                        diridx[2] = 1;
                        break;
                }

                for(int k = 0; k < npoints; ++k)
                for(int j = 0; j < npoints; ++j) {
                    // compute the DLT of grid_function
                    for(int a = 0; a < npoints; ++a) {
                        scratch[a] = 0;
                        for(int i = 0; i < npoints; ++i) {
                            scratch[a] += dltop[a][i]*grid_function[
                                i*stride_grid_function[diridx[0]] +
                                j*stride_grid_function[diridx[1]] +
                                k*stride_grid_function[diridx[2]]];
                        }
                    }
                    // apply filter to the DLT
                    for(int a = 0; a < npoints; ++a) {
                        scratch[a] = scratch[a]*filter[a];
                    }
                    // compute the IDLT of scratch
                    for(int i = 0; i < npoints; ++i) {
                        CCTK_REAL & u = grid_function[
                            i*stride_grid_function[diridx[0]] +
                            j*stride_grid_function[diridx[1]] +
                            k*stride_grid_function[diridx[2]]];
                        u = 0;
                        for(int a = 0; a < npoints; ++a) {
                            u += idltop[i][a]*scratch[a];
                        }
                    }
                }
            }
        }

        //! computes a smoothness indicator of the given grid function
        /*!
         *  This is a measure of how smooth is a function in a scale
         *  from 0 to 1 (1 for smooth functions).
         *
         *  The idea is to use this measure to hybridize the numerical
         *  scheme that we use with some more dissipative scheme
         *
         *  Currently we only return 0 or 1, but this may change in
         *  the future.
         *
         *  NOTE: Currently this is not working very well
         *  TODO: This has to be improved
         */
        static CCTK_REAL smoothness_indicator(
                //! [in] target *restricted* grid function
                CCTK_REAL const * grid_function,
                //! [in] stride used to access the grid_function
                int const stride_grid_function[3]) {
            CCTK_REAL scratch[npoints];
            for(int dir = 0; dir < 3; ++dir) {
                // Direction map
                int diridx[3];
                switch(dir) {
                    case 0:
                        diridx[0] = 0;
                        diridx[1] = 1;
                        diridx[2] = 2;
                        break;
                    case 1:
                        diridx[0] = 1;
                        diridx[1] = 0;
                        diridx[2] = 2;
                        break;
                    case 2:
                        diridx[0] = 2;
                        diridx[1] = 0;
                        diridx[2] = 1;
                        break;
                }

                for(int k = 0; k < npoints; ++k)
                for(int j = 0; j < npoints; ++j) {
                    // compute the DLT of grid_function
                    for(int a = 0; a < npoints; ++a) {
                        scratch[a] = 0;
                        for(int i = 0; i < npoints; ++i) {
                            scratch[a] += dltop[a][i]*grid_function[
                                i*stride_grid_function[diridx[0]] +
                                j*stride_grid_function[diridx[1]] +
                                k*stride_grid_function[diridx[2]]];
                        }
                    }
                    // look at the decay of the coefficients
                    for(int a = 0; a < npoints; ++a) {
                        scratch[a] = std::abs(scratch[a]);
                    }
                    if(scratch[npoints-1] > std::max(scratch[0], scratch[1])) {
                        return 0;
                    }
                }
            }
            return 1.0;
        }

        //! Computes the prolongation of a grid function on a whole element
        /*!
         *  The prolongation is computed in a tensor product way
         *  \f[
                u^{abc} = P^c_k P^b_j P^a_i u^{ijk}
            \f]
         *  where \f$P\f$ is the prolongation matrix
         *
         *  \see GNIGrid::prolongation
         */
        static void prolongate_full(
                //! [in] target *restricted* grid function to prolongate
                CCTK_REAL const * coarse_grid_function,
                //! [in] stride of the coarse grid function array
                int const stride_coarse_grid_function[3],
                //! [out] target *restricted* prolongated grid function
                CCTK_REAL * prolonged_grid_function,
                //! [in] stride of the prolongated grid function array
                int const stride_prolonged_grid_function[3]) {
            CCTK_REAL scratch[npoints];

            // u^{ajk} = P^a_i u^{ijk}
            for(int k = 0; k <   npoints; ++k)
            for(int j = 0; j <   npoints; ++j)
            for(int a = 0; a < 2*npoints; ++a) {
                CCTK_REAL & Pu = prolonged_grid_function[
                      a*stride_prolonged_grid_function[0] +
                    2*j*stride_prolonged_grid_function[1] +
                    2*k*stride_prolonged_grid_function[2]];
                Pu = 0;
                for(int i = 0; i < npoints; ++i) {
                    Pu += prolongation[a][i]*coarse_grid_function[
                        i*stride_coarse_grid_function[0] +
                        j*stride_coarse_grid_function[1] +
                        k*stride_coarse_grid_function[2]];
                }
            }

            // u^{abk} = P^b_j u^{ajk}
            for(int k = 0; k <   npoints; ++k)
            for(int a = 0; a < 2*npoints; ++a) {
                for(int j = 0; j < npoints; ++j) {
                    scratch[j] = prolonged_grid_function[
                          a*stride_prolonged_grid_function[0] +
                        2*j*stride_prolonged_grid_function[1] +
                        2*k*stride_prolonged_grid_function[2]];
                }
                for(int b = 0; b < 2*npoints; ++b) {
                    CCTK_REAL & Pu = prolonged_grid_function[
                          a*stride_prolonged_grid_function[0] +
                          b*stride_prolonged_grid_function[1] +
                        2*k*stride_prolonged_grid_function[2]];
                    Pu = 0;
                    for(int j = 0; j < npoints; ++j) {
                        Pu += prolongation[b][j]*scratch[j];
                    }
                }
            }

            // u^{abc} = P^c_k u^{abk}
            for(int b = 0; b < 2*npoints; ++b)
            for(int a = 0; a < 2*npoints; ++a) {
                for(int k = 0; k < npoints; ++k) {
                    scratch[k] = prolonged_grid_function[
                          a*stride_prolonged_grid_function[0] +
                          b*stride_prolonged_grid_function[1] +
                        2*k*stride_prolonged_grid_function[2]];
                }
                for(int c = 0; c < 2*npoints; ++c) {
                    CCTK_REAL & Pu = prolonged_grid_function[
                        a*stride_prolonged_grid_function[0] +
                        b*stride_prolonged_grid_function[1] +
                        c*stride_prolonged_grid_function[2]];
                    Pu = 0;
                    for(int k = 0; k < npoints; ++k) {
                        Pu += prolongation[c][k]*scratch[k];
                    }
                }
            }
        }
        //! \brief Computes the prolongation of a grid function on a face
        //! of the element
        /*!
         *  The prolongation is computed in a tensor product way
         *  \f[
                u^{ab} = P^b_j P^a_i u^{ij}
            \f]
         *  where \f$P\f$ is the prolongation matrix
         *
         *  NOTE: The input and output arrays should be offsetted to point
         *  to the first point in the ``compactified'' element
         *
         *  \see GNIGrid::prolongation
         */
        static void prolongate_2D(
                //! [in] target *restricted* grid function to prolongate
                CCTK_REAL const * coarse_grid_function,
                //! [in] stride of the coarse grid function array
                int const stride_coarse_grid_function[2],
                //! [out] target *restricted* prolongated grid function
                CCTK_REAL * prolonged_grid_function,
                //! [in] stride of the prolongated grid function array
                int const stride_prolonged_grid_function[2]) {
            CCTK_REAL scratch[npoints];

            // u^{aj} = P^a_i u^{ij}
            for(int a = 0; a < 2*npoints; ++a)
            for(int j = 0; j <   npoints; ++j) {
                CCTK_REAL & Pu = prolonged_grid_function[
                      a*stride_prolonged_grid_function[0] +
                    2*j*stride_prolonged_grid_function[1]];
                Pu = 0;
                for(int i = 0; i < npoints; ++i) {
                    Pu += prolongation[a][i]*coarse_grid_function[
                        i*stride_coarse_grid_function[0] +
                        j*stride_coarse_grid_function[1]];
                }
            }

            // u^{ab} = P^b_j u^{aj}
            for(int a = 0; a < 2*npoints; ++a) {
                for(int j = 0; j < npoints; ++j) {
                    scratch[j] = prolonged_grid_function[
                          a*stride_prolonged_grid_function[0] +
                        2*j*stride_prolonged_grid_function[1]];
                }
                for(int b = 0; b < 2*npoints; ++b) {
                    CCTK_REAL & Pu = prolonged_grid_function[
                        a*stride_prolonged_grid_function[0] +
                        b*stride_prolonged_grid_function[1]];
                    Pu = 0;
                    for(int j = 0; j < npoints; ++j) {
                        Pu += prolongation[b][j]*scratch[j];
                    }
                }
            }
        }


        //! Computes the restriction of a grid function on a whole element
        /*!
         *  The restriction is computed in a tensor product way
         *  \f[
                u^{ijk} = R^i_a R^j_b R^k_c u^{abc}
            \f]
         *  where \f$R\f$ is the restriction operator
         *
         *  \see GNIGrid::restriction
         */
        static void restrict_full(
                //! [in] target *restricted* grid function to restrict
                CCTK_REAL const * grid_function,
                //! [in] stride of the grid function array
                int const stride_grid_function[3],
                //! [out] target *restricted* restricted grid function
                CCTK_REAL * restricted_grid_function,
                //! [in] stride of the restricted grid function array
                int const stride_restricted_grid_function[3]) {
            CCTK_REAL * rest = NULL;
            try {
                rest = new CCTK_REAL[utils::pow<3>(2*npoints)];
            }
            catch(std::bad_alloc & e) {
#pragma omp critical
                CCTK_WARN(CCTK_WARN_ABORT, "Out of memory!");
            }
            int const stride_rest[3] = {1, 2*npoints, 2*npoints*2*npoints};
            CCTK_REAL scratch[2*npoints];

            // Copy the grid function to the workspace
            for(int c = 0; c < 2*npoints; ++c)
            for(int b = 0; b < 2*npoints; ++b)
            for(int a = 0; a < 2*npoints; ++a) {
                rest[
                    a*stride_rest[0] +
                    b*stride_rest[1] +
                    c*stride_rest[2]] =
                        grid_function[
                            a*stride_grid_function[0] +
                            b*stride_grid_function[1] +
                            c*stride_grid_function[2]];
            }

            // u^{abk} = R^k_c u^{abc}
            for(int b = 0; b < 2*npoints; ++b)
            for(int a = 0; a < 2*npoints; ++a)
            for(int k = 0; k <   npoints; ++k) {
                CCTK_REAL & Ru = rest[
                      a*stride_rest[0] +
                      b*stride_rest[1] +
                    2*k*stride_rest[2]];
                Ru = 0;
                for(int c = 0; c < 2*npoints; ++c) {
                    Ru += restriction[k][c]*grid_function[
                        a*stride_grid_function[0] +
                        b*stride_grid_function[1] +
                        c*stride_grid_function[2]];
                }
            }

            // u^{ajk} = R^j_b u^{abk}
            for(int a = 0; a < 2*npoints; ++a)
            for(int k = 0; k <   npoints; ++k) {
                for(int b = 0; b < 2*npoints; ++b) {
                    scratch[b] = rest[
                          a*stride_rest[0] +
                          b*stride_rest[1] +
                        2*k*stride_rest[2]];
                }
                for(int j = 0; j < npoints; ++j) {
                    CCTK_REAL & Ru = rest[
                          a*stride_rest[0] +
                        2*j*stride_rest[1] +
                        2*k*stride_rest[2]];
                    Ru = 0;
                    for(int b = 0; b < 2*npoints; ++b) {
                        Ru += restriction[j][b]*scratch[b];
                    }
                }
            }

            // u^{ijk} = R^i_a u^{ajk}
            for(int k = 0; k < npoints; ++k)
            for(int j = 0; j < npoints; ++j) {
                for(int a = 0; a < 2*npoints; ++a) {
                    scratch[a] = rest[
                          a*stride_rest[0] +
                        2*j*stride_rest[1] +
                        2*k*stride_rest[2]];
                }
                for(int i = 0; i < npoints; ++i) {
                    CCTK_REAL & Ru = rest[
                        2*i*stride_rest[0] +
                        2*j*stride_rest[1] +
                        2*k*stride_rest[2]];
                    Ru = 0;
                    for(int a = 0; a < 2*npoints; ++a) {
                        Ru += restriction[i][a]*scratch[a];
                    }
                }
            }

            // Copy back the results
            for(int k = 0; k < npoints; ++k)
            for(int j = 0; j < npoints; ++j)
            for(int i = 0; i < npoints; ++i) {
                restricted_grid_function[
                    i*stride_restricted_grid_function[0] +
                    j*stride_restricted_grid_function[1] +
                    k*stride_restricted_grid_function[2]] =
                        rest[
                            2*i*stride_rest[0] +
                            2*j*stride_rest[1] +
                            2*k*stride_rest[2]];
            }

            delete[] rest;
        }
        //! \brief Computes the prolongation of a grid function on a face
        //! of the element
        /*!
         *  The restriction is computed in a tensor product way
         *  \f[
                u^{ij} = R^i_a R^j_b u^{ab}
            \f]
         *  where \f$R\f$ is the restriction operator
         *
         *  NOTE: The input and output arrays should be offsetted to point
         *  to the first point in the ``compactified'' element
         *
         *  \see GNIGrid::restriction
         */
        static void restrict_2D(
                //! [in] target *restricted* grid function to restrict
                CCTK_REAL const * grid_function,
                //! [in] stride of the grid function array
                int const stride_grid_function[2],
                //! [out] target *restricted* restricted grid function
                CCTK_REAL * restricted_grid_function,
                //! [in] stride of the restricted grid function array
                int const stride_restricted_grid_function[2]) {
            CCTK_REAL * rest;
            try {
                rest = new CCTK_REAL[utils::pow<3>(2*npoints)];
            }
            catch(std::bad_alloc & e) {
#pragma omp critical
                CCTK_WARN(CCTK_WARN_ABORT, "Out of memory!");
            }
            int const stride_rest[2] = {1, npoints};
            CCTK_REAL scratch[2*npoints];

            // Copy the grid function to the workspace
            for(int b = 0; b < 2*npoints; ++b)
            for(int a = 0; a < 2*npoints; ++a) {
                rest[a*stride_rest[0] + b*stride_rest[1]] =
                    grid_function[a*stride_grid_function[0] +
                        b*stride_grid_function[1]];
            }

            // u^{aj} = R^j_b u^{ab}
            for(int a = 0; a < 2*npoints; ++a)
            for(int j = 0; j <   npoints; ++j) {
                CCTK_REAL & Ru = rest[a*stride_rest[0] + 2*j*stride_rest[1]];
                Ru = 0;
                for(int b = 0; b < 2*npoints; ++b) {
                    Ru += restriction[j][b]*grid_function[
                        a*stride_grid_function[0] +
                        b*stride_grid_function[1]];
                }
            }

            // u^{ij} = R^i_a u^{aj}
            for(int j = 0; j < npoints; ++j) {
                for(int a = 0; a < 2*npoints; ++a) {
                    scratch[a] = rest[a*stride_rest[0] + 2*j*stride_rest[0]];
                }
                for(int i = 0; i < npoints; ++i) {
                    CCTK_REAL & Ru = rest[2*i*stride_rest[0] +
                        2*j*stride_rest[1]];
                    Ru = 0;
                    for(int a = 0; a < 2*npoints; ++a) {
                        Ru += restriction[i][a]*scratch[a];
                    }
                }
            }

            // Copy back the results
            for(int j = 0; j < npoints; ++j)
            for(int i = 0; i < npoints; ++i) {
                restricted_grid_function[
                    i*stride_restricted_grid_function[0] +
                    j*stride_restricted_grid_function[1]] =
                        rest[2*i*stride_rest[0] + 2*j*stride_rest[1]];
            }

            delete[] rest;
        }

        //! collocation points in \f$[-1,1]\f$
        static CCTK_REAL const node[npoints];
        //! integration weights in \f$[-1,1]\f$
        static CCTK_REAL const weight[npoints];
        //! discrete derivative operator (i.e. the pseudo-spectral derivative)
        /*!
         *  \f[
                u'(x_i) \approx \sum l_k'(x_i) u(x_k) =: \sum_k D_{ik} u(x_k).
            \f]
         */
        static CCTK_REAL const diffop[npoints][npoints];
        //! pre-computed interpolation coefficients
        /*!
         *  \f[
                q_i = \prod_{k\neq i} \frac{1}{x_i - x_k}.
            \f]
         *
         *  So that the Lagrange polynomials can be efficiently computed as
         *  \f[
               l_i(x) = \prod_{k\neq i} \frac{x - x_k}{x_i - x_k} =
                    q_i \prod_{k\neq i} [x - x_k].
            \f]
         */
        static CCTK_REAL const icoeff[npoints];
        //! pre-computed discrete-Legendre-transform matrix
        /*!
         *  \f[
                \hat{u}_{ijk} = \sum_{abc} \mathrm{DLT}_{ia}\ \mathrm{DLT}_{jb}\
                    \mathrm{DLT}_{kc}\ u_{abc}.
            \f]
         */
        static CCTK_REAL const dltop[npoints][npoints];
        //! pre-computed inverse-discrete-Legendre-transform matrix
        /*!
         *  \f[
                u_{ijk} = \sum_{abc} \mathrm{IDLT}_{ia}\ \mathrm{IDLT}_{jb}\
                    \mathrm{IDLT}_{kc}\ \hat{u}_{abc}.
            \f]
         */
        static CCTK_REAL const idltop[npoints][npoints];
        //! prolongation operator
        /*!
         *  Prolongation operator for fixed ratio, 2, AMR
         *  \f[
                u^a = P^a_i u^i = l_i(x^a) u^i
         *  \f]
         *  where \f$u^a\f$ are the collocation coefficients on the fine mesh
         *  and \f$u^i\f$ are the collocation coefficients on the coarse mesh.
         */
        static CCTK_REAL const prolongation[2*npoints][npoints];
        //! restriction operator
        /*!
         *  Restriction operator for fixed ratio, 2, AMR
         *  \f[
                u^i =  R^i_a u^a
         *  \f]
         *  where \f$u^a\f$ are the collocation coefficients on the fine mesh
         *  and \f$u^i\f$ are the collocation coefficients on the coarse mesh.
         *
         *  The restriction operator is computed as pseudo-inverse of the
         *  prolongation operator, i.e. so that \f$ R P = \mathrm{Id} \f$ in the
         *  least-square sense.
         */
        static CCTK_REAL const restriction[npoints][2*npoints];
};

} // namespace

#endif
