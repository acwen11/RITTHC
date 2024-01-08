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


#ifndef HRSCC_OBSERVER_HH
#define HRSCC_OBSERVER_HH

#include <cctk.h>

namespace hrscc {

//! used to comunicate data to/from the physics driver
/*!
 *  This class stores a copy of the data at the current location.
 *
 *  Depending on the situation not all of the data might be meaningfull.
 */
template<typename claw_t>
class Observer {
    public:
        typedef claw_t claw;

        //! cactus grid hierarchy
        cGH const * const cctkGH;

        //! conserved variables at the point
        CCTK_REAL conserved[claw::nequations];
        //! primitive variables at the point
        CCTK_REAL primitive[claw::nequations];
        //! rhs at the point
        CCTK_REAL rhs[claw::nequations];
        //! external fields at the point
        CCTK_REAL field[claw::nexternal];

        //! bitmasks at the point
        CCTK_INT bitmask[claw::nbitmasks];

        //! array for the fluxes at the point
        CCTK_REAL flux[3][claw::nequations];
        //! array of eigenvalue at the point
        CCTK_REAL eigenvalue[claw::nequations];
        //! \brief array of left eigenvectors at the point stored in a
        //! row-major order
        CCTK_REAL left_eigenvector[claw::nequations][claw::nequations];
        //! \brief array of right eigenvectors at the point stored in a
        //! row-major order
        CCTK_REAL right_eigenvector[claw::nequations][claw::nequations];

        //! current location on the grid
        CCTK_REAL x, y, z;

        //! current position in the grid (local index)
        int i, j, k;
        int ijk;

        //! signals if the current position is shifted from the grid point
        bool shift[3];

        //! The constructor
        Observer(
                //! [in] pointer to the Cactus hierarchy
                cGH const * const cctkGH,
                //! [in] coordinates
                CCTK_REAL const * const coordinates[3],
                //! [in] real grid variables
                CCTK_REAL * const grid_variable[claw::nvariables],
                //! [in] bitmasks
                CCTK_INT * const bitmask[claw::nbitmasks]
                ): cctkGH(cctkGH), i(-1), j(-1), k(-1), ijk(-1) {
            for(int a = 0; a < 3; ++a) {
                _M_grid_coordinates[a] = coordinates[a];
            }

            for(int v = 0; v < claw::nvariables; ++v) {
                _M_grid_variable[v] = grid_variable[v];
            }
            for(int v = 0; v < claw::nbitmasks; ++v) {
                _M_grid_bitmask[v] = bitmask[v];
            }

            for(int a = 0; a < claw::nequations; ++a) {
                flux[0][a] = flux[1][a] = flux[2][a] = 0;
                eigenvalue[a] = 0;
                for(int b = 0; b < claw::nequations; ++b) {
                    left_eigenvector[a][b] = 0;
                    right_eigenvector[a][b] = 0;
                }
            }

            i = j = k = ijk = -1;
            shift[0] = shift[1] = shift[2] = false;
        }

        //! Moves the observer to a given grid point
        void jump_to_location(
                //! [in] x-index
                int i_,
                //! [in] y-index
                int j_,
                //! [in] z-index
                int k_
                ) {
            assert(i_ >= 0 && i_ < cctkGH->cctk_lsh[0]);
            assert(j_ >= 0 && j_ < cctkGH->cctk_lsh[1]);
            assert(k_ >= 0 && k_ < cctkGH->cctk_lsh[2]);

            i = i_;
            j = j_;
            k = k_;

            shift[0] = false;
            shift[1] = false;
            shift[2] = false;

            ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            x = _M_grid_coordinates[0][ijk];
            y = _M_grid_coordinates[1][ijk];
            z = _M_grid_coordinates[2][ijk];

            int idx = 0;
            for(int v = 0; v < claw::nequations; ++v) {
                conserved[v] = _M_grid_variable[idx][ijk];
                ++idx;
            }
            for(int v = 0; v < claw::nequations; ++v) {
                primitive[v] = _M_grid_variable[idx][ijk];
                ++idx;
            }
            for(int v = 0; v < claw::nequations; ++v) {
                rhs[v] = _M_grid_variable[idx][ijk];
                ++idx;
            }
            for(int v = 0; v < claw::nexternal; ++v) {
                field[v] = _M_grid_variable[idx][ijk];
                ++idx;
            }

            for(int v = 0; v < claw::nbitmasks; ++v) {
                bitmask[v] = _M_grid_bitmask[v][ijk];
            }
        }

        //! \brief jump to the given cell interface and reconstruct the fields
        //! with a linear interpolation
        /*!
         *  The bitmask at the interface is computed as the bitwise AND of
         *  the bitmask at the two neighbouring points
         */
        void linterp(
                //! [in] x-index
                int i_,
                //! [in] y-index
                int j_,
                //! [in] z-index
                int k_,
                //! [in] shift of half grid point in the x-direction
                bool shift_x,
                //! [in] shift of half grid point in the y-direction
                bool shift_y,
                //! [in] shift of half grid point in the z-direction
                bool shift_z
                ) {
            assert(i_ >= 0 && i_ < cctkGH->cctk_lsh[0]);
            assert(j_ >= 0 && j_ < cctkGH->cctk_lsh[1]);
            assert(k_ >= 0 && k_ < cctkGH->cctk_lsh[2]);
            assert(i_ + shift_x >= 0 && i_ + shift_x < cctkGH->cctk_lsh[0]);
            assert(j_ + shift_y >= 0 && j_ + shift_y < cctkGH->cctk_lsh[1]);
            assert(k_ + shift_z >= 0 && k_ + shift_z < cctkGH->cctk_lsh[2]);

            i = i_;
            j = j_;
            k = k_;

            shift[0] = shift_x;
            shift[1] = shift_y;
            shift[2] = shift_z;

            ijk = -1;
            int const ijk1 = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int const ijk2 = CCTK_GFINDEX3D(cctkGH, i+shift_x, j+shift_y,
                    k+shift_z);

            x = ldexp(_M_grid_coordinates[0][ijk1] +
                    _M_grid_coordinates[0][ijk2], -1);
            y = ldexp(_M_grid_coordinates[1][ijk1] +
                    _M_grid_coordinates[1][ijk2], -1);
            z = ldexp(_M_grid_coordinates[2][ijk1] +
                    _M_grid_coordinates[2][ijk2], -1);

            int idx = 0;
            for(int v = 0; v < claw::nequations; ++v) {
                conserved[v] = ldexp(_M_grid_variable[idx][ijk1] +
                        _M_grid_variable[idx][ijk2], -1);
                ++idx;
            }
            for(int v = 0; v < claw::nequations; ++v) {
                primitive[v] = ldexp(_M_grid_variable[idx][ijk1] +
                        _M_grid_variable[idx][ijk2], -1);
                ++idx;
            }
            for(int v = 0; v < claw::nequations; ++v) {
                rhs[v] = ldexp(_M_grid_variable[idx][ijk1] +
                        _M_grid_variable[idx][ijk2], -1);
                ++idx;
            }
            for(int v = 0; v < claw::nexternal; ++v) {
                field[v] = ldexp(_M_grid_variable[idx][ijk1] +
                        _M_grid_variable[idx][ijk2], -1);
                ++idx;
            }

            for(int v = 0; v < claw::nbitmasks; ++v) {
                bitmask[v] = _M_grid_bitmask[v][ijk1] &
                    _M_grid_bitmask[v][ijk2];
            }
        }

        //! writes values of the variables into the grid functions
        /*!
         *  NOTE: this does not make sense if the observer is not positioned
         *  on a grid point
         */
        void record() {
            assert(ijk >= 0);

            int idx = 0;
            for(int v = 0; v < claw::nequations; ++v) {
                _M_grid_variable[idx][ijk] = conserved[v];
                ++idx;
            }
            for(int v = 0; v < claw::nequations; ++v) {
                _M_grid_variable[idx][ijk] = primitive[v];
                ++idx;
            }
            for(int v = 0; v < claw::nequations; ++v) {
                _M_grid_variable[idx][ijk] = rhs[v];
                ++idx;
            }
            for(int v = 0; v < claw::nexternal; ++v) {
                _M_grid_variable[idx][ijk] = field[v];
                ++idx;
            }

            for(int v = 0; v < claw::nbitmasks; ++v) {
                _M_grid_bitmask[v][ijk] = bitmask[v];
            }
        }

        //! consistency check on the given eigenvectors
        CCTK_REAL check_eigenvectors() const {
            CCTK_REAL err = 0.0;
            CCTK_REAL LRab;

            for(int a = 0; a < claw::nequations; ++a) {
                for(int b = 0; b < claw::nequations; ++b) {
                    LRab = 0;
                    for(int c = 0; c < claw::nequations; ++c) {
                        LRab += left_eigenvector[a][c]*right_eigenvector[c][b];
                    }
                    err = std::max(err, std::abs(LRab -
                                static_cast<CCTK_REAL>(a == b)));
                }
            }
            return err;
        }

        //! reset the eigenvalues
        /*!
         *  Sets all the eigenvalues to zero
         *
         *  This is usefull in situations where the spectral decomposition of
         *  the system becomes singular (typically when unphysical values
         *  appear in the calculation)
         */
        void reset_eigenvalues() {
            for(int a = 0; a < claw::nequations; ++a) {
                eigenvalue[a] = 0;
            }
        }

        //! reset the eigenvectors
        /*!
         *  Overwrite the eigenvectors with the identity matrix.
         *
         *  This is usefull in situations where the spectral decomposition of
         *  the system becomes singular (typically when unphysical values
         *  appear in the calculation)
         */
        void reset_eigenvectors() {
            for(int a = 0; a < claw::nequations; ++a) {
                for(int b = 0; b < claw::nequations; ++b) {
                    left_eigenvector[a][b]  = static_cast<CCTK_REAL>(a == b);
                    right_eigenvector[a][b] = static_cast<CCTK_REAL>(a == b);
                }
            }
        }
    private:
        CCTK_REAL const * _M_grid_coordinates[3];
        CCTK_REAL * _M_grid_variable[claw::nvariables];
        CCTK_INT  * _M_grid_bitmask[claw::nbitmasks];
};

} // namespace

#endif
