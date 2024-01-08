//  AdvectHRSC: solves the advection equation using HRSCCore and Cactus
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


#ifndef ADV_CLAW_HH
#define ADV_CLAW_HH

#include "hrscc.hh"

#include "adv_traits.hh"

class LinearAdvection: public hrscc::CLaw<LinearAdvection> {
    public:
        typedef hrscc::CLaw<LinearAdvection> claw;

        LinearAdvection();

        inline void prim_to_all(
                hrscc::Observer<claw> & // observer
                ) const {}

        template<hrscc::policy::direction_t dir>
        inline void fluxes(
                hrscc::Observer<claw> & observer
                ) const {
            observer.flux[dir][0] = observer.field[dir]*observer.primitive[0];
        }

        template<hrscc::policy::direction_t dir>
        inline void eigenvalues(
                hrscc::Observer<claw> & observer
                ) const {
            observer.eigenvalue[0] = observer.field[dir];
        }

        template<hrscc::policy::direction_t dir>
        inline void eig(
                hrscc::Observer<claw> & observer
                ) const {
            observer.eigenvalue[0] = observer.field[dir];
            observer.left_eigenvector[0][0] = 1.0;
            observer.right_eigenvector[0][0] = 1.0;
        }
};

#endif
