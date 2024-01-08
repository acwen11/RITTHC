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


#ifndef HRSCC_CLAW_HH
#define HRSCC_CLAW_HH

#include <cctk.h>

#include <hrscc_observer.hh>
#include <hrscc_traits.hh>
#include <hrscc_typedefs.hh>

namespace hrscc {

//! conservation law prototype
/*!
 *  This class provides the specifications for the physics driver. The final
 *  user of the HRSCCore thorn is supposed to create a class \e MyLaw
 *  \code
 *      class MyLaw: public hrscc::CLaw<MyLaw>;
 *  \endcode
 *  which overloads all the virtual methods present here.
 *
 *  The user should also provide the required traits of his/her class by
 *  specialising hrscc::traits
 *
 *  \tparam derived_t derived class (for static polymorphism)
 */
template<typename derived_t>
class CLaw {
    public:
        typedef derived_t derived;

        //! number of equations in our conservation law
        enum {nequations = traits<derived>::nequations};
        //! number of external fields directly needed during the evolution
        enum {nexternal  = traits<derived>::nexternal};
        //! total number of CCTK_REAL variables used for the evolution
        enum {nvariables = 3*nequations + nexternal};
        //! number of bitmask fields used by the physics driver
        enum {nbitmasks  = traits<derived>::nbitmasks};

        //! \brief if true this means that we have an exact conservation law,
        //! ie with no sources
        static bool const pure = traits<derived>::pure;

        //! Cactus variable indices of the conserved quantities
        static int conserved_idx[nequations];
        //! Cactus variable indices of the primitive quantities
        static int primitive_idx[nequations];
        //! Cactus variable indices of the RHS variables
        static int rhs_idx[nequations];
        //! Cactus variable indices for the extra fields
        static int field_idx[nexternal];
        //! Cactus variable indices of the bitmasks used by the physics driver
        static int bitmask_idx[nbitmasks];
        //! Cactus variable indices for the numerical fluxes in each direction
        /*!
         *  Each point stores the flux flowing to its left
         */
        static int num_flux_idx[3*nequations];

        //! Lower bound for the conserved variables (for the positivity limiter)
        static CCTK_REAL conserved_lbound[nequations];

        //! compute all the variables from the primitive ones in a point
        inline void prim_to_all(
                //! [in, out] local value of all the variables
                Observer<CLaw<derived> > & observer
                ) const {
            static_cast<derived const*>(this)->prim_to_all(observer);
        }
        //! compute all the variables from the conserved ones in a point
        inline void cons_to_all(
                //! [in, out] local value of all the variables
                Observer<CLaw<derived> > & observer
                ) const {
            static_cast<derived const*>(this)->cons_to_all(observer);
        }
        //! compute the fluxes in the given direction in a point
        template<policy::direction_t direction>
        inline void fluxes(
                //! [in, out] local value of all the variables
                Observer<CLaw<derived> > & observer
                ) const {
            static_cast<derived const*>(this)->
              template fluxes<direction>(observer);
        }
        //! compute the eigenvalues in the given direction in a point
        template<policy::direction_t direction>
        inline void eigenvalues(
                //! [in, out] local value of all the variables
                Observer<CLaw<derived> > & observer
                ) const {
            static_cast<derived const*>(this)->
              template eigenvalues<direction>(observer);
        }
        //! \brief compute both the eigenvalues and the (left and right)
        //! eigenvectors in a point
        template<policy::direction_t direction>
        inline void eig(
                //! [in, out] local value of all the variables
                Observer<CLaw<derived> > & observer
                ) const {
            static_cast<derived const*>(this)->
              template eig<direction>(observer);
        }
};

} // namespace

#endif
