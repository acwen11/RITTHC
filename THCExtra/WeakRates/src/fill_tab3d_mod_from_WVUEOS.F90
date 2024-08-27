#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
 

! #########################################################
! TABLE UNITS:
!        density              g/cm^3
!        temperature          MeV
!        ye                   number fraction per baryon
! #########################################################

INTEGER function fill_tab3d_mod_from_WVUEOS()
  ! This routine copies EOS table constants from WVU_EOS.
  ! Some quantities in the module are no longer used as they are only needed for
  ! the internal weak table reader.

  use table3d_mod

  implicit none

  call WVU_EOS_share_nuc_eos(nrho, ntemp, nye,&
                              eos_rhomin, eos_rhomax, eos_tempmin, eos_tempmax,&
                              eos_yemin, eos_yemax, dlrho, dltemp, dye)

  ! WeakRates uses log base 10
  eos_lrhomin = eos_rhomin * log10(exp(1.0))
  eos_lrhomax = eos_rhomax * log10(exp(1.0))
  eos_ltempmin = eos_tempmin * log10(exp(1.0))
  eos_ltempmax = eos_tempmax * log10(exp(1.0))
  dlrho = dlrho * log10(exp(1.0))
  dltemp = dltemp * log10(exp(1.0))

  ! Hard coding mineps = 0 adjusted value!!!
  mass_fact = 930.17637269

  fill_tab3d_mod_from_WVUEOS = 0
end function fill_tab3d_mod_from_WVUEOS
