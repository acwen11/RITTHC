#include "cctk.h"
#include "cctk_Parameters.h"

INTEGER FUNCTION enforceTableBounds(rho, temp, ye)
    use table3d_mod
    implicit none
    real*8, intent(inout) :: rho, temp, ye
    DECLARE_CCTK_PARAMETERS

    enforceTableBounds = 0

    ! TODO: ACTUALLY CHECK!!!
    ! rho is in cgs, temp is in MeV
    if(rho.lt.0.00010486587514604532) then
       rho = 0.00010486587514604532
       if(use_rho_min_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif
    if(rho.gt.1.0486587514603878e16) then
       rho = 1.0486587514603878e16
       if(use_rho_max_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif

    if(temp.lt.0.0008557232527644771) then
       temp = 0.0008557232527644771
       if(use_temp_min_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif
    if(temp.gt.250.6109253032114) then
       temp = 250.6109253032114
       if(use_temp_max_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif

    if(ye.lt.0.05) then
       ye = 0.05
       if(use_ye_min_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif
    if(ye.gt.0.56) then
       ye = 0.56
       if(use_ye_max_ext.eq.0) then
          enforceTableBounds = -1
       endif
    endif

    return
end function enforceTableBounds
