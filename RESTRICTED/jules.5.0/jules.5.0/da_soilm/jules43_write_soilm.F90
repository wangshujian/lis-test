!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: jules43_write_soilm
! \label{jules43_write_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 28Feb2017: Mahdi Navari; Modified for JULES 4.3 
!
! !INTERFACE:
subroutine jules43_write_soilm(ftn,n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use jules43_lsmMod
  use LIS_historyMod, only : LIS_writevar_restart
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  integer :: t
  real, allocatable    :: tmp(:)

  allocate(tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules43_struc(n)%jules43(t)%smcl(1)
  enddo  
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules43_struc(n)%jules43(t)%smcl(2)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules43_struc(n)%jules43(t)%smcl(3)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = jules43_struc(n)%jules43(t)%smcl(4)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)
  deallocate(tmp)

end subroutine jules43_write_soilm

