!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module jules43_dasoilm_Mod
!BOP
!
! !MODULE: jules43_dasoilm_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! Sujay Kumar; Initial Code
! 16Feb2017: Mahdi Navari; Modified for JULES 4.3 
!
! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: jules43_dasoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: jules43_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer              	:: nbins
     integer             	:: ntimes
     integer           		:: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: jules43_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: jules43_dasoilm_init
! \label{jules43_dasoilm_init}
! 
! !INTERFACE:
  subroutine jules43_dasoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
    integer                :: n 
    character*100          :: modelcdffile(LIS_rc%nnest)
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(jules43_dasm_struc)) then 
       allocate(jules43_dasm_struc(LIS_rc%nnest))
    endif
    
!TBD: SVK
#if 0 
    if(LIS_rc%dascaloption(k).eq."Linear scaling") then 
       call ESMF_ConfigFindLabel(LIS_config,"JULES.4.3 soil moisture CDF file:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'JULES.4.3 soil moisture CDF file: not defined')
       enddo
       
       do n=1,LIS_rc%nnest
       
!Hardcoded for now.
          jules43_dasm_struc(n)%nbins = 100
          
          call LIS_getCDFattributes(modelcdffile(n),&
               jules43_dasm_struc(n)%ntimes, ngrid)
          
          allocate(jules43_dasm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), jules43_dasm_struc(n)%ntimes, &
               jules43_dasm_struc(n)%nbins))
          allocate(jules43_dasm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), jules43_dasm_struc(n)%ntimes, &
               jules43_dasm_struc(n)%nbins))
          
          call LIS_readCDFdata(n,&
               jules43_dasm_struc(n)%nbins, &
               jules43_dasm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               jules43_dasm_struc(n)%model_xrange,&
               jules43_dasm_struc(n)%model_cdf)
       enddo
    endif
#endif

  end subroutine jules43_dasoilm_init
end module jules43_dasoilm_Mod
