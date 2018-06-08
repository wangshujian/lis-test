!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: jules43_qc_soilmobs
! \label{jules43_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 16Feb2017: Mahdi Navari; Modified for JULES 4.3 
!
! !INTERFACE:
subroutine jules43_qc_soilmobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use jules43_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil moisture observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen or (3) ground is fully or partially covered 
!  with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_sm_field

  real, pointer            :: smobs(:)
  integer                  :: t,j,pft
  integer                  :: gid
  integer                  :: status
  real                     :: lat,lon
  real                     :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! total soil moisture
  real                     :: smc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc4(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: sthu1(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! unfrozen soil moisture
  real                     :: sthu2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthu3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthu4(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: sthf1(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! frozen soil moisture
  real                     :: sthf2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthf3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sthf4(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real			   :: t_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: frac_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN  
  real 		 	   :: sliq_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: nsnow_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN 
  real 		 	   :: p_s_smvcwt_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN
  real 		 	   :: p_s_smvcst_tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)) ! MN

  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: sneqv_obs(LIS_rc%obs_ngrid(k))
  real                     :: nsnow_obs(LIS_rc%obs_ngrid(k))!sca_obs
  real                     :: shdfac_obs(LIS_rc%obs_ngrid(k))
  real                     :: t1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcwlt_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcmax_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthu1_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthu2_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthu3_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthu4_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf1_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf2_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf3_obs(LIS_rc%obs_ngrid(k))
  real                     :: sthf4_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))

  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in jules43_qc_soilmobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status) ! [m3/m3] ?
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in jules43_qc_soilmobs")

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
!     In Jules unfrozn soil moisture content are stored as a independent variable
!     so we can use that variable instead of computing that using smc and sthu  
!     total soil moisture ?
#if 0 
     smc1(t) = jules43_struc(n)%jules43(t)%smcl(1) ! [kg/m2]
     smc2(t) = jules43_struc(n)%jules43(t)%smcl(2)
     smc3(t) = jules43_struc(n)%jules43(t)%smcl(3)
     smc4(t) = jules43_struc(n)%jules43(t)%smcl(4)
     ! Unfrozen soil moisture
     sthu1(t) = jules43_struc(n)%jules43(t)%p_s_sthu(1) ![-]    Noah -->sh2o
     sthu2(t) = jules43_struc(n)%jules43(t)%p_s_sthu(2)
     sthu3(t) = jules43_struc(n)%jules43(t)%p_s_sthu(3)
     sthu4(t) = jules43_struc(n)%jules43(t)%p_s_sthu(4)
#endif

     ! Frozen soil moisture
     sthf1(t) = jules43_struc(n)%jules43(t)%p_s_sthf(1) ![-]    
     sthf2(t) = jules43_struc(n)%jules43(t)%p_s_sthf(2)
     sthf3(t) = jules43_struc(n)%jules43(t)%p_s_sthf(3)
     sthf4(t) = jules43_struc(n)%jules43(t)%p_s_sthf(4)

     stc1(t) = jules43_struc(n)%jules43(t)%t_soil(1) ! stc
     stc2(t) = jules43_struc(n)%jules43(t)%t_soil(2)
     stc3(t) = jules43_struc(n)%jules43(t)%t_soil(3)
     stc4(t) = jules43_struc(n)%jules43(t)%t_soil(4)

!     vegt(t) = jules43_struc(n)%jules43(t)%vegt  
  enddo

! MN ? considering using "LIS_convertPatchSpaceToObsSpace" is this loop correct?  
   do j=1, LIS_domain(n)%grid(k)%ntiles
      t = LIS_domain(n)%grid(k)%subgrid_tiles(j)
      pft= int(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt)
      vegt(t) = pft  
      t_tmp(t)  = jules43_struc(n)%jules43(t)%tstar_tile(pft)
      frac_tmp(t)= jules43_struc(n)%jules43(t)%frac(pft)
      sliq_tmp(t)= sum(jules43_struc(n)%jules43(t)%sliq(pft,:))
      nsnow_tmp(t)= jules43_struc(n)%jules43(t)%nsnow(pft)
      p_s_smvcwt_tmp(t)= jules43_struc(n)%jules43(t)%p_s_smvcwt(1)
      p_s_smvcst_tmp(t)= jules43_struc(n)%jules43(t)%p_s_smvcst(1)
   enddo


  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       jules43_struc(n)%jules43(:)%rainf,&
       rainf_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sliq_tmp,&
       sneqv_obs)  ! sneqv
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       nsnow_tmp ,&
       nsnow_obs) ! noah33_struc(n)%noah(:)%sca 
	! ????????????????????
	! MN: jules does not store the snow cover fraction. We use this 
	! variable to determine presence of snow on ground. I think  
	! we can use any variable that shows the presence of snow on ground
	! here  number of snow layer was replaced with "sca"
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
        frac_tmp,&
       shdfac_obs) ! noah33_struc(n)%noah(:)%shdfac
	! ????????????????????
	! I think we can also use LIS_surface(n,LIS_rc%lsm_index)%tile(t)%fgrd
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       t_tmp,&
       t1_obs) ! noah33_struc(n)%noah(:)%t1
! tstar_tile is a function of pft (plant functional types)   
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       p_s_smvcst_tmp,&
       smcmax_obs) !smcmax Here I assumed p_s_smvcst for all layers are the same
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       p_s_smvcwt_tmp,&
       smcwlt_obs) !smcwlt  Here I assumed p_s_smvcwt for all layers are the same
! ? why this "jules43_struc(n)%jules43(:)%p_s_smvcwt(1)" does not work
#if 0 
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       smc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc2,&
       smc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc3,&
       smc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc4,&
       smc4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu1,&
       sthu1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu2,&
       sthu2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu3,&
       sthu3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthu4,&
       sthu4_obs)
#endif 

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf1,&
       sthf1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf2,&
       sthf2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf3,&
       sthf3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sthf4,&
       sthf4_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc1,&
       stc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc2,&
       stc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc3,&
       stc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc4,&
       stc4_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       vegt,&
       vegt_obs)


  do t = 1,LIS_rc%obs_ngrid(k)
     if(smobs(t).ne.LIS_rc%udef) then 
        if(rainf_obs(t).gt.3E-6) then 
           smobs(t) = LIS_rc%udef
! use sthf(1, :) Frozen soil moisture content of the layers as a 
! fraction of saturation.
!        elseif(abs(smc1_obs(t)- &
!             sthu1_obs(t)).gt.0.0001) then
!           smobs(t) = LIS_rc%udef
!        elseif(abs(smc2_obs(t)- &
!             sthu2_obs(t)).gt.0.0001) then
!           smobs(t) = LIS_rc%udef
!        elseif(abs(smc3_obs(t)- &
!             sthu3_obs(t)).gt.0.0001) then
!           smobs(t) = LIS_rc%udef
!        elseif(abs(smc4_obs(t)- &
!             sthu4_obs(t)).gt.0.0001) then
!           smobs(t) = LIS_rc%udef
        elseif(sthf1_obs(t).gt.0.0001) then 
           smobs(t) = LIS_rc%udef! check the 0.0001 againt wilting point should be smaller than wilting point
        elseif(sthf2_obs(t).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(sthf3_obs(t).gt.0.0001) then
           smobs(t) = LIS_rc%udef
        elseif(sthf4_obs(t).gt.0.0001) then
           smobs(t) = LIS_rc%udef
! layer temperature --> jules43_struc(n)%jules43(t)%t_soil(i) 
        elseif(stc1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc2_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc3_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
        elseif(stc4_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef
! skin temperature 
! what is Tile surface temperatures (K) jules43_struc(n)%jules43(t)%tstar_tile(pft)
! is that skin temperature
        elseif(t1_obs(t).le.LIS_CONST_TKFRZ) then
           smobs(t) = LIS_rc%udef

!I got this from namelist (jules_surface_types.nml)
!ice=9,
!lake=7,
!nnvg=4,
!npft=5,
!soil=8,
!urban=6,)
! I got this from website (comparing the number with the following surface type makes sense. 
!Five Plant Functional Types (PFTs)
!Broadleaf trees
!Needle leaf trees
!C3 (temperate) grass
!C4 (tropical) grass
!Shrubs
!Four non-vegetation types
!Urban
!Inland water
!Bare soil
!Land-ice
!for Noah we just compare against forest hear we have to consider the pft of 1,2,6,7,9 ?

        elseif(vegt_obs(t).eq.1) then !Broadleaf trees
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).eq.2) then !Needle leaf trees
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).eq.6) then !Urban
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).eq.7) then !Inland water
           smobs(t) = LIS_rc%udef
        elseif(vegt_obs(t).eq.9) then !Land-ice
           smobs(t) = LIS_rc%udef



        elseif(sneqv_obs(t).gt.0.001) then 
           smobs(t) = LIS_rc%udef
        elseif(nsnow_obs(t).gt.0) then   ! MN number of snow layer  
           smobs(t) = LIS_rc%udef
        elseif(shdfac_obs(t).gt.0.7) then 
           smobs(t) = LIS_rc%udef        
!too close to the tails, could be due to scaling, so reject. 
        elseif(smcmax_obs(t)-smobs(t).lt.0.02) then ! MN: smbos comes from PILDAS [m3/m3]
           smobs(t) = LIS_rc%udef
        elseif(smobs(t) - smcwlt_obs(t).lt.0.02) then 
           smobs(t) = LIS_rc%udef
        endif
! apply the unit conversion. ! In the EnKF equation X(+) = X(-) + K(y-Hx), Units of the terms y and ! Hx must be the same. we just need to change the unit of y.  Unit of Hx ! term from JULES is [kg/m2] therefore we need to change the unit of y to! be consistent with Hx. 
!Question?!for the unit conversion FROM [M3/M3] to [kg/m2] we need layer thickness.  
!Here I have used the thickness of the first layer. (which is 10 cm in both PILDAS and JULES)    
       smobs = smobs / (1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000)   ! [m3w/m3s] / [1/m1s][m3w/kg] --> [kg/m2s]

     endif
  enddo

end subroutine jules43_qc_soilmobs

