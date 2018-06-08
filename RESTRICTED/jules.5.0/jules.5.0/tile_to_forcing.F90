subroutine tile_to_forcing(n, t)
  use jules50_lsmMod
  use forcing
  use LIS_coreMod, only : LIS_rc
  use LIS_FORC_AttributesMod  
  use LIS_constantsMod, only : LIS_CONST_TKFRZ  
  use debug_latlon 
  use LIS_coreMod, only : LIS_rc
  implicit none 
  integer, intent(in) :: n, t

  qw_1_ij(1,1) = jules50_struc(n)%jules50(t)%qair/jules50_struc(n)%forc_count                       !  Total water content (Kg/Kg)
  tl_1_ij(1,1) = jules50_struc(n)%jules50(t)%tair/jules50_struc(n)%forc_count                       ! Ice/liquid water temperature (k)
  u_1_ij(1,1)  = jules50_struc(n)%jules50(t)%wind_e/jules50_struc(n)%forc_count                     ! W'ly wind component (m s-1)
  v_1_ij(1,1)  = jules50_struc(n)%jules50(t)%wind_n/jules50_struc(n)%forc_count                     ! S'ly wind component (m s-1)
  pstar_ij(1,1)   = jules50_struc(n)%jules50(t)%psurf/jules50_struc(n)%forc_count                   ! Surface pressure (Pascals)
  u_0_ij(1,1)  = 0.0                                                    ! W'ly component of surface current (m s-1)
  v_0_ij(1,1)  = 0.0                                                    ! S'ly component of surface current (m s-1)
  
  if(jules50_struc(n)%jules50(t)%rainf .ne. LIS_rc%udef) then
    ls_rain_ij(1,1) = jules50_struc(n)%jules50(t)%rainf/jules50_struc(n)%forc_count                 ! Large-scale rain (kg m-2/s)
  else
    ls_rain_ij(1,1) = 0.0
  endif

  if ( LIS_FORC_CRainf%selectOpt .eq. 1) then 
    if(jules50_struc(n)%jules50(t)%rainf_c .ne. LIS_rc%udef) then
      con_rain_ij(1,1) = jules50_struc(n)%jules50(t)%rainf_c/jules50_struc(n)%forc_count            ! Convective rain (kg m-2/s)
    else
      con_rain_ij(1,1) = 0.0
    endif
  else
    con_rain_ij(1,1) = 0.0
  endif

  if( LIS_FORC_Snowf%selectOpt .eq. 1) then
    if(jules50_struc(n)%jules50(t)%snowf .ne. LIS_rc%udef) then
      ls_snow_ij(1,1) = jules50_struc(n)%jules50(t)%snowf/jules50_struc(n)%forc_count               ! Large-scale snowfall (kg m-2/s)
      ! ls_rain_ij(1,1) = 0.0 
    else
      ls_snow_ij(1,1) = 0.0
    endif
    con_snow_ij(1,1) = 0.0 ! Convective snowfall (kg m-2/s)
  else
    if (jules50_struc(n)%jules50(t)%tair < LIS_CONST_TKFRZ) then  
      ls_snow_ij(1,1) = ls_rain_ij(1,1)
      ls_rain_ij(1,1) = 0
      ! con_snow_ij(1,1) = con_rain(t,1)
      con_snow_ij(1,1) = 0
    else
      ls_snow_ij(1,1) = 0.0
      con_snow_ij(1,1) = 0.0 
    endif
  endif  

  sw_down_ij(1,1) = jules50_struc(n)%jules50(t)%swdown/jules50_struc(n)%forc_count                        ! Surface downward SW radiation (W m-2)
  lw_down_ij(1,1) = jules50_struc(n)%jules50(t)%lwdown/jules50_struc(n)%forc_count                        ! Surface downward LW radiation (W m-2)                
  !diurnal_temperature_range(1,1)                      ! diurnal temperature range (K), used when l_dailydisagg=T
  !diff_rad(1,1)       ! Input diffuse radiation (W m-2)
end subroutine tile_to_forcing

