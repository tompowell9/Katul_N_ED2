Module phenology_jeong
  implicit none
  
  real, parameter :: a = 65.9
  real, parameter :: b = 443.7
  real, parameter :: c = -0.034
  real, parameter :: d = 3.5

  real, parameter :: flush_timescale = 14.  ! days
  real, parameter :: color_timescale = 7.  ! days
  real  :: N_want
  real  :: N_have
  real  :: N_for_leaf

Contains
  
  subroutine update_jeong_phenology(month, year, doy, cpoly, isi, lat)

    use ed_state_vars, only: polygontype, sitetype, patchtype
    use pft_coms, only: phenology, c2n_leaf, q, qsw, l2n_stem, c2n_stem, &
         c2n_storage, C_resorption_factor, N_resorption_factor
    use decomp_coms, only: f_labile, f_fast
    use ed_therm_lib, only: calc_veg_hcap, update_veg_energy_cweh
    use ed_max_dims, only: n_pft
    use allometry, only: area_indices, ed_biomass, dbh2bl
    use consts_coms, only: pio180, pi1
	use nutrient_constants, only: nstorage_max_factor
!    use damm_model, only: activation_energy, activation_energy_max, ae_phen_amp

    implicit none

    type(polygontype), target :: cpoly
    real, intent(in) :: lat
    integer, intent(in) :: doy, isi, month, year
    type(sitetype), pointer :: csite
    integer :: ipa, ipft
    type(patchtype), pointer :: cpatch
    logical, dimension(n_pft) :: leaf_out_cold
    real :: budburst_threshold, color_threshold
    integer :: ico
    real :: bl_max, delta_bleaf, old_hcapveg, old_leaf_hcap, old_wood_hcap
    integer :: doy_use, iday
    real :: daylength, daylength_solstice, acosarg, solar_declination
    real :: day_angle, daylength_winter_solstice, oak_lai, pine_lai
    real, dimension(n_pft) :: max_glf, min_glf

    ! Current declination angle
    print*,'updating phenology'
    doy_use=doy
    if(mod(year,4) == 0 .and. doy > 59)doy_use = doy-1
    day_angle = 2. * pi1 * (doy_use - 1) / 365.
    solar_declination = 0.006918 - 0.399912 * cos(day_angle) + &
         0.070257 * sin(day_angle) - 0.006758*cos(2.*day_angle) +  &
         0.000907 * sin(2.*day_angle) - 0.002697 * cos(3.*day_angle) + &
         0.001480 * sin(3.*day_angle)
    acosarg = -tan(lat*pio180)*tan(solar_declination)
    if(acosarg < -1.)acosarg = -1.
    if(acosarg > 1.)acosarg = 1.
    daylength = 24./pi1 * acos(acosarg)

    ! Declination angle at solstice
    doy_use=172
    day_angle = 2. * pi1 * (doy_use - 1) / 365.
    solar_declination = 0.006918 - 0.399912 * cos(day_angle) + &
         0.070257 * sin(day_angle) - 0.006758*cos(2.*day_angle) +  &
         0.000907 * sin(2.*day_angle) - 0.002697 * cos(3.*day_angle) + &
         0.001480 * sin(3.*day_angle)
    acosarg = -tan(lat*pio180)*tan(solar_declination)
    if(acosarg < -1.)acosarg = -1.
    if(acosarg > 1.)acosarg = 1.
    daylength_solstice = 24./pi1 * acos(acosarg)

    ! Loop over patches.
    csite => cpoly%site(isi)
    do ipa = 1, csite%npatches
       cpatch => csite%patch(ipa)

       ! First, update degree days for spring phenology.
       if(month < 12)then
          if(csite%avg_daily_temp(ipa) > 278.15)then
             csite%sum_dgd(ipa) = csite%sum_dgd(ipa) +   &
                  (csite%avg_daily_temp(ipa)-278.15)
          else
             csite%sum_chd(ipa) = csite%sum_chd(ipa) + 1.0
          endif
       else
          csite%sum_dgd(ipa) = 0.0
          csite%sum_chd(ipa) = 0.0
       endif

       ! Update degree days for fall phenology.
       if(daylength < 15.)then
          csite%sum_cdd(ipa) = csite%sum_cdd(ipa) + &
               min(csite%avg_daily_temp(ipa)-283.15,0.)
       else
          csite%sum_cdd(ipa) = 0.0
       endif

       ! Reset litter inputs.
       csite%fsc_in(ipa) = 0.0
       csite%fsn_in(ipa) = 0.0
       csite%ssc_in(ipa) = 0.0
       csite%ssl_in(ipa) = 0.0
       csite%slsn_in(ipa) = 0.0!ATT
       csite%slsc_in(ipa) = 0.0!ATT
       
       !loop through and check if evergreen, if evergreen set min_glf to 1.0
       do ipft = 1, n_pft
          if (phenology(ipft) == 0)then
             min_glf(ipft) = 1.0
          else 
             min_glf(ipft) = 0.0
          endif
       enddo
       max_glf(:) = 1.0
       ! Determine if there is budburst.
       do ipft = 1, n_pft
          leaf_out_cold(ipft) = .false.
          ! Check if it is currently dormant.
          if(csite%phenophase(ipft,ipa) == 1)then!phenophase 1 = dormant
             ! Compute threshold.
             budburst_threshold = a + b * exp(c*csite%sum_chd(ipa)) +   &
                  d * (-10.0 + 5.5)
             if(csite%sum_dgd(ipa) >= budburst_threshold)then
                ! If threshold is exceeded, set flags.
                leaf_out_cold(ipft) = .true.
                csite%phenophase(ipft,ipa) = 2 !phenophase 2 = leaves flushing
                csite%green_leaf_factor(ipft,ipa) = min_glf(ipft)
             endif
          endif
       enddo
       ! If leaves are flushing, update green_leaf_factor.
       do ipft = 1, n_pft
          if(csite%phenophase(ipft,ipa) == 2)then
             csite%sum_cdd(ipa) = 0.0
             ! Update factor.
             if(phenology(ipft) == 0)then
                csite%green_leaf_factor(ipft,ipa) =   &
                     csite%green_leaf_factor(ipft,ipa) + 1./flush_timescale*0.25
             else
                csite%green_leaf_factor(ipft,ipa) =   &
                     csite%green_leaf_factor(ipft,ipa) + 1./flush_timescale
             endif
             ! If green_leaf_factor is at maximum, set new phenophase.
             if(csite%green_leaf_factor(ipft,ipa) >= max_glf(ipft))then
                csite%green_leaf_factor(ipft,ipa) = max_glf(ipft)
                csite%phenophase(ipft,ipa) = 3 !phenophase 3 = full leav cover
             endif
          endif
       enddo
        ! Determine if there is leaf coloration.
       write(*,*) "Our sum_cdd right now is: ", csite%sum_cdd(ipa)
       do ipft = 1, n_pft
          ! Check if it is currently green.
          if(csite%phenophase(ipft,ipa) == 3)then
             ! Compute threshold.
             color_threshold = -25.
             csite%green_leaf_factor(ipft,ipa) = max_glf(ipft)
             if(csite%sum_cdd(ipa) <= color_threshold)then
                ! If threshold is exceeded, set flags.
                if (phenology(ipft) /= 0)then
                   csite%phenophase(ipft,ipa) = 4
                endif
             endif
          endif
       enddo
       ! If leaves are coloring, update green_leaf_factor.
       do ipft = 1, n_pft
          if(csite%phenophase(ipft,ipa) == 4)then !phenophse 4 = leaves are coloring
             csite%sum_dgd(ipa) = 0.0
             csite%sum_chd(ipa) = 0.0
             ! Update factor.
             if(phenology(ipft) == 0)then
                csite%green_leaf_factor(ipft,ipa) =   &
                     csite%green_leaf_factor(ipft,ipa) - 1./color_timescale * 0.25
             else
                csite%green_leaf_factor(ipft,ipa) =   &
                     csite%green_leaf_factor(ipft,ipa) - 1./color_timescale
             endif
             ! If green_leaf_factor is at zero, set new phenophase.
             if(csite%green_leaf_factor(ipft,ipa) <= min_glf(ipft))then
                csite%green_leaf_factor(ipft,ipa) = min_glf(ipft)
                csite%phenophase(ipft,ipa) = 1
             endif
          endif
       enddo
       ! Update individual cohorts. phenophase 1=bare; 2=flushing; 3= full cover 4=coloring

       do ico = 1, cpatch%ncohorts
          cpatch%leaf_drop(ico) = 0.0
          ipft = cpatch%pft(ico)
          
          delta_bleaf = -9.9e9
          !dropping leaves but not resorbed yet
          if(csite%phenophase(ipft,ipa) == 1)then
             if(phenology(ipft) == 0)then
                bl_max = csite%green_leaf_factor(ipft,ipa) *   &
                     dbh2bl(cpatch%dbh(ico),ipft)
                delta_bleaf = cpatch%bleaf(ico) - bl_max
          !      cpatch%bleaf(ico) = min(cpatch%bleaf(ico),bl_max)
             elseif(phenology(ipft) == 2)then
           !     cpatch%phenology_status(ico) = 2
                bl_max = 0.0
                delta_bleaf = cpatch%bleaf(ico)
           !     cpatch%bleaf(ico) = 0.0
             endif
          endif

          if(csite%phenophase(ipft,ipa) == 4)then
             bl_max = csite%green_leaf_factor(ipft,ipa) *   &
                  dbh2bl(cpatch%dbh(ico),ipft)
             delta_bleaf = cpatch%bleaf(ico) - bl_max
!             if(delta_bleaf >= 0.0)then
!                cpatch%phenology_status(ico) = 0
!                cpatch%bleaf(ico) = bl_max
!             endif
          endif
         ! print*,'delta_bleaf',delta_bleaf
          
		 if(delta_bleaf < -1e-6) then!ATT
                      cpatch%phenology_status(ico) = 1
            else if (delta_bleaf > 1e-6) then
                      cpatch%phenology_status(ico) = -1 
            else
                      cpatch%phenology_status(ico) = 0
          end if!end ATT

          if(delta_bleaf > 0.0)then
               cpatch%leaf_drop(ico) = (1.0 - C_resorption_factor(ipft)) * delta_bleaf

               csite%fsc_in(ipa) = csite%fsc_in(ipa)                                       &
                                 + cpatch%nplant(ico) * delta_bleaf              &
                                 * f_labile(ipft) * f_fast(ipft) 						 &
								 * (1. - C_resorption_factor(ipft))
               csite%fsn_in(ipa) = csite%fsn_in(ipa)                                       &
                                 + cpatch%nplant(ico) * delta_bleaf              &
                                 * f_labile(ipft) * f_fast(ipft) / c2n_leaf(ipft)  &
								 * (1. - N_resorption_factor(ipft))
               csite%ssc_in(ipa) = csite%ssc_in(ipa)                                       &
                                 + cpatch%nplant(ico) * delta_bleaf              &
                                 * (1.0-f_labile(ipft)) 						 &
								 * (1. - C_resorption_factor(ipft))
               csite%ssl_in(ipa) = csite%ssl_in(ipa)                                       &
                                 + cpatch%nplant(ico) * delta_bleaf              &
                                 * (1.0 - f_labile(ipft)) * l2n_stem / c2n_stem(ipft) &
								 * (1. - C_resorption_factor(ipft))
			   csite%slsc_in(ipa) = csite%slsc_in(ipa)	&
			   					+ (1. - f_fast(ipft)) * delta_bleaf * f_labile(ipft) &
								* (1. - C_resorption_factor(ipft)) *	cpatch%nplant(ico)
			   csite%slsn_in(ipa) = csite%slsn_in(ipa)	&
			   					+ (1. - f_fast(ipft)) * delta_bleaf * f_labile(ipft) / c2n_leaf(ipft) &
								* (1. - N_resorption_factor(ipft)) *	cpatch%nplant(ico)
               !----- Adjust plant carbon pools. ------------------------------------------!
               cpatch%balive(ico)   = cpatch%balive(ico) - delta_bleaf
               cpatch%bstorage(ico) = cpatch%bstorage(ico) + C_resorption_factor(ipft)      &
                                    * delta_bleaf
			   cpatch%nstorage(ico) = cpatch%nstorage(ico) + N_resorption_factor(ipft)		&
			   						* delta_bleaf / c2n_leaf(ipft)
			   ! check whether nstorage reaches maximum
			   if (cpatch%nstorage(ico) > cpatch%nstorage_min(ico) *  nstorage_max_factor) then
			   	csite%fsn_in(ipa) = csite%fsn_in(ipa) + (cpatch%nstorage(ico) - &
						cpatch%nstorage_min(ico) * nstorage_max_factor) * f_fast(ipft) * cpatch%nplant(ico)
               
			   	csite%slsn_in(ipa) = csite%slsn_in(ipa) + (cpatch%nstorage(ico) - &
						cpatch%nstorage_min(ico) * nstorage_max_factor) * (1. - f_fast(ipft)) * cpatch%nplant(ico)
				cpatch%nstorage(ico) = cpatch%nstorage_min(ico) * nstorage_max_factor
			   endif
                cpatch%bleaf(ico)     = cpatch%bleaf(ico) - delta_bleaf               
               cpatch%cb(13,ico)     = cpatch%cb(13,ico)     - cpatch%leaf_drop(ico)
               cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) - cpatch%leaf_drop(ico)
          endif
          
          if(leaf_out_cold(ipft))then
             cpatch%phenology_status(ico) = 1
             cpatch%bleaf(ico) = (csite%green_leaf_factor(ipft,ipa) *  &
                  cpatch%balive(ico) /   &
                  (1.+qsw(ipft)*cpatch%hite(ico)+q(ipft)))
          endif
          call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico), &
               cpatch%bdead(ico), cpatch%balive(ico), cpatch%dbh(ico), &
               cpatch%hite(ico), cpatch%pft(ico), cpatch%sla(ico), &
               cpatch%lai(ico), cpatch%wpa(ico), cpatch%wai(ico), &
               cpatch%crown_area(ico), cpatch%bsapwood(ico))
          cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico), cpatch%balive(ico), &
               cpatch%bleaf(ico), cpatch%pft(ico), cpatch%hite(ico), &
               cpatch%bstorage(ico), cpatch%bsapwood(ico))

          old_leaf_hcap = cpatch%leaf_hcap(ico)
          old_wood_hcap = cpatch%wood_hcap(ico)
          call calc_veg_hcap(cpatch%bleaf(ico), &
               cpatch%bdead(ico), cpatch%bsapwood(ico), cpatch%nplant(ico), &
               cpatch%pft(ico), cpatch%leaf_hcap(ico), cpatch%wood_hcap(ico))
          call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap, old_wood_hcap)
          call is_resolvable(csite,ipa,ico,csite%green_leaf_factor(:,ipa))
       enddo
    enddo


    ! Declination angle at winter solstice
!    doy_use=355
!    day_angle = 2. * pi1 * (doy_use - 1) / 365.
!    solar_declination = 0.006918 - 0.399912 * cos(day_angle) + &
!         0.070257 * sin(day_angle) - 0.006758*cos(2.*day_angle) +  &
!         0.000907 * sin(2.*day_angle) - 0.002697 * cos(3.*day_angle) + &
!         0.001480 * sin(3.*day_angle)
!    acosarg = -tan(lat*pio180)*tan(solar_declination)
!    if(acosarg < -1.)acosarg = -1.
!    if(acosarg > 1.)acosarg = 1.
!    daylength_winter_solstice = 24./pi1 * acos(acosarg)

    do ipa = 1, csite%npatches
       do ipft = 1, n_pft
          if(phenology(ipft) == 2)then
             csite%leaf_aging_factor(ipft,ipa) =   &
                  min(1.,(daylength/daylength_solstice)**2)
          else
             csite%leaf_aging_factor(ipft,ipa) = 1.0
          endif
       enddo
    enddo


!    activation_energy = (activation_energy_max - ae_phen_amp) +   &
!         max(0.0,daylength-daylength_winter_solstice) * ae_phen_amp / &
!         (daylength_solstice-daylength_winter_solstice)
!print*,csite%green_leaf_factor(6,1),csite%green_leaf_factor(10,1),csite%leaf_aging_factor(6,1),csite%leaf_aging_factor(10,1)
!print*,csite%green_leaf_factor(6,1),csite%green_leaf_factor(10,1),pine_lai, oak_lai
!print*,doy,activation_energy
 
    return
  end subroutine update_jeong_phenology
  
end Module phenology_jeong
