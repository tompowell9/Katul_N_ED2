!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main driver for the longer-term vegetation dynamics.  This    !
! has become a file by itself to reduce the number of sub-routines that are doubled        !
! between ED-2.1 stand alone and the coupled model.                                        !
!------------------------------------------------------------------------------------------!
subroutine vegetation_dynamics(new_month,new_year)
   use grid_coms        , only : ngrids
   use ed_misc_coms     , only : current_time, iyeara           & ! intent(in)
                               , dtlsm                  & ! intent(in)
                               , frqsum                 & ! intent(in)
                               , ied_init_mode          ! ! intent(in)
   use disturb_coms     , only : include_fire           ! ! intent(in)
   use disturbance_utils, only : apply_disturbances     & ! subroutine
                               , site_disturbance_rates ! ! subroutine
   use fuse_fiss_utils  , only : fuse_patches           ! ! subroutine
   use ed_state_vars    , only : edgrid_g               & ! intent(inout)
                               , edtype,polygontype,sitetype,patchtype                 ! ! variable type
   use growth_balive    , only : dbalive_dt             & ! subroutine
                               , dbalive_dt_eq_0        ! ! subroutine
   use consts_coms      , only : day_sec                & ! intent(in)
                               , yr_day                 ! ! intent(in)
   use mem_polygons     , only : maxpatch               ! ! intent(in)
   use pft_coms, only: c2n_leaf, c2n_storage, c2n_recruit, &
                       c2n_stem, c2n_slow, c2n_structural, &
                       c2n_storage,n_pft,include_pft !AF 

  implicit none
  !----- Arguments. ----------------------------------------------------------------------!
  logical     , intent(in)   :: new_month
  logical     , intent(in)   :: new_year
  !----- Local variables. ----------------------------------------------------------------!
  type(edtype), pointer      :: cgrid
  real                       :: tfact1
  real                       :: tfact2
  integer                    :: doy
  integer                    :: ip
  integer                    :: isite
  integer                    :: ifm
  real                       :: oldmsn!JL!
  real                       :: newmsn !JL!
  integer                    :: i !ATT
  !----- External functions. -------------------------------------------------------------!
  integer     , external     :: julday
  !---------------------------------------------------------------------------------------!

  real(kind=8) :: oldpn, oldsn, newpn, newsn, oldtn, newtn,new_balive,new_bdead,new_bstorage,old_balive,old_bdead,old_bstorage


  type(polygontype), pointer :: cpoly
  type(sitetype), pointer :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipy, isi, ipa, ico

  oldpn=0.;oldtn=0.;oldsn=0.;newpn=0.;newsn=0.;newtn=0.;newmsn=0.;oldmsn=0. !old and newmsn were added JL
  old_balive = 0.;old_bdead = 0.; old_bstorage = 0.;  new_balive = 0.; new_bdead = 0.; new_bstorage = 0. 



   !----- Find the day of year. -----------------------------------------------------------!
   doy = julday(current_time%month, current_time%date, current_time%year)
  
   !----- Time factor for normalizing daily variables updated on the DTLSM step. ----------!
   tfact1 = dtlsm / day_sec
   !----- Time factor for averaging dailies. ----------------------------------------------!
   tfact2 = 1.0 / yr_day

   !----- Apply events. -------------------------------------------------------------------!
   call prescribed_event(current_time%year,doy)

  
   !---------------------------------------------------------------------------------------!
   !   Loop over all domains.                                                              !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids

      cgrid => edgrid_g(ifm) 

      do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            do ipa=1,csite%npatches
               oldsn = oldsn + (csite%slow_soil_N(ipa) + csite%fast_soil_N(ipa) +           &
                       csite%structural_soil_C(ipa)/c2n_structural +                        &
                       csite%mineralized_soil_N(ipa)) * csite%area(ipa)!ATT
               cpatch => csite%patch(ipa)
               !AF: Change the PFTS in the line of code below to optimize for your specific site and trial
               !add do loop below to remove pft dependence
               do i=1,n_pft
                  if(include_pft(i)) then
                     oldpn = oldpn + csite%area(ipa) * (csite%repro(i,ipa)/c2n_recruit(i))
                  end if
               end do
               oldmsn = oldmsn + csite%mineralized_soil_N(ipa) * csite%area(ipa)!JL!
               do ico=1,cpatch%ncohorts
                    oldpn = oldpn + csite%area(ipa) * cpatch%nplant(ico) *                    &
                            (cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico)) + cpatch%bdead(ico) &
                            /c2n_stem(cpatch%pft(ico))+cpatch%nstorage(ico)) !jl changed bstorage to nstorage
                    old_balive = old_balive + csite%area(ipa) * cpatch%nplant(ico) * cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))
                    old_bdead = old_bdead + csite%area(ipa) * cpatch%nplant(ico) * cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))
                    old_bstorage = old_bstorage + csite%area(ipa) * cpatch%nplant(ico) * cpatch%nstorage(ico)

              enddo
            enddo
         enddo
      enddo
!      print*,'DN, INITIAL', 'SoilN',oldsn,'PlantN',oldpn, 'ForestN',oldsn+oldpn,'MSN',oldmsn,'BaliveN', &
!              old_balive,'BdeadN',old_bdead,'NstorageN',old_bstorage

      !------------------------------------------------------------------------------------!
      !     The following block corresponds to the daily time-step.                        !
      !------------------------------------------------------------------------------------!
      !----- Standardise the fast-scale uptake and respiration, for growth rates. ---------!
      call normalize_ed_daily_vars(cgrid, tfact1)
	  call normalize_daily_vars(cgrid)  ! update leaf water potential
      !----- Update phenology and growth of live tissues. ---------------------------------!
      select case (ied_init_mode)
      case (-8)
         !----- Special case, in which we don't solve the actual vegetation dynamics. -----!
         call phenology_driver_eq_0(cgrid,doy,current_time%month, tfact1)
         call dbalive_dt_eq_0(cgrid,tfact2)
      case default
         call phenology_driver(cgrid,doy,current_time%month, tfact1)
		!	print*,'leave pheno, enter dbalive_dt'
         call dbalive_dt(cgrid,tfact2)
      end select
!	  print*,'LAI',cgrid%polygon(1)%site(1)%patch(1)%lai,'bstorage',cgrid%polygon(1)%site(1)%patch(1)%bstorage
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     The following block corresponds to the monthly time-step:                      !
      !------------------------------------------------------------------------------------!
      if (new_month) then

         !----- Update the mean workload counter. -----------------------------------------!
         call update_workload(cgrid)

         !----- Update the growth of the structural biomass. ------------------------------!
         call structural_growth(cgrid, current_time%month)

         !----- Update soil colums with appropriate organic layer thickness ATT
!         call update_organic_layer(cgrid) 

         !----- Solve the reproduction rates. ---------------------------------------------!
         call reproduction(cgrid,current_time%month)

         !----- Update the fire disturbance rates. ----------------------------------------!
         if (include_fire /= 0) then
            call fire_frequency(current_time%month,cgrid)
         end if

         !----- Update the disturbance rates. ---------------------------------------------!
         call site_disturbance_rates(current_time%month, current_time%year, cgrid)
      endif

      !------  update dmean and mmean values for NPP allocation terms ---------------------!
      call normalize_ed_dailyNPP_vars(cgrid)
      
      !------------------------------------------------------------------------------------!
      !     This should be done every day, but after the longer-scale steps.  We update    !
      ! the carbon and nitrogen pools, and re-set the daily variables.                     !
      !------------------------------------------------------------------------------------!
!			print*,'update C_and_N'
      call update_C_and_N_pools(cgrid)

      call zero_ed_daily_vars(cgrid)
!	  call reset_daily_vars(cgrid)  ! this should be called after output daily
!	  vars
      !------------------------------------------------------------------------------------!



     do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            do ipa=1,csite%npatches
               newsn = newsn + (csite%slow_soil_N(ipa) + csite%fast_soil_N(ipa) +           &
                       csite%structural_soil_C(ipa) / c2n_structural +                      &
                       csite%mineralized_soil_N(ipa)) * csite%area(ipa)!ATT

               cpatch => csite%patch(ipa)
               !AF: Change the PFTS in the line of code below to optimize for your specific site and trial
               do i=1,n_pft
                  if(include_pft(i)) then
                     newpn = newpn + csite%area(ipa) * (csite%repro(i,ipa)/c2n_recruit(i))!AF.
                  end if
               end do
               newmsn = newmsn + csite%mineralized_soil_N(ipa) * csite%area(ipa)!JL! 
               do ico=1,cpatch%ncohorts
                  newpn = newpn + csite%area(ipa) * cpatch%nplant(ico) * (cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico)) &
                            + cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))+cpatch%nstorage(ico))
        
                  new_balive = new_balive + csite%area(ipa) * cpatch%nplant(ico) * cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))
                  new_bdead = new_bdead + csite%area(ipa) * cpatch%nplant(ico) * cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))
                  new_bstorage = new_bstorage + csite%area(ipa) * cpatch%nplant(ico) * cpatch%nstorage(ico)
               enddo
            enddo
         enddo
      	  cgrid%ForestN(ipy)           = newsn+newpn !JL!
	      cgrid%PlantN(ipy)            = newpn !JL!
    	  cgrid%SoilN(ipy)             = newsn !JL!
      enddo
!      print*,'DN, FINAL','SoilN',newsn,'PlantN',newpn,'ForestN',newsn+newpn,'MSN', newmsn,'Balive',new_balive,'Bdead',  &
!              new_bdead,'Nstorage',new_bstorage ! 
     

         !----- This is actually the yearly time-step, apply the disturbances. ------------!
         if (new_month .and. new_year) then
            !comment out this update if you wish to use the regular ED2IN treefall disturbance
            !this subroutine ramps up (linearly) disturbance rate over a specified interval
            !call update_treefall_rate(iyeara, current_time%year)
     
            call apply_disturbances(cgrid)
  
         end if

      !------------------------------------------------------------------------------------!
      !      Fuse patches last, after all updates have been applied.  This reduces the     !
      ! number of patch variables that actually need to be fused.                          !
      !------------------------------------------------------------------------------------!
      if(new_year) then
         if (maxpatch >= 0) call fuse_patches(cgrid,ifm)
      end if
      !------------------------------------------------------------------------------------!



      !----- Recalculate the AGB and basal area at the polygon level. ---------------------!
      call update_polygon_derived_props(cgrid)
      call print_C_and_N_budgets(cgrid)
      !------------------------------------------------------------------------------------!


   end do

   return
end subroutine vegetation_dynamics
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
subroutine update_organic_layer(cgrid)
  use ed_state_vars, only: edtype, polygontype, sitetype, patchtype
  use grid_coms, only: nzg
  use soil_coms, only: slz, dslz,soil
  use consts_coms, only: alli, cice,cliq,t3ple,tsupercool
  use nutrient_constants, only: organic_soil_texture, organic_matter_density
  implicit none
  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype), pointer :: csite
  type(patchtype), pointer :: cpatch
  integer :: ipy,isi,ipa, k
  real :: new_litter_depth
  real :: new_layer_thickness, dryhcap, current_depth
  integer :: npix, npix_lev, offset, nshift, new_org_layer_index,  &
       old_org_layer_index
  real, allocatable, dimension(:) :: soil_tempk_old_pix,  &
       soil_tempk_new_pix, soil_psi_old_pix, soil_psi_new_pix,  &
       soil_fl_new_pix, soil_fl_old_pix
  real, dimension(nzg) :: new_psi, new_fl

  do ipy = 1, cgrid%npolygons
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1, cpoly%nsites
        csite => cpoly%site(isi)

        do ipa = 1, csite%npatches
           cpatch => csite%patch(ipa)

!           print*,'updating litter'

           ! Calculate current organic layer depth.  If there isn't one, then
           ! this will be set to nzg+1.
           old_org_layer_index = nzg + 1
           do k = nzg,1,-1
              ! need to make ntext_soil patch level if planning to add disturbance
              if(csite%ntext_soil(k,ipa) == organic_soil_texture)then
                 old_org_layer_index = k
              else
                 exit
              endif
           enddo

           ! Calculate the new litter layer depth.
           new_litter_depth = (csite%fast_soil_C(ipa) + csite%slow_soil_C(ipa))/organic_matter_density  ! Units: [m]  
!           print*,'new litter depth', new_litter_depth

           ! Calculate the index corresponding to the new organic layer depth
           new_org_layer_index = nzg + 1
           current_depth = 0.5 * dslz(nzg)
           if(new_litter_depth > current_depth)then
              new_org_layer_index = nzg
           endif
           do k = nzg-1,1,-1
              current_depth = current_depth + 0.5*dslz(k+1)+0.5*dslz(k)
              if(new_litter_depth > current_depth)then
                 new_org_layer_index = k
              endif
           enddo
!           print*,'new_org_layer_index', new_org_layer_index

           ! Don't let index change by more than one soil layer.  This 
           ! is mainly for computational stability.
           if(new_org_layer_index > old_org_layer_index)then
              new_org_layer_index = old_org_layer_index + 1
           endif
           if(new_org_layer_index < old_org_layer_index)then
              new_org_layer_index = old_org_layer_index - 1
           endif

           ! Adding a new litter layer
           if(new_org_layer_index < old_org_layer_index)then

              ! Additional litter thickness
              new_layer_thickness = dslz(old_org_layer_index-1)

              ! Subdivide the soil column into centimeter-sized pixels.
              npix = nint(-slz(1)*100.)

              ! Allocate temperature arrays corresponding to the pixels.
              allocate(soil_tempk_old_pix(npix))
              allocate(soil_tempk_new_pix(npix))
              allocate(soil_psi_old_pix(npix))
              allocate(soil_psi_new_pix(npix))
              allocate(soil_fl_old_pix(npix))
              allocate(soil_fl_new_pix(npix))

              ! Fill up the arrays
              offset = 0
              do k = 1, nzg
                 npix_lev = nint(dslz(k)*100.)
                 soil_tempk_old_pix((1+offset):(npix_lev+offset)) =  &
                      csite%soil_tempk(k,ipa)
                 soil_fl_old_pix((1+offset):(npix_lev+offset)) =  &
                      csite%soil_fracliq(k,ipa)
                 soil_psi_old_pix((1+offset):(npix_lev+offset)) =  &
                      soil(csite%ntext_soil(k,ipa))%slpots /   &
                      (csite%soil_water(k,ipa)/  &
                      soil(csite%ntext_soil(k,ipa))%slmsts)**  &
                      soil(csite%ntext_soil(k,ipa))%slbs
                 offset = offset + npix_lev
              enddo

              ! Size of the shift in pixel space
              nshift = nint(new_layer_thickness*100.)

              ! Shift the pixels
              do k = 1, npix-nshift
                 soil_tempk_new_pix(k) = soil_tempk_old_pix(k+nshift)
                 soil_psi_new_pix(k) = soil_psi_old_pix(k+nshift)
                 soil_fl_new_pix(k) = soil_fl_old_pix(k+nshift)
              enddo

              ! New pixels are the same as the pixels that were previously
              ! on top.
              do k = npix-nshift+1,npix
                 soil_tempk_new_pix(k) = soil_tempk_new_pix(npix-nshift)
                 soil_psi_new_pix(k) = soil_psi_new_pix(npix-nshift)
                 soil_fl_new_pix(k) = soil_fl_new_pix(npix-nshift)
              enddo

              ! Average over pixels to re-map to the SLZ space.
              offset = 0
              do k = 1, nzg
                 npix_lev = nint(dslz(k)*100.)
                 csite%soil_tempk(k,ipa) =   &
                      sum(soil_tempk_new_pix((1+offset):(npix_lev+offset))) / &
                      npix_lev
                 new_fl(k) =   &
                      sum(soil_fl_new_pix((1+offset):(npix_lev+offset))) / &
                      npix_lev
                 new_psi(k) =   &
                      sum(soil_psi_new_pix((1+offset):(npix_lev+offset))) / &
                      npix_lev
                 ! Set the new soil texture
                 csite%ntext_soil(old_org_layer_index-1,ipa) = organic_soil_texture
!                 print*,'ntextsoil', csite%ntext_soil(old_org_layer_index-1,ipa),'nzg',k
                 ! Compute the new soil water
                 csite%soil_water(k,ipa) =   &
                      (soil(csite%ntext_soil(k,ipa))%slpots / new_psi(k))**   &
                      (1./soil(csite%ntext_soil(k,ipa))%slbs) *  & 
                      soil(csite%ntext_soil(k,ipa))%slmsts
                 offset = offset + npix_lev
              enddo

              ! Get the new soil energy and fracliq.
              do k = 1, nzg
                 dryhcap = soil(csite%ntext_soil(k,ipa))%slcpd
                 if(csite%soil_tempk(k,ipa) < 273.15)then
                    csite%soil_fracliq(k,ipa) = 0.
                    csite%soil_energy(k,ipa) = csite%soil_tempk(k,ipa) * &
                         (cice * csite%soil_water(k,ipa) * 1000. + dryhcap)
                 elseif(csite%soil_tempk(k,ipa) > 273.15)then
                    csite%soil_fracliq(k,ipa) = 1.
                    csite%soil_energy(k,ipa) = csite%soil_tempk(k,ipa) * &
                         (cliq * csite%soil_water(k,ipa) * 1000. + dryhcap) - &
                         csite%soil_water(k,ipa)*1000. *cliq *tsupercool
                 else
                    csite%soil_fracliq(k,ipa) = new_fl(k)
                    csite%soil_energy(k,ipa) = new_fl(k) * &
                         csite%soil_water(k,ipa) * alli + (dryhcap + &
                         csite%soil_water(k,ipa)*1000.*cice)*t3ple
                 endif
              enddo

              deallocate(soil_tempk_old_pix)
              deallocate(soil_tempk_new_pix)
              deallocate(soil_psi_old_pix)
              deallocate(soil_psi_new_pix)
              deallocate(soil_fl_old_pix)
              deallocate(soil_fl_new_pix)
           endif

           !===========================================================

           ! Removing a litter layer
           if(new_org_layer_index > old_org_layer_index)then

              ! Thickness of the eliminated layer
              new_layer_thickness = dslz(old_org_layer_index)

              ! Subdivide the soil column into centimeter-sized pixels.
              npix = nint(-slz(1)*100.)

              ! Allocate temperature arrays corresponding to the pixels.
              allocate(soil_tempk_old_pix(npix))
              allocate(soil_tempk_new_pix(npix))
              allocate(soil_psi_old_pix(npix))
              allocate(soil_psi_new_pix(npix))
              allocate(soil_fl_old_pix(npix))
              allocate(soil_fl_new_pix(npix))

              offset = 0
              do k = 1, nzg
                 npix_lev = nint(dslz(k)*100.)
                 soil_tempk_old_pix((1+offset):(npix_lev+offset)) =  &
                      csite%soil_tempk(k,ipa)
                 soil_fl_old_pix((1+offset):(npix_lev+offset)) =  &
                      csite%soil_fracliq(k,ipa)
                 soil_psi_old_pix((1+offset):(npix_lev+offset)) =  &
                      soil(csite%ntext_soil(k,ipa))%slpots /   &
                      (csite%soil_water(k,ipa)/  &
                      soil(csite%ntext_soil(k,ipa))%slmsts)**  &
                      soil(csite%ntext_soil(k,ipa))%slbs
                 offset = offset + npix_lev
              enddo

              nshift = nint(new_layer_thickness*100.)

              do k = nshift+1,npix
                 soil_tempk_new_pix(k) = soil_tempk_old_pix(k-nshift)
                 soil_psi_new_pix(k) = soil_psi_old_pix(k-nshift)
                 soil_fl_new_pix(k) = soil_fl_old_pix(k-nshift)
              enddo

              do k = 1, nshift
                 soil_tempk_new_pix(k) = soil_tempk_old_pix(1)
                 soil_psi_new_pix(k) = soil_psi_old_pix(1)
                 soil_fl_new_pix(k) = soil_fl_old_pix(1)
              enddo

              offset = 0
              do k = 1, nzg
                 npix_lev = nint(dslz(k)*100.)
                 csite%soil_tempk(k,ipa) =   &
                      sum(soil_tempk_new_pix((1+offset):(npix_lev+offset))) / &
                      npix_lev
                 new_fl(k) =   &
                      sum(soil_fl_new_pix((1+offset):(npix_lev+offset))) / &
                      npix_lev
                 new_psi(k) =   &
                      sum(soil_psi_new_pix((1+offset):(npix_lev+offset))) / &
                      npix_lev
                 csite%ntext_soil(old_org_layer_index,ipa) =   &
                      csite%ntext_soil(old_org_layer_index-1,ipa)
                 csite%soil_water(k,ipa) =   &
                      (soil(csite%ntext_soil(k,ipa))%slpots / new_psi(k))**   &
                      (1./soil(csite%ntext_soil(k,ipa))%slbs) *  & 
                      soil(csite%ntext_soil(k,ipa))%slmsts

                 offset = offset + npix_lev
              enddo

              do k = 1, nzg
                 dryhcap = soil(csite%ntext_soil(k,ipa))%slcpd
                 if(csite%soil_tempk(k,ipa) < 273.15)then
                    csite%soil_fracliq(k,ipa) = 0.
                    csite%soil_energy(k,ipa) = csite%soil_tempk(k,ipa) * &
                         (cice * csite%soil_water(k,ipa) * 1000. + dryhcap)
                 elseif(csite%soil_tempk(k,ipa) > 273.15)then
                    csite%soil_fracliq(k,ipa) = 1.
                    csite%soil_energy(k,ipa) = csite%soil_tempk(k,ipa) * &
                         (cliq * csite%soil_water(k,ipa) * 1000. + dryhcap) - &
                         csite%soil_water(k,ipa)*1000. *cliq *tsupercool
                 else
                    csite%soil_fracliq(k,ipa) = new_fl(k)
                    csite%soil_energy(k,ipa) = new_fl(k) * &
                         csite%soil_water(k,ipa) * alli + (dryhcap + &
                         csite%soil_water(k,ipa)*1000.*cice)*t3ple
                 endif
              enddo

              deallocate(soil_tempk_old_pix)
              deallocate(soil_tempk_new_pix)
              deallocate(soil_psi_old_pix)
              deallocate(soil_psi_new_pix)
              deallocate(soil_fl_old_pix)
              deallocate(soil_fl_new_pix)
           endif

        enddo

     enddo

  enddo
!  print*,'ntextsoil', csite%ntext_soil

  return
end subroutine update_organic_layer

!==========================================================================================!
!==========================================================================================!

subroutine update_treefall_rate(first_year, current_year)
  use disturb_coms       , only: treefall_disturbance_rate
  use nutrient_constants , only: dist_start                    & 
                               , dist_end                      & 
                               , max_treefall_disturbance_rate                               
  implicit none

  integer, intent(in) :: first_year, current_year

  treefall_disturbance_rate = (min(dist_end,max(dist_start,current_year-first_year)) - dist_start) &
                               * (max_treefall_disturbance_rate/(dist_end-dist_start))
  return
end subroutine update_treefall_rate
!==========================================================================================!
!==========================================================================================!
