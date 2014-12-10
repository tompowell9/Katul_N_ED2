!==========================================================================================!
!==========================================================================================!
module growth_balive
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will update the alive biomass, compute the respiration terms  !
   ! other than leaf respiration and keep nitrogen budget.                                                          !
   ! IMPORTANT: The order of the operations here affect the C/N budgets, so don't change   !
   !            the order of the operations unless you really know what you are doing.     !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt(cgrid, tfact)
      use ed_state_vars    , only : edtype                 & ! structure
                                  , polygontype            & ! structure
                                  , sitetype               & ! structure
                                  , patchtype              ! ! structure
      use pft_coms         , only : q                      & ! intent(in)
                                  , qsw                    & ! intent(in)
                                  , c2n_storage            & ! intent(in)
                                  , c2n_leaf               & ! intent(in)              !JL!   
                                  , c2n_stem               & ! intent(in)              !JL!                
                                  , growth_resp_factor     & ! intent(in)
                                  , storage_turnover_rate  & ! intent(in)
                                  , mort3                  & ! intent(in)
                                  , root_beta              &
                                  , plant_N_supply_scale   &
                                  , phenology              ! ! intent(in)
      use physiology_coms  , only : N_plant_lim            ! ! intent(in)
      use grid_coms        , only : nzg                    ! ! intent(in)
      use ed_therm_lib     , only : calc_veg_hcap          & ! function
                                  , update_veg_energy_cweh ! ! function
      use allometry        , only : area_indices           & ! subroutine
                                  , ed_biomass             &! ! function
                                  , dbh2bl                 ! !JL
      use mortality        , only : mortality_rates        ! ! subroutine
      use phenology_coms   , only : theta_crit             ! ! intent(in)
      use soil_coms        , only : slz					  &
	  							  , nlsl
      use nutrient_constants, only: nstorage_max_factor
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: temp_dep
      real                          :: total_nl_broot
      real			    :: nl_broot
      real      		    :: nl_fraction
      real                          :: N_demand
      real                          :: N_supply
      real                          :: old_nstorage

      !------------------------------------------------------------------------------------!
!print*,'enter dbalive_dt'
      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

			   ! Calculate total broot in the nutrient layer
			   ! This will be used to calculate plant available nitrogen later
			   total_nl_broot = 0.
			   do ico = 1,cpatch%ncohorts
			  		nl_fraction = (1 - &  
							root_beta(cpatch%pft(ico)) ** &
							(-slz(nlsl) / &
							 (-slz(cpatch%krdepth(ico))) &
							 )) 
			   		total_nl_broot = total_nl_broot + &
							cpatch%broot(ico) *  nl_fraction * cpatch%nplant(ico)
					! the unit is kgC/m2

			   enddo

               !----- Loop over cohorts. --------------------------------------------------!
               do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Update the elongation factor. ------------------------------------!
!                  select case (phenology(ipft))
!                  case (4)
!                     cpatch%elongf(ico) = max(0.0, min(1.0, cpatch%paw_avg(ico)            &
!                                                          / theta_crit))
!                  case default
!                     cpatch%elongf(ico) = 1.0
!
!                  end select
!				These are not compatible with XXT's TDF phenology scheme


!if(ico == 2) print*,'bstorage before growth balive',cpatch%bstorage(ico)
                  
                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc
				  ! not sure whether we should keep it because now bleaf, bsap,
				  ! broot might be decoupled.....
                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs using actual pools.                      !
                  !------------------------------------------------------------------------!
!print*,'plant_main'
                  call plant_maintenance(cpatch,ico,cpatch%broot(ico),cpatch%bleaf(ico)    &
                                        ,tfact,daily_C_gain,csite%avg_daily_temp(ipa))

                 !----- Subtract maintenance costs from pools. ---------------------------!
                  cpatch%balive(ico)    = cpatch%balive(ico)                               &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%bleaf(ico)     = cpatch%bleaf(ico)                                &
                                        - cpatch%leaf_maintenance(ico)       
                  cpatch%broot(ico)     = cpatch%broot(ico)                                &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                                &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)                            &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)

                  !------------------------------------------------------------------------!
                  !     The commented line is an experimental and arbitrary test, borrowed !
                  ! from maintainence temperature dependency. [[MCD]]                      !
                  !------------------------------------------------------------------------!
                  ! temp_dep = 1.0                                                         &
                  !          / ( 1.0  + exp( 0.4 * (278.15 - csite%avg_daily_temp(ipa))))
                  temp_dep = 1.0
                  !------------------------------------------------------------------------!

                  cpatch%storage_respiration(ico) = cpatch%bstorage(ico)                   &
                                                  * storage_turnover_rate(ipft)            &
                                                  * tfact * temp_dep

                  cpatch%bstorage(ico) = cpatch%bstorage(ico)                              &
                                         - cpatch%storage_respiration(ico)

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances. Track    !
                  !      carbon available for growth to use for nitrogen calculations      !
                  !------------------------------------------------------------------------!
!print*,'plant_carbon_balance'
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_max)
!                     cpatch%carbon_balance(ico) = carbon_balance !JL!

                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------! 
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  cpatch%growth_respiration(ico) = max(0.0, daily_C_gain                   &
                                                          * growth_resp_factor(ipft))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the "virtual" leaf respiration for the coming day.            !
                  !------------------------------------------------------------------------!
				  !not sure what it does but generally storage_turnover_rate is 0
				  !
                  cpatch%vleaf_respiration(ico) = (1.0-csite%green_leaf_factor(ipft,ipa))  &
                                                * salloci * cpatch%balive(ico)             &
                                                * storage_turnover_rate(ipft)              &
                                                * tfact * temp_dep
                  !------------------------------------------------------------------------!

	        ! put resorption and unresorbed C& N in litter()



                  if (N_plant_lim == 1) then

                     !-----------------------------------------------------------------------!
                     ! XXT & ATT
                     ! Calculate fsn according to carbon_balance_pot
                     ! ----------------------------------------------------------------------!
                     nl_fraction = (1 - &  
                          root_beta(ipft) ** &
                          (-slz(nlsl) / &
                          (-slz(cpatch%krdepth(ico))) &
                          )) 
                     nl_broot = cpatch%broot(ico) * nl_fraction
                     
!                     N_supply = cpatch%nstorage(ico) + &
!                          nl_broot / total_nl_broot  * &  ! kgC plant-1 / kgC 	m-2
!                          csite%mineralized_soil_N(ipa)
                     ! Assumption about N demand
                     ! Here we assume the potential demand of trees is the nitrogen
                     ! required to put all the carbon balance into leaf/root which
                     ! is nitrogen rich.

					 ! XXT calculated N_uptake in rk4_driver coupled with water
					 cpatch%nitrogen_uptake(ico) = min(nl_broot / total_nl_broot *		 &
												 csite%mineralized_soil_N(ipa),&
							 					cpatch%nstorage_min(ico) * &
												nstorage_max_factor -  &
												cpatch%nstorage(ico))
					 csite%total_plant_nitrogen_uptake(ipa) =  &
					 csite%total_plant_nitrogen_uptake(ipa) + &
				 	cpatch%nitrogen_uptake(ico) * cpatch%nplant(ico)

					cpatch%nstorage(ico) = cpatch%nstorage(ico) + cpatch%nitrogen_uptake(ico)
					 N_supply = cpatch%nstorage(ico)
                     
                     N_demand = max(0.0,carbon_balance_pot / c2n_leaf(ipft)) ! in case cb is negative

                     cpatch%fsn(ico) = N_supply / (N_demand + N_supply)    ! used  for the next day
					 
					 ! consider extreme case
					 if (N_demand + N_supply == 0.) then
					 	cpatch%fsn(ico) = 0.0
					 endif

                     !-----------------------------------------------------------------------!
                     ! XXT & ATT														
                     ! Update nstorage according to soil nitrogen
                     ! can't exceed max_nstorage
                     ! ----------------------------------------------------------------------!
!                     old_nstorage = cpatch%nstorage(ico) 
!                     cpatch%nstorage(ico) = min(nstorage_max_factor * cpatch%nstorage_min(ico),N_supply)
!	if(ico == 1) print*,'N_supply', N_supply,'cpatch%nstorage',cpatch%nstorage(ico),'fsn',cpatch%fsn(ico)
!                     csite%total_plant_nitrogen_uptake(ipa) =  & 
!                          csite%total_plant_nitrogen_uptake(ipa) + &
!                          (cpatch%nstorage(ico) - old_nstorage) * cpatch%nplant(ico)
	if(isnan(csite%total_plant_nitrogen_uptake(ipa))) then
		print*,'nan plant update ico',ico,'old_nstorage',old_nstorage,'N_supply',N_supply,&
				'nstorage',cpatch%nstorage(ico),'nlsl',nlsl
	endif
                  else  ! No nitrogen limitation
                     cpatch%fsn(ico) = 1.0
                     cpatch%nstorage(ico) = 10. * cpatch%nstorage_min(ico)
                     ! give the plant plenty of nitrogen so that they are not
                     ! limited
                  endif

				  cpatch%nitrogen_uptake(ico) = 0.0  ! zero N_uptake

                  !------------------------------------------------------------------------!
                  !      Allocate plant carbon balance to balive and bstorage.             !
                  !------------------------------------------------------------------------!
                  balive_in = cpatch%balive(ico)
                
!print*,'alloc_plant_c'
                  call alloc_plant_c_balance(csite,ipa,ico,salloc,salloci,carbon_balance   &
                                            ,csite%green_leaf_factor(ipft,ipa))
                  !------------------------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
!print*,'mortality'
                  call mortality_rates(cpatch,ipa,ico,csite%avg_daily_temp(ipa)            &
                                      ,csite%age(ipa))
                  dndt = - sum(cpatch%mort_rate(:,ico)) * cpatch%nplant(ico) * tfact

                  !------- Update monthly mortality rate [plants/m2/month]. ---------------!
                  cpatch%monthly_dndt(ico) = cpatch%monthly_dndt(ico) + dndt

              
                   !----- Updating LAI, WPA, and WAI. --------------------------------------!
                  call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)                   &
                                   ,cpatch%bdead(ico),cpatch%balive(ico),cpatch%dbh(ico)   &
                                   ,cpatch%hite(ico) ,cpatch%pft(ico),cpatch%sla(ico)      &
                                   ,cpatch%lai(ico),cpatch%wpa(ico),cpatch%wai(ico)        &
                                   ,cpatch%crown_area(ico),cpatch%bsapwood(ico))

                  !----- Update above-ground biomass. -------------------------------------!
                  cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)        &
                                              ,cpatch%bleaf(ico),cpatch%pft(ico)           &
                                              ,cpatch%hite(ico),cpatch%bstorage(ico)       &
                                              ,cpatch%bsapwood(ico))

                  !------------------------------------------------------------------------!
                  !     It is likely that biomass has changed, therefore, update           !
                  ! vegetation energy and heat capacity.                                   !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap         = cpatch%leaf_hcap(ico)
                  old_wood_hcap         = cpatch%wood_hcap(ico)

                  call calc_veg_hcap(cpatch%bleaf(ico) ,cpatch%bdead(ico)                  &
                                    ,cpatch%bsapwood(ico),cpatch%nplant(ico)               &
                                    ,cpatch%pft(ico)                                       &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))

 
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)

                 !----- Update the stability status. -------------------------------------!
                  call is_resolvable(csite,ipa,ico,csite%green_leaf_factor(:,ipa))

!if(ico == 2) print*,'bstorage after growth balive',cpatch%bstorage(ico)
               end do
                   !----- Update litter. ----------------------------------------------------!
!print*,'litter'
                    call litter(csite,ipa)
                   !------------------------------------------------------------------------!
                     
               !----- Update patch LAI, WAI, height, roughness... -------------------------!
               call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,ipa)

               !----- Recalculate storage terms (for budget assessment). ------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
            end do
         end do
      end do
!print*,'leave growth_balive'
      return
   end subroutine dbalive_dt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will compute the respiration terms other than leaf                !
   ! respiration, plus the carbon balance and maintenance costs but without                !
   ! updating the pools.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt_eq_0(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           & ! subroutine
                                 , ed_biomass             ! ! function
      use mortality       , only : mortality_rates        ! ! subroutine
      use phenology_coms  , only : theta_crit             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake
      real                          :: N_uptake_pot
      !------------------------------------------------------------------------------------!

!print*,'enter growth_balive'
      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Loop over cohorts. --------------------------------------------------!
               do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Update the elongation factor. ------------------------------------!
                  select case (phenology(ipft))
                  case (4)
                     cpatch%elongf(ico) = max(0.0, min(1.0, cpatch%paw_avg(ico)/theta_crit))
                  case default
                     cpatch%elongf(ico) = 1.0
                  end select

                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0

                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc
                  
                  !----- Leaf and root biomass. -------------------------------------------!
                  bl = cpatch%bleaf(ico)
                  br = cpatch%broot(ico)

                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs.                                         !
                  !------------------------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain               &
                                        ,csite%avg_daily_temp(ipa))

                  !----- Subtract maintenance costs from balive. --------------------------!
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                                &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)                            &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_max)

                  !------------------------------------------------------------------------!
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  cpatch%growth_respiration(ico)  = max(0.0, daily_C_gain                  &
                                                           * growth_resp_factor(ipft))
                  cpatch%storage_respiration(ico) = cpatch%bstorage(ico)                   &
                                                  * storage_turnover_rate(ipft) * tfact
                  cpatch%vleaf_respiration(ico) =                                          &
                                        (1.0 - csite%green_leaf_factor(ipft,ipa))          &
                                      / (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))     &
                                      * cpatch%balive(ico) * storage_turnover_rate(ipft)   &
                                      * tfact

                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
!                  if (N_plant_lim == 1) then
 !                    call potential_N_uptake(cpatch,ico,salloc,salloci,balive_in           &
  !                                          ,carbon_balance_pot,N_uptake_pot               &
   !                                         ,csite%green_leaf_factor(ipft,ipa))
    !              end if

                  !------------------------------------------------------------------------!
                  !  Increment the [kgN/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                 ! csite%total_plant_nitrogen_uptake(ipa) =                                 &
                 !                                 csite%total_plant_nitrogen_uptake(ipa)  &
                 !                              + nitrogen_uptake * cpatch%nplant(ico)

                  !----- Calculate plant N limitation factor. -----------------------------!
                  if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0) then
                     cpatch%fsn(ico) = 1.0
                  else
                     nitrogen_supply = plant_N_supply_scale * br                           &
                                     * csite%mineralized_soil_N(ipa)
                     cpatch%fsn(ico) = nitrogen_supply / (nitrogen_supply + N_uptake_pot)
                  end if
                  
                  !------------------------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ipa,ico,csite%avg_daily_temp(ipa)            &
                                      ,csite%age(ipa))
               end do

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
            end do
         end do
      end do

      return
   end subroutine dbalive_dt_eq_0
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will transfer some of the stored carbon to balive in order to put  !
   ! the plant back on allometry.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine transfer_C_from_storage(cpatch,ico,salloc,salloci,nitrogen_uptake            &
                                     ,N_uptake_pot)
      use ed_state_vars , only : patchtype
      use pft_coms      , only : c2n_leaf    & ! intent(in)
                               , c2n_storage & ! intent(in)
                               , c2n_stem    & ! intent(in)
                               , q           & ! intent(in)
                               , qsw         ! ! intent(in)
      use decomp_coms   , only : f_labile    ! ! intent(in)
      use allometry     , only : dbh2bl      ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(inout) :: N_uptake_pot
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: off_allometry_cb
      real                           :: increment
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Only do the transfer if leaves exist.                                          !
      !------------------------------------------------------------------------------------!
      if (cpatch%phenology_status(ico) == 2) return
     
      !----- Alias for pft type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)
     
      !----- Determine how much biomass we need to go back to allometry. ------------------!
      off_allometry_cb = dbh2bl(cpatch%dbh(ico),ipft) * salloc - cpatch%balive(ico)

      !----- If plants have storage, transfer it to balive. -------------------------------!
      increment            = max(0.0,min(max(0.0, off_allometry_cb),cpatch%bstorage(ico)))
      cpatch%balive(ico)   = cpatch%balive(ico)   + increment
      cpatch%bstorage(ico) = cpatch%bstorage(ico) - increment

      !----- Compute sapwood and fine root biomass. ---------------------------------------!
      cpatch%broot(ico)    = q(ipft) * cpatch%balive(ico) * salloci
      cpatch%bsapwood(ico) = qsw(ipft) * cpatch%hite(ico) * cpatch%balive(ico) * salloci

      !------------------------------------------------------------------------------------!
      !      N uptake is required since c2n_leaf < c2n_storage.  Units are kgN/plant/day.  !
      !------------------------------------------------------------------------------------!
      nitrogen_uptake = increment * (        f_labile(ipft)  / c2n_leaf(ipft)              &
                                    + (1.0 - f_labile(ipft)) / c2n_stem(ipft)              &
                                    -  1.0 / c2n_storage)
      N_uptake_pot    = nitrogen_uptake

      return
   end subroutine transfer_C_from_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain,tempk)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : phenology          & ! intent(in)
                              , root_turnover_rate & ! intent(in)
                              , leaf_turnover_rate ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: br
      real           , intent(in)    :: bl
      real           , intent(in)    :: tfact
      real           , intent(in)    :: tempk
      real           , intent(out)   :: daily_C_gain
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: maintenance_temp_dep
      !------------------------------------------------------------------------------------!

      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Get the temperature dependence. ---------------------------------------------!
      if (phenology(ipft) == 0) then   ! Evergreen
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk))) 
!		 don't know where it comes from...
      else
         maintenance_temp_dep = 1.0
      end if

      !----- Calculate maintenance demand (kgC/plant/year). -------------------------------!
      cpatch%root_maintenance(ico) = root_turnover_rate(ipft) * br * maintenance_temp_dep

	  if (phenology(ipft) /= 3)then  
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl * maintenance_temp_dep
      else   ! Kim's radiation phenology
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl                      &
                                      * cpatch%turnover_amp(ico) * maintenance_temp_dep
      end if


      !----- Convert units of maintenance to [kgC/plant/day]. -----------------------------!
      cpatch%leaf_maintenance(ico) = cpatch%leaf_maintenance(ico) * tfact
      cpatch%root_maintenance(ico) = cpatch%root_maintenance(ico) * tfact


      !----- Compute daily C uptake [kgC/plant/day]. --------------------------------------!
      if(cpatch%nplant(ico) > tiny(1.0)) then
         daily_C_gain = umol_2_kgC * day_sec * ( cpatch%today_gpp(ico)                     &
                                               - cpatch%today_leaf_resp(ico)               &
                                               - cpatch%today_root_resp(ico))              &
                                             / cpatch%nplant(ico)
!	if (ico == 10) print*,'npp/gpp', daily_C_gain / umol_2_kgC / day_sec  * &
!	cpatch%nplant(ico) / cpatch%today_gpp(ico)
      else
         daily_C_gain = 0.0
      end if

      return


   end subroutine plant_maintenance
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance             &
                                   ,carbon_balance_pot,carbon_balance_max)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : growth_resp_factor ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      use ed_max_dims  , only : n_pft              ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ipa
      integer                  , intent(in)  :: ico
      real                     , intent(in)  :: daily_C_gain
      real                     , intent(out) :: carbon_balance
      real                     , intent(out) :: carbon_balance_pot
      real                     , intent(out) :: carbon_balance_max
      !----- Local variables. -------------------------------------------------------------!
      real                                   :: daily_C_gain_pot
      real                                   :: daily_C_gain_max
      real                                   :: growth_respiration_pot
      real                                   :: growth_respiration_max
      integer                                :: ipft
      !----- Local constants. -------------------------------------------------------------!
      logical                  , parameter   :: print_debug = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical, dimension(n_pft), save        :: first_time  = .true.
      !------------------------------------------------------------------------------------!

      !----- Alias for PFT type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Calculate actual daily carbon balance: kgC/plant/day. -----------------------!
      carbon_balance = daily_C_gain - cpatch%growth_respiration(ico)                       &
                                    - cpatch%vleaf_respiration(ico)

      if (cpatch%nplant(ico) > tiny(1.0)) then

         !---------------------------------------------------------------------------------!
         !      Calculate potential carbon balance (used for nitrogen demand function).    !
         ! [kgC/plant/day].                                                                !
         !---------------------------------------------------------------------------------!
         daily_C_gain_pot       = umol_2_kgC * day_sec * ( cpatch%today_gpp_pot(ico)       &
                                                         - cpatch%today_leaf_resp(ico)     &
                                                         - cpatch%today_root_resp(ico))    &
                                                       / cpatch%nplant(ico)
         growth_respiration_pot = max(0.0, daily_C_gain_pot * growth_resp_factor(ipft))

         carbon_balance_pot = daily_C_gain_pot - growth_respiration_pot			           &
                                               - cpatch%vleaf_respiration(ico)
         
       !----- Calculate maximum carbon balance (used for mortality). --------------------!
         daily_C_gain_max       = umol_2_kgC * day_sec * ( cpatch%today_gpp_max(ico)       &
                                                         - cpatch%today_leaf_resp(ico)     &
                                                         - cpatch%today_root_resp(ico) )   &
                                                       / cpatch%nplant(ico)
         growth_respiration_max = max(0.0, daily_C_gain_max * growth_resp_factor(ipft))
         carbon_balance_max     = daily_C_gain_max - growth_respiration_max                &
                                                   - cpatch%vleaf_respiration(ico)
      else
         carbon_balance_max = 0.0
         carbon_balance_pot = 0.0
      end if

      !----- Carbon balances for mortality. -----------------------------------------------!
      cpatch%cb(13,ico)     = cpatch%cb(13,ico)     + carbon_balance
      cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) + carbon_balance_max

      if (print_debug) then

         if (first_time(ipft)) then
            first_time(ipft) = .false.
            write (unit=30+ipft,fmt='(a10,15(1x,a12))')                                    &
                '      TIME','       PATCH','      COHORT','      NPLANT','    CB_TODAY'   &
                            ,' GROWTH_RESP','  VLEAF_RESP','   TODAY_GPP','TODAY_GPPMAX'   &
                            ,'  TODAY_LEAF','  TODAY_ROOT',' CBMAX_TODAY','          CB'   &
                            ,'       CBMAX','  LEAF_MAINT','  ROOT_MAINT'
         end if

         write(unit=30+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),13(1x,es12.5))')               &
              current_time%month,'/',current_time%date,'/',current_time%year               &
             ,ipa,ico,cpatch%nplant(ico),carbon_balance,cpatch%growth_respiration(ico)     &
             ,cpatch%vleaf_respiration(ico),cpatch%today_gpp(ico)                          &
             ,cpatch%today_gpp_max(ico),cpatch%today_leaf_resp(ico)                        &
             ,cpatch%today_root_resp(ico),carbon_balance_max,cpatch%cb(13,ico)             &
             ,cpatch%cb_max(13,ico),cpatch%leaf_maintenance(ico)                           &
             ,cpatch%root_maintenance(ico)
      end if

      return

   end subroutine plant_carbon_balances
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_plant_c_balance(csite,ipa,ico,salloc,salloci,carbon_balance            &
                                   ,green_leaf_factor)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , sla          & ! intent(in)
                               , q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     &! ! intent(in)
                               , f_fast  
      use allometry     , only : dbh2bl       &! ! function
	                       , dbh2bs 
      use nutrient_constants, only : nstorage_max_factor
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: carbon_balance
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bl_max,br_max,bsapwood_max
      real                           :: balive_max
      real                           :: increment
      real                           :: total_burnoff
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwood
      real                           :: available_carbon
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwood
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwood
      real		             :: sapwood_burnoff, extra_nitrogen
      real                           :: nitrogen_demand, carbon_demand
      logical                        :: on_allometry

      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 
	  
	  ! There's no need to have different branches with different leaf fullness
	  ! as before

      if(carbon_balance + cpatch%bstorage(ico) > 0.)then
      
         available_carbon = cpatch%bstorage(ico) + carbon_balance
!         if (cpatch%phenology_status(ico) /= -1) then 
            ! when trees are growing
            !------------------------------------------------------------------------------!
            !     Maximum bleaf that the allometric relationship would allow.  If the      !
            ! plant is drought stress (elongf < 1), we do not allow the plant to get back  !
            ! to full allometry.                                                           !
            !------------------------------------------------------------------------------!
            bl_max     = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor                  &
                 * cpatch%elongf(ico)
            br_max     = dbh2bl(cpatch%dbh(ico),ipft) * q(ipft) * (	cpatch%elongf(ico) + 1.0) / 2.0
            bsapwood_max = dbh2bs(cpatch%dbh(ico),ipft)
			
            delta_bleaf = max(0.,bl_max - cpatch%bleaf(ico))
            delta_broot = max(0.,br_max - cpatch%broot(ico))
            delta_bsapwood = max(0.,bsapwood_max - cpatch%bsapwood(ico))

            carbon_demand = delta_bleaf + delta_broot + delta_bsapwood
	!		carbon_demand can't be negative
            
            ! first determine the actual carbon demand and nitrogen demand
            if (available_carbon >= carbon_demand) then ! enough carbon to grow to on-allometry 
               nitrogen_demand = carbon_demand / c2n_leaf(ipft)
            else ! not enough carbon grow to on-allometry
               carbon_demand = available_carbon
               nitrogen_demand = carbon_demand / c2n_leaf(ipft)
            endif

            ! then update the bstorage and nstorage
            if (cpatch%nstorage(ico) >= nitrogen_demand) then
               cpatch%bstorage(ico) = available_carbon - carbon_demand 
               cpatch%nstorage(ico) = cpatch%nstorage(ico) - nitrogen_demand
            else ! not enough nitrogen
               carbon_demand = carbon_demand * (cpatch%nstorage(ico) /	nitrogen_demand)
               cpatch%nstorage(ico) = 0.0
               cpatch%bstorage(ico) = available_carbon - carbon_demand
            endif

            ! let broot and bleaf grow first
			if (delta_bleaf == 0. .and. delta_broot == 0.) then
				f_bleaf = 0.
				f_broot = 0.
			else
            	f_bleaf = delta_bleaf / (delta_bleaf +delta_broot)
	            f_broot = delta_broot / (delta_bleaf +delta_broot)
			endif
            
            tr_bleaf = min(delta_bleaf,carbon_demand * f_bleaf)
            tr_broot = min(delta_broot,carbon_demand * f_broot)
            tr_bsapwood = min(delta_bsapwood,carbon_demand - tr_bleaf - tr_broot)
!if(ico == 1) print*,'delta_bleaf',delta_bleaf,'delta_broot',delta_broot
!if(ico == 1) print*,'delta_sapwood',delta_bsapwood,'carbon_demand',carbon_demand
!if(ico == 1) print*,'tr_bleaf',tr_bleaf,'tr_broot',tr_broot,'tr_bsapwood',tr_bsapwood
            cpatch%bleaf(ico) = cpatch%bleaf(ico) +  tr_bleaf
            cpatch%broot(ico) = cpatch%broot(ico) +  tr_broot
            cpatch%bsapwood(ico) = cpatch%bsapwood(ico) + tr_bsapwood
            !finally update balive
            cpatch%balive(ico) = cpatch%balive(ico) + carbon_demand

            cpatch%today_nppleaf(ico) = tr_bleaf * cpatch%nplant(ico)
            cpatch%today_nppfroot(ico) = tr_broot * cpatch%nplant(ico)
            cpatch%today_nppsapwood(ico) = tr_bsapwood * cpatch%nplant(ico)
            cpatch%today_nppdaily(ico) = carbon_balance * cpatch%nplant(ico)

            ! update phenology_status
			! the phenology_status should be updated in phenology_drive
!            if (abs(cpatch%bleaf(ico) - bl_max) < 1.0e-8) then
!               if (cpatch%elongf(ico) * csite%green_leaf_factor(ipft,ipa) > 0.9) then
!                  cpatch%phenology_status(ico) = 0
!               else
!                  cpatch%phenology_status(ico) = 1
!               endif
!            elseif (cpatch%bleaf(ico) > 0.) then
!               cpatch%phenology_status(ico) = 1
!            else
!               cpatch%phenology_status(ico) = 2
!            endif
!         else  ! The canopy is full
!            cpatch%bstorage(ico) = available_carbon
!            cpatch%today_nppleaf(ico) = 0.
!            cpatch%today_nppfroot(ico) = 0.
!            cpatch%today_nppsapwood(ico) = 0.
!            cpatch%today_nppdaily(ico) = carbon_balance * cpatch%nplant(ico)
!         endif


      else   ! CB + bstorage < 0
         increment =  cpatch%bstorage(ico) + carbon_balance
         !----- Use Storage pool first then take out of balive. ------------------------!
         total_burnoff            =  - increment
         cpatch%bstorage(ico) = 0.0
         
         ! distribute increment according to bleaf, broot
         ! if burn off bleaf and broot is not enough, burn off sapwood
         sapwood_burnoff = max(0.0,total_burnoff - (cpatch%bleaf(ico) + &
              cpatch%broot(ico)))  
		 if (cpatch%bleaf(ico) == 0. .and. cpatch%broot(ico) == 0.) then
		 	f_bleaf = 0.
			f_broot = 0.
		 else
	         f_bleaf = cpatch%bleaf(ico) / (cpatch%bleaf(ico) + cpatch%broot(ico))
    	     f_broot = cpatch%broot(ico) / (cpatch%bleaf(ico) + cpatch%broot(ico))
		 endif

		 
         cpatch%bleaf(ico)    = max(0.0,cpatch%bleaf(ico) - f_bleaf *   (total_burnoff - sapwood_burnoff))
         cpatch%broot(ico)    = max(0.0,cpatch%broot(ico) - f_broot *   (total_burnoff - sapwood_burnoff))

         cpatch%bsapwood(ico) = max(0.0,cpatch%bsapwood(ico) -   sapwood_burnoff)
 
         cpatch%balive(ico)   = max(0.,cpatch%balive(ico) - total_burnoff)

!         if (cpatch%phenology_status(ico) >= 0) then
!            !We only need to update phenology_status when its bigger than 0	
!            if (cpatch%bleaf(ico) > 0 )then
!               cpatch%phenology_status(ico) = 1
!            else
!               print*,'ico',ico,'All leaves burn'
!              cpatch%phenology_status(ico) = 2  ! no leaf
!            endif
!            
!         endif
			   
         ! The burn-off carbon is treated as CO2
         ! The nitrogen associated will go into nstorage or into soil if
         ! nstorage is full
         extra_nitrogen = total_burnoff / c2n_leaf(ipft) 
         if (cpatch%nstorage(ico) + extra_nitrogen <= cpatch%nstorage_min(ico) * nstorage_max_factor) then
            cpatch%nstorage(ico) = cpatch%nstorage(ico) + extra_nitrogen
            extra_nitrogen = 0.
         else
            extra_nitrogen = extra_nitrogen - (cpatch%nstorage_min(ico)					&
                 * nstorage_max_factor - cpatch%nstorage(ico))
            cpatch%nstorage(ico) = cpatch%nstorage_min(ico) * nstorage_max_factor
         endif

         ! the extra_nitrogen goes to soil pools
         csite%slsn_in(ipa) = csite%slsn_in(ipa) +  &
              extra_nitrogen * (1.-f_fast(ipft)) &
              * cpatch%nplant(ico)

         csite%fsn_in(ipa) = csite%fsn_in(ipa) +  &
              extra_nitrogen * (f_fast(ipft)) &
              * cpatch%nplant(ico)
		
         cpatch%today_nppleaf(ico) = 0.
         cpatch%today_nppfroot(ico) = 0.
         cpatch%today_nppsapwood(ico) = 0.
         cpatch%today_nppdaily(ico) = carbon_balance * cpatch%nplant(ico)
      endif
         
      return
   end subroutine alloc_plant_c_balance
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   subroutine litter(csite,ipa)

      use ed_state_vars      , only : patchtype & ! structure
                             , sitetype  ! ! structure
      use pft_coms           , only : c2n_leaf  & ! intent(in)
                             , c2n_stem  & ! intent(in)
                             , l2n_stem  & ! intent(in)
                             , N_resorption_factor &
                             , C_resorption_factor
      use decomp_coms        , only : f_labile & ! ! intent(in)
                             , f_fast
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: plant_litter
      real                         :: plant_litter_f
      real                         :: plant_litter_s
       !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !      Add fine root and leaf turnover to the litter.                                !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)
         !------------------------------------------------------------------------!
         ! Calculate C & N resorption in plant maintenance
    	 ! Currently both C and N storage don't have a maximum
		 ! constraint....          
		 !
         !------------------------------------------------------------------------! 
		 cpatch%nstorage(ico) = cpatch%nstorage(ico) +			       	  &
		 		(cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) &
				 / c2n_leaf(ipft) * N_resorption_factor(ipft)
		 cpatch%bstorage(ico) = cpatch%bstorage(ico) +				  &
				 (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) &
                                 * C_resorption_factor(ipft)

         !------------------------------------------------------------------------! 
         ! Unresorbed maintenance goes to fast/slow soil pools according to f_fast!  
         !------------------------------------------------------------------------!
 	      csite%fsn_in(ipa) = csite%fsn_in(ipa) + f_fast(ipft) * &
			  (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) / &
			   c2n_leaf(ipft) * (1. - N_resorption_factor(ipft)) *   cpatch%nplant(ico)
		
		  csite%slsn_in(ipa) = csite%slsn_in(ipa) + (1. - f_fast(ipft)) * &
			  (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) / &
			  c2n_leaf(ipft) * (1. - N_resorption_factor(ipft)) *   cpatch%nplant(ico)

		  csite%fsc_in(ipa) = csite%fsc_in(ipa) + f_fast(ipft) * &
		  	  (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) * &
			  f_labile(ipft)  * &
			  (1. - C_resorption_factor(ipft)) * cpatch%nplant(ico)

		  csite%slsc_in(ipa) = csite%slsc_in(ipa) + (1. - f_fast(ipft)) * &
		  	  (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) * &
			  f_labile(ipft) * &
			  (1. - C_resorption_factor(ipft)) * cpatch%nplant(ico)

		  csite%ssc_in(ipa) = csite%ssc_in(ipa)  + (1. - f_labile(ipft)) * &
		  	  (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) * &
			  (1. - C_resorption_factor(ipft)) * cpatch%nplant(ico)

		  csite%ssl_in(ipa) = csite%ssl_in(ipa)  + (1. - f_labile(ipft)) * &
		  	  (cpatch%leaf_maintenance(ico)+cpatch%root_maintenance(ico)) * &
			  (1. - C_resorption_factor(ipft)) * cpatch%nplant(ico) * &
			  l2n_stem / c2n_stem(ipft)
         
      end do
      return
   end subroutine litter
   !=======================================================================================!
   !=======================================================================================!
  
 
end module growth_balive
!==========================================================================================!
!==========================================================================================!
