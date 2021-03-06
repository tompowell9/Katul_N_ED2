!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the structural growth of plants.                        !
!                                                                                          !
! IMPORTANT: Do not change the order of operations below unless you know what you are      !
!            doing.  Changing the order can affect the C/N budgets.                        !
!------------------------------------------------------------------------------------------!
subroutine structural_growth(cgrid, month)
   use ed_state_vars , only : edtype                 & ! structure
                            , polygontype            & ! structure
                            , sitetype               & ! structure
                            , patchtype              ! ! structure
   use pft_coms      , only : q                      & ! intent(in)
                            , qsw                    & ! intent(in)
                            , seedling_mortality     & ! intent(in)
                            , c2n_leaf               & ! intent(in)
                            , c2n_storage            & ! intent(in)
                            , c2n_recruit            & ! intent(in)
                            , c2n_stem               & ! intent(in)
                            , l2n_stem               & ! intent(in)
                            , negligible_nplant      & ! intent(in)
                            , agf_bs                 ! ! intent(in)
  use physiology_coms , only : N_plant_lim           ! ! intent(in)                     !JL 
  use decomp_coms   , only : f_labile                & ! intent(in)
                           , f_fast
  use ed_max_dims   , only : n_pft                  & ! intent(in)
                            , n_dbh                  ! ! intent(in)
   use ed_therm_lib  , only : calc_veg_hcap          & ! function
                            , update_veg_energy_cweh ! ! function
   use allometry     , only:  dbh2bl                 ! ! function
   use ed_misc_coms , only : current_time       ! ! intent(in)
   use nutrient_constants, only: organic_matter_density &
                               , crit_olayer_min        &
                               , crit_olayer_max        &
                               , aspen_seedling_mort    &
                               , blackspruce_seedling_mort
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: month
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer    :: cpoly
   type(sitetype)   , pointer    :: csite
   type(patchtype)  , pointer    :: cpatch
   integer                       :: ipy
   integer                       :: isi
   integer                       :: ipa
   integer                       :: ico
   integer                       :: ilu
   integer                       :: ipft
   integer                       :: update_month
   integer                       :: imonth
   real                          :: salloc
   real                          :: salloci
   real                          :: balive_in
   real                          :: bdead_in
   real                          :: hite_in
   real                          :: dbh_in
   real                          :: nplant_in
   real                          :: bstorage_in
   real                          :: agb_in
   real                          :: ba_in
   real                          :: f_bseeds
   real                          :: f_bdead
   real                          :: balive_mort_litter
   real                          :: bstorage_mort_litter
   real                          :: struct_litter
   real                          :: mort_litter
   real                          :: seed_litter
   real                          :: cb_act
   real                          :: cb_max
   real                          :: old_leaf_hcap
   real                          :: old_wood_hcap
   real				 :: ba_mort
   real				 :: total_ba_mort
   real                          :: organic_layer_depth !ATT
   real                          :: mort_fac !ATT
   real                          :: nstorage_mort_litter                               !JL
   real                          :: nstorage_in                                        !JL
   real                          :: deadwood_ratio !XXT
   real                          :: available_bstorage !XXT the carbon available for structural growth
   real				 :: available_nstorage
   real                          :: struct_nitrogen_demand
   real                          :: actual_available_bstorage

   !---------------------------------------------------------------------------------------!

   real :: n1, n2,nup,nloss,esum
   esum=0.
   ba_mort = 0.0
   total_ba_mort = 0.0

   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      !----- Initialization. --------------------------------------------------------------!
      cpoly%basal_area(:,:,:) = 0.0
      cpoly%agb(:,:,:)        = 0.0

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            ilu = csite%dist_type(ipa)
            n1=0.;n2=0.;nup=0.;nloss=0.

            !ATT 
            !calculate organic layer depth 
            organic_layer_depth = (csite%fast_soil_C(ipa) + csite%slow_soil_C(ipa))/organic_matter_density

            cohortloop: do ico = 1,cpatch%ncohorts
               !----- Assigning an alias for PFT type. ------------------------------------!
               ipft    = cpatch%pft(ico)
               !Calculate mort_fac. This modifies seedling mortality for aspen or black 
               !spruce trees based on organic layer depth. Seedling mortality starts
               !increaseing linearly at a depth crit_olayer_min and stops increasing 
               !at crit_olayer_max
               mort_fac = 1.0
               select case(ipft)
               case (21)
                  if (organic_layer_depth >= crit_olayer_min) then
                     mort_fac = (seedling_mortality(21) + (1 - seedling_mortality(21))*aspen_seedling_mort   &
                          *(min(organic_layer_depth,crit_olayer_max) - crit_olayer_min)/           &
                          (crit_olayer_max - crit_olayer_min)) &
                          /seedling_mortality(21)
                  end if
               case (22)
                  if (organic_layer_depth >= crit_olayer_min) then
                     mort_fac = (seedling_mortality(22) + (1 - seedling_mortality(22))*blackspruce_seedling_mort   &!changed from 0.75
                          *(min(organic_layer_depth,crit_olayer_max) - crit_olayer_min)/           &
                          (crit_olayer_max - crit_olayer_min)) &
                          /seedling_mortality(22)
                  end if
               end select
              !end ATT  

               salloc  = 1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico)
               salloci = 1.0 / salloc

               !----- Remember inputs in order to calculate increments later on. ----------!
               balive_in   = cpatch%balive(ico)
               bdead_in    = cpatch%bdead(ico)
               hite_in     = cpatch%hite(ico)
               dbh_in      = cpatch%dbh(ico)
               nplant_in   = cpatch%nplant(ico)
               bstorage_in = cpatch%bstorage(ico)
               nstorage_in = cpatch%nstorage(ico)                                      !JL!
               agb_in      = cpatch%agb(ico)
               ba_in       = cpatch%basarea(ico)

                                                                                        !JL!
              n1=n1+csite%area(ipa)*cpatch%nplant(ico)*(cpatch%balive(ico)/c2n_leaf(ipft)+ &
                    cpatch%nstorage(ico)+cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico)))   

               !---------------------------------------------------------------------------!
               !    Apply mortality, and do not allow nplant < negligible_nplant (such a   !
               ! sparse cohort is about to be terminated, anyway).                         !
               ! NB: monthly_dndt may be negative.                                         !
               !---------------------------------------------------------------------------!
               cpatch%monthly_dndt(ico) = max(cpatch%monthly_dndt(ico)                     &
                                             ,negligible_nplant(ipft) - cpatch%nplant(ico))
               !cpatch%monthly_dndt(ico) = 0.
               cpatch%nplant(ico)       = cpatch%nplant(ico) + cpatch%monthly_dndt(ico)
	       ba_mort = ba_mort + cpatch%monthly_dndt(ico) * cpatch%basarea(ico)* csite%area(ipa) * cpoly%area(isi) * 12

               !----- Calculate litter owing to mortality. --------------------------------!
               balive_mort_litter   = - cpatch%balive(ico)   * cpatch%monthly_dndt(ico)
               bstorage_mort_litter = - cpatch%bstorage(ico) * cpatch%monthly_dndt(ico)
               nstorage_mort_litter = - cpatch%nstorage(ico) * cpatch%monthly_dndt(ico)
               struct_litter        = - cpatch%bdead(ico)    * cpatch%monthly_dndt(ico)
               mort_litter          = balive_mort_litter + bstorage_mort_litter            &
                                    + struct_litter

               !----- Reset monthly_dndt. -------------------------------------------------!
               cpatch%monthly_dndt(ico) = 0.0

               available_bstorage = max(0.0, cpatch%bstorage(ico) - cpatch%bstorage_min(ico))  ! extra carbon
               available_nstorage = max(0.0, cpatch%nstorage(ico) - cpatch%nstorage_min(ico))  ! extra carbon
	           

               call plant_structural_allocation(cpatch%pft(ico),cpatch%hite(ico)           &!ATT
                    ,cgrid%lat(ipy),month                          &
                    ,cpatch%phenology_status(ico),f_bseeds,f_bdead)

               ! Determine the actual carbon available for growth
               deadwood_ratio = cpatch%bdead(ico) / (cpatch%bdead(ico) +cpatch%bsapwood(ico))   ! deadwood/sapwood ratio
               struct_nitrogen_demand = available_bstorage * f_bdead * deadwood_ratio /  c2n_stem(ipft) + &
                    available_bstorage * f_bdead * (1 - deadwood_ratio) / c2n_leaf(ipft) + &
                    available_bstorage * f_bseeds /	c2n_recruit(ipft)
               if (available_bstorage > 0. .and. struct_nitrogen_demand > 0.) then            
                    actual_available_bstorage = available_bstorage * &
                       min(1.0,available_nstorage / struct_nitrogen_demand)  ! consider nitrogen limitation
               else
                  actual_available_bstorage = 0.   ! otherwise actual_vailable_bstorage will have nan values...
               endif

               !----- Grow plants; bdead gets fraction f_bdead of bstorage. ---------------!
			   if(current_time%year > 2008 .or. .true. ) then  !start growing from 2009
	               cpatch%bdead(ico) = cpatch%bdead(ico) + f_bdead * actual_available_bstorage * deadwood_ratio
			   endif

               if(isnan(cpatch%bdead(ico))) then
                  print*,'ico',ico,'f_bdead',f_bdead,'available_nstorage',available_nstorage,'struct_nitrogen_demand',struct_nitrogen_demand,&
                       'dead_wood_ratio',deadwood_ratio,'f_bseeds',f_bseeds,'available_bstorage',available_bstorage
               endif
               !------ NPP allocation to wood and course roots in KgC /m2 -----------------!
               cpatch%today_NPPwood(ico) = agf_bs * f_bdead * actual_available_bstorage * deadwood_ratio         &
                                          * cpatch%nplant(ico)
               cpatch%today_NPPcroot(ico) = (1. - agf_bs) * f_bdead * actual_available_bstorage * deadwood_ratio &
                                          * cpatch%nplant(ico)
                                          
               !---------------------------------------------------------------------------!
               !      Calculate total seed production and seed litter.  The seed pool gets !
               ! a fraction f_bseeds of bstorage.                                          !
               !---------------------------------------------------------------------------!
               cpatch%bseeds(ico) = f_bseeds * actual_available_bstorage
               
               cpatch%today_NPPseeds(ico) = f_bseeds * actual_available_bstorage                &
                                          * cpatch%nplant(ico)
               
               seed_litter        = cpatch%bseeds(ico) * cpatch%nplant(ico)                &
                                  * mort_fac * seedling_mortality(ipft)
               if(isnan(seed_litter)) then
                  print*,'bseeds',cpatch%bseeds(ico),'nplant',cpatch%nplant(ico),'mort_fac',mort_fac,'seedling_mortality',seedling_mortality(ipft),'actual_available_bstorage',actual_available_bstorage,'available_bstorage',available_bstorage,'available_nstorage',available_nstorage,'struct_nitrogen_demand',struct_nitrogen_demand
               endif

               nloss=nloss+csite%area(ipa)*seed_litter/c2n_recruit(ipft)               

			   ! update bstorage and nstorage
			   if (current_time%year > 2008 .or. .true.) then
               cpatch%bstorage(ico) = cpatch%bstorage(ico)  -	 & 
					actual_available_bstorage * (f_bdead * deadwood_ratio + f_bseeds)
               cpatch%nstorage(ico) = max(cpatch%nstorage(ico)  -	 & 
					actual_available_bstorage * (f_bdead * deadwood_ratio / c2n_stem(ipft) + f_bseeds / c2n_recruit(ipft)),0.)
			   else
			   cpatch%bstorage(ico) = cpatch%bstorage(ico) - actual_available_bstorage
			   endif
if(cpatch%nstorage(ico) < 0.)then
print*,'7 Negative nstorage, nstorage',cpatch%nstorage(ico),'actual_available_bstorage',actual_available_bstorage
endif

if (ipa == 1 .and. ico == 1) then
	print*,'actual_available_bstorage',actual_available_bstorage,'f_bdead',f_bdead,'deadwood_ratio',deadwood_ratio
	print*,'N limitation',available_nstorage / struct_nitrogen_demand 
	print*,'phenology_status',cpatch%phenology_status(ico)
endif

               n2=n2+csite%area(ipa)*cpatch%nplant(ico)                                    &
                 * (cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))                           &
                  +cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))+cpatch%nstorage(ico)        &
                  +cpatch%bseeds(ico)/c2n_recruit(ipft))                                 !JL

               csite%slsc_in(ipa) = csite%slsc_in(ipa) + (f_labile(ipft) * balive_mort_litter &!ATT
                                 + bstorage_mort_litter + seed_litter) * (1 - f_fast(ipft))
               csite%slsn_in(ipa) = csite%slsn_in(ipa)                                       &
                                 + (f_labile(ipft) * balive_mort_litter / c2n_leaf(ipft)    &
                                 + seed_litter / c2n_recruit(ipft)                         &
                                 + nstorage_mort_litter) * (1 - f_fast(ipft)) 
               csite%fsc_in(ipa) = csite%fsc_in(ipa) + (f_labile(ipft) * balive_mort_litter &
                                 + bstorage_mort_litter + seed_litter ) *  f_fast(ipft)

               csite%fsn_in(ipa) = csite%fsn_in(ipa)                                       &
                                 + (f_labile(ipft) * balive_mort_litter / c2n_leaf(ipft)    &
                                 + seed_litter / c2n_recruit(ipft)                         &
                                 + nstorage_mort_litter) * f_fast(ipft)  
if(isnan(csite%slsn_in(ipa)))then
print*,'balive_mort_litter',balive_mort_litter,'seed_litter',seed_litter,'nstorage_mort_litter',nstorage_mort_litter,'bstorage_mort_litter',bstorage_mort_litter
endif

               csite%ssc_in(ipa) = csite%ssc_in(ipa) + struct_litter                       &
                                 + (1.0 - f_labile(ipft)) * balive_mort_litter
               csite%ssl_in(ipa) = csite%ssl_in(ipa)                                       &
                                 + ( (1.0 - f_labile(ipft)) * balive_mort_litter           &
                                    + struct_litter ) * l2n_stem / c2n_stem(cpatch%pft(ico))
               


               !----- Calculate the derived cohort properties. ----------------------------!
               call update_derived_cohort_props(cpatch,ico                                 &
                                                  ,csite%green_leaf_factor(ipft,ipa)       &
                                                  ,cpoly%lsl(isi))

				! Update bstorage_min and nstorage_min after updating dbh
               !---------------------------------------------------------------------------!
               ! MLO. We now update the heat capacity and the vegetation internal energy.  !
               !      Since no energy or water balance is done here, we simply update the  !
               !      energy in order to keep the same temperature and water as before.    !
               !      Internal energy is an extensive variable, we just account for the    !
               !      difference in the heat capacity to update it.                        !
               !---------------------------------------------------------------------------!
               old_leaf_hcap = cpatch%leaf_hcap(ico)
               old_wood_hcap = cpatch%wood_hcap(ico)
               call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico),cpatch%bsapwood(ico) &
                                 ,cpatch%nplant(ico),cpatch%pft(ico)                       &
                                 ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
               call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
               call is_resolvable(csite,ipa,ico,csite%green_leaf_factor(:,ipa))
               !---------------------------------------------------------------------------!


               !----- Update annual average carbon balances for mortality. ----------------!
               update_month = month - 1
               if(update_month == 0) update_month = 12
               cpatch%cb(update_month,ico)     = cpatch%cb(13,ico)
               cpatch%cb_max(update_month,ico) = cpatch%cb_max(13,ico)

               !----- If monthly files are written, save the current carbon balance. ------!
               if (associated(cpatch%mmean_cb)) then
                  cpatch%mmean_cb(ico)         = cpatch%cb(13,ico)
               end if

               !----- Reset the current month integrator. ---------------------------------!
               cpatch%cb(13,ico)               = 0.0
               cpatch%cb_max(13,ico)           = 0.0

               !----- Compute the relative carbon balance. --------------------------------!
               cb_act = 0.0
               cb_max = 0.0
               do imonth = 1,12
                  cb_act = cb_act + cpatch%cb(imonth,ico)
                  cb_max = cb_max + cpatch%cb_max(imonth,ico)
               end do
               if(cb_max > 0.0)then
                  cpatch%cbr_bar(ico) = cb_act / cb_max
               else
                  cpatch%cbr_bar(ico) = 0.0
               end if
               !---------------------------------------------------------------------------!


               !----- Update interesting output quantities. -------------------------------!
               call update_vital_rates(cpatch,ico,ilu,dbh_in,bdead_in,balive_in,hite_in    &
                                      ,bstorage_in,nplant_in,agb_in,ba_in,mort_litter      &
                                      ,csite%area(ipa),cpoly%basal_area(:,:,isi)           &
                                      ,cpoly%agb(:,:,isi),cpoly%basal_area_growth(:,:,isi) &
                                      ,cpoly%agb_growth(:,:,isi)                           &
                                      ,cpoly%basal_area_mort(:,:,isi)                      &
                                      ,cpoly%agb_mort(:,:,isi))


            end do cohortloop

            !----- Age the patch if this is not agriculture. ------------------------------!
            if (csite%dist_type(ipa) /= 1) csite%age(ipa) = csite%age(ipa) + 1.0/12.0
            esum=esum+n2-nup+nloss-n1
         end do patchloop
      end do siteloop
   end do polyloop
   total_ba_mort = sum(cpoly%basal_area_mort)
  ! print*,'ba_mort',ba_mort,'basal_area_mort',total_ba_mort
 !  print *,'organic_layer_depth', organic_layer_depth,'npatches',csite%npatches
  ! print*,'fast_soil_C', csite%fast_soil_C,'slow_soil_C', csite%slow_soil_C!ATT
  ! print*,'esum',esum
!   print*,'END STRUCT'

   return
end subroutine structural_growth
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute some structural growth variables without really         !
! updating the plant structure.  This should be used for test purposes only.               !
!------------------------------------------------------------------------------------------!
subroutine structural_growth_eq_0(cgrid, month)
   use ed_state_vars , only : edtype                 & ! structure
                            , polygontype            & ! structure
                            , sitetype               & ! structure
                            , patchtype              ! ! structure
   use pft_coms      , only : q                      & ! intent(in)
                            , qsw                    & ! intent(in)
                            , seedling_mortality     & ! intent(in)
                            , c2n_leaf               & ! intent(in)
                            , c2n_storage            & ! intent(in)
                            , c2n_recruit            & ! intent(in)
                            , c2n_stem               & ! intent(in)
                            , l2n_stem               ! ! intent(in)
   use decomp_coms   , only : f_labile               ! ! intent(in)
   use ed_max_dims   , only : n_pft                  & ! intent(in)
                            , n_dbh                  ! ! intent(in)
   use ed_therm_lib  , only : calc_veg_hcap          & ! function
                            , update_veg_energy_cweh ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: month
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer    :: cpoly
   type(sitetype)   , pointer    :: csite
   type(patchtype)  , pointer    :: cpatch
   integer                       :: ipy
   integer                       :: isi
   integer                       :: ipa
   integer                       :: ico
   integer                       :: ilu
   integer                       :: ipft
   integer                       :: update_month
   integer                       :: imonth
   real                          :: salloc
   real                          :: salloci
   real                          :: balive_in
   real                          :: bdead_in
   real                          :: hite_in
   real                          :: dbh_in
   real                          :: nplant_in
   real                          :: bstorage_in
   real                          :: agb_in
   real                          :: ba_in
   real                          :: f_bseeds
   real                          :: f_bdead
   real                          :: balive_mort_litter
   real                          :: bstorage_mort_litter
   real                          :: nstorage_mort_litter
   real                          :: struct_litter
   real                          :: mort_litter
   real                          :: seed_litter
   real                          :: net_seed_N_uptake
   real                          :: net_stem_N_uptake
   real                          :: cb_act
   real                          :: cb_max
   real                          :: old_leaf_hcap
   real                          :: old_wood_hcap
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      !----- Initialization. --------------------------------------------------------------!
      cpoly%basal_area(:,:,:) = 0.0
      cpoly%agb(:,:,:)        = 0.0

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            ilu = csite%dist_type(ipa)

            cohortloop: do ico = 1,cpatch%ncohorts
               !----- Assigning an alias for PFT type. ------------------------------------!
               ipft    = cpatch%pft(ico)

               salloc  = 1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico)
               salloci = 1.0 / salloc

               !----- Remember inputs in order to calculate increments later on. ----------!
               balive_in   = cpatch%balive(ico)
               bdead_in    = cpatch%bdead(ico)
               hite_in     = cpatch%hite(ico)
               dbh_in      = cpatch%dbh(ico)
               nplant_in   = cpatch%nplant(ico)
               bstorage_in = cpatch%bstorage(ico)
               agb_in      = cpatch%agb(ico)
               ba_in       = cpatch%basarea(ico)

               !----- Reset monthly_dndt. -------------------------------------------------!
               cpatch%monthly_dndt(ico) = 0.0

               !---------------------------------------------------------------------------!
               !      Rebalance the plant nitrogen uptake considering the actual alloc-    !
               ! ation to structural growth.  This is necessary because c2n_stem does not  !
               ! necessarily equal c2n_storage.                                            !
               !---------------------------------------------------------------------------!
               net_stem_N_uptake = (cpatch%bdead(ico) - bdead_in) * cpatch%nplant(ico)     &
                                 * ( 1.0 / c2n_stem(cpatch%pft(ico)) - 1.0 / c2n_storage)
               
               !---------------------------------------------------------------------------!
               !      Calculate total seed production and seed litter.  The seed pool gets !
               ! a fraction f_bseeds of bstorage.                                          !
               !---------------------------------------------------------------------------!
               cpatch%bseeds(ico) = 0.
               seed_litter        = 0.
               
               !---------------------------------------------------------------------------!
               !      Rebalance the plant nitrogen uptake considering the actual alloc-    !
               ! ation to seeds.  This is necessary because c2n_recruit does not have to   !
               ! be equal to c2n_storage.                                                  !
               !---------------------------------------------------------------------------!
               net_seed_N_uptake = cpatch%bseeds(ico) * cpatch%nplant(ico)                 &
                                 * (1.0 / c2n_recruit(ipft) - 1.0 / c2n_storage)

               !----- Finalize litter inputs. ---------------------------------------------!
               !Black Spruce litter inputs go into a slower decomposing pool
               select case (cpatch%pft(ico))!ATT
                  case (22)
                     csite%slsc_in(ipa) = csite%slsc_in(ipa) + f_labile(ipft) * balive_mort_litter &!ATT
                                 + bstorage_mort_litter + seed_litter
                     csite%slsn_in(ipa) = csite%slsn_in(ipa)                                       &
                                 + f_labile(ipft) * balive_mort_litter / c2n_leaf(ipft)    &!ATT
                                 + nstorage_mort_litter                                     &
                                 + seed_litter / c2n_recruit(ipft)
                  case default
                     csite%fsc_in(ipa) = csite%fsc_in(ipa) + f_labile(ipft) * balive_mort_litter &
                                 + bstorage_mort_litter + seed_litter 
                     csite%fsn_in(ipa) = csite%fsn_in(ipa)                                       &
                                 + f_labile(ipft) * balive_mort_litter / c2n_leaf(ipft)          &
                                 + nstorage_mort_litter                                          &
                                 + seed_litter / c2n_recruit(ipft)
                  end select!ATT
               !csite%fsc_in(ipa) = csite%fsc_in(ipa) + f_labile(ipft) * balive_mort_litter &
               !                  + bstorage_mort_litter + seed_litter
               !csite%fsn_in(ipa) = csite%fsn_in(ipa)                                       &
               !                  + f_labile(ipft) * balive_mort_litter / c2n_leaf(ipft)    &
               !                  + bstorage_mort_litter/ c2n_storage                       &
               !                  + seed_litter / c2n_recruit(ipft)
               csite%ssc_in(ipa) = csite%ssc_in(ipa) + struct_litter                       &
                                 + (1.0 - f_labile(ipft)) * balive_mort_litter
               csite%ssl_in(ipa) = csite%ssl_in(ipa) +                                     &
                                   ((1.0 - f_labile(ipft)) * balive_mort_litter            &
                                    + struct_litter ) * l2n_stem / c2n_stem(ipft)
               csite%total_plant_nitrogen_uptake(ipa) =                                    &
                      csite%total_plant_nitrogen_uptake(ipa) + net_seed_N_uptake           &
                    + net_stem_N_uptake

               !----- Update annual average carbon balances for mortality. ----------------!
               update_month = month - 1
               if(update_month == 0) update_month = 12
               cpatch%cb(update_month,ico)     = cpatch%cb(13,ico)
               cpatch%cb_max(update_month,ico) = cpatch%cb_max(13,ico)
               cpatch%cb(13,ico)               = 0.0
               cpatch%cb_max(13,ico)           = 0.0
               cb_act = 0.0
               cb_max = 0.0
               do imonth = 1,12
                  cb_act = cb_act + cpatch%cb(imonth,ico)
                  cb_max = cb_max + cpatch%cb_max(imonth,ico)
               end do
               if(cb_max > 0.0)then
                  cpatch%cbr_bar(ico) = cb_act / cb_max
               else
                  cpatch%cbr_bar(ico) = 0.0
               end if

               !----- Update interesting output quantities. -------------------------------!
               call update_vital_rates(cpatch,ico,ilu,dbh_in,bdead_in,balive_in,hite_in    &
                                      ,bstorage_in,nplant_in,agb_in,ba_in,mort_litter      &
                                      ,csite%area(ipa),cpoly%basal_area(:,:,isi)           &
                                      ,cpoly%agb(:,:,isi),cpoly%basal_area_growth(:,:,isi) &
                                      ,cpoly%agb_growth(:,:,isi)                           &
                                      ,cpoly%basal_area_mort(:,:,isi)                      &
                                      ,cpoly%agb_mort(:,:,isi))

            end do cohortloop

         end do patchloop
      end do siteloop
   end do polyloop


   return
end subroutine structural_growth_eq_0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will decide the partition of storage biomass into seeds and dead     !
! (structural) biomass.                                                                    !
!------------------------------------------------------------------------------------------!
subroutine plant_structural_allocation(ipft,hite,lat,month,phen_status,f_bseeds,f_bdead)
   use pft_coms      , only : phenology    & ! intent(in)
                            , repro_min_h  & ! intent(in)
                            , r_fract      & ! intent(in)
                            , st_fract     & ! intent(in)
                            , hgt_max      & ! intent(in)
                            , is_grass     ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer        , intent(in)  :: ipft
   integer        , intent(in)  :: month
   real           , intent(in)  :: hite
   real           , intent(in)  :: lat
   integer        , intent(in)  :: phen_status
   real           , intent(out) :: f_bseeds
   real           , intent(out) :: f_bdead
   !----- Local variables -----------------------------------------------------------------!
   logical                      :: late_spring
   !---------------------------------------------------------------------------------------!


   !----- Check whether this is late spring... --------------------------------------------!
   !late_spring = (lat >= 0.0 .and. month == 6) .or. (lat < 0.0 .and. month == 12) 
   late_spring = (lat >= 0.0 .and. (month == 6)) .or. (lat < 0.0 .and. month == 12)!ATT

   !---------------------------------------------------------------------------------------!
   !      Calculate fraction of bstorage going to bdead and reproduction.  First we must   !
   ! make sure that the plant should do something here.  A plant should not allocate any-  !
   ! thing to reproduction or growth if it is not the right time of year (for cold         !
   ! deciduous plants), or if the plants are actively dropping leaves or off allometry.    !
   !---------------------------------------------------------------------------------------!
   if ((phen_status == 0))then
      if (is_grass(ipft) .and. hite >= hgt_max(ipft)) then
         !---------------------------------------------------------------------------------!
         !    Grasses have reached the maximum height, stop growing in size and send       !
         ! everything to reproduction.                                                     !
         !---------------------------------------------------------------------------------!
         f_bseeds = 1.0 - st_fract(ipft)
      elseif (hite <= repro_min_h(ipft)) then
         !----- The plant is too short, invest as much as it can in growth. ---------------!
         f_bseeds = 0.0
      else
         !----- Plant is with a certain height, use prescribed reproduction rate. ---------!
         f_bseeds = r_fract(ipft) 
         !f_bseeds = 0.0
      end if
      f_bdead  = 1.0 - st_fract(ipft) - f_bseeds
      !f_bdead = 0.0
   else
      f_bdead  = 0.0
      f_bseeds = 0.0
   end if 
   !---------------------------------------------------------------------------------------!
          
   return
end subroutine plant_structural_allocation
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will assign values derived from the basic properties of a given      !
! cohort.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine update_derived_cohort_props(cpatch,ico,green_leaf_factor,lsl)

   use ed_state_vars , only : patchtype           ! ! structure
   use pft_coms      , only : phenology           & ! intent(in)
                            , q                   & ! intent(in)
                            , qsw                 &
							, c2n_leaf			  &
							, SLA! ! intent(in)
   use allometry     , only : bd2dbh              & ! function
                            , dbh2h               & ! function
                            , dbh2bl              & ! function
                            , dbh2krdepth         & ! function
                            , ed_biomass          & ! function
                            , area_indices        ! ! subroutine
   use consts_coms   , only : pio4                ! ! intent(in)
   use phenology_coms, only: theta_crit   ! ! intent(in)
   
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(patchtype), target     :: cpatch
   integer        , intent(in) :: ico
   real           , intent(in) :: green_leaf_factor
   integer        , intent(in) :: lsl
   !----- Local variables -----------------------------------------------------------------!
   real :: bl
   real :: bl_max
   !---------------------------------------------------------------------------------------!

   !----- Gett DBH and height from structural biomass. ------------------------------------!
   cpatch%dbh(ico)  = bd2dbh(cpatch%pft(ico), cpatch%bdead(ico))
   cpatch%hite(ico) = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))
   
   ! ---- XXT Update the sla according to Height   ---------- 
   ! From Loyld et al. 2010 Optimisation of photosynthetic carbon
   ! gain and within-canopy gradients of associated foliar traits
   ! for Amazon forest trees. Biogesciences
   select case (cpatch%pft(ico))
     case (21,22,2,3,4)
        cpatch%hite_coef(ico) = 1.0
        cpatch%sla(ico) = SLA(cpatch%pft(ico))

     case default
         cpatch%hite_coef(ico) = max(1.0,min(1.40, & 
					 1.0 + max(0.0,min((cpatch%hite(ico)-5.)/15.0,1.0)) * 0.40))
         cpatch%sla(ico) = SLA(cpatch%pft(ico)) *  &
                max(0.65, min(1.00, &
                        1.00 - max(0.0,min((cpatch%hite(ico) - 5.)/15.0,1.0)) * 0.35))
                
   end select

   ! update bstorage_min and nstorage_min
   cpatch%bstorage_min(ico) = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico)) * (1 + q(cpatch%pft(ico)))
   if (cpatch%bstorage_min(ico) < 1.0e-8) then
      cpatch%bstorage_min(ico) = 0.
   endif
   cpatch%nstorage_min(ico) = cpatch%bstorage_min(ico) / c2n_leaf(cpatch%pft(ico))

   !----- Check the phenology status and whether it needs to change. ----------------------!
!   select case (cpatch%phenology_status(ico))
!   case (0,1)
!
!      select case (phenology(cpatch%pft(ico)))
!      case (4)
!         cpatch%elongf(ico)  = max(0.0,min (1.0, cpatch%paw_avg(ico)/theta_crit))
!      case default
!         cpatch%elongf(ico)  = 1.0
!      end select
!
!      bl_max = dbh2bl(cpatch%dbh(ico),cpatch%pft(ico))                                     &
!             * green_leaf_factor * cpatch%elongf(ico)
!      !------------------------------------------------------------------------------------!
!      !     If LEAF biomass is not the maximum, set it to 1 (leaves partially flushed),    !
!      ! otherwise, set it to 0 (leaves are fully flushed).                                 !
!      !------------------------------------------------------------------------------------!
!      if (cpatch%bleaf(ico) < bl_max) then
!         cpatch%phenology_status(ico) = 1
!      else
!         cpatch%phenology_status(ico) = 0
!      end if
!
!   end select
      
   !----- Update LAI, WPA, and WAI --------------------------------------------------------!
   call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bdead(ico)                &
              ,cpatch%balive(ico),cpatch%dbh(ico), cpatch%hite(ico),cpatch%pft(ico)        &
              ,cpatch%sla(ico),cpatch%lai(ico),cpatch%wpa(ico),cpatch%wai(ico)             &
              ,cpatch%crown_area(ico),cpatch%bsapwood(ico))

   !----- Finding the new basal area and above-ground biomass. ----------------------------!
   cpatch%basarea(ico) = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)                
   cpatch%agb(ico)     = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico),cpatch%bleaf(ico) &
                                   ,cpatch%pft(ico),cpatch%hite(ico) ,cpatch%bstorage(ico) &
                                   ,cpatch%bsapwood(ico))

   !----- Update rooting depth ------------------------------------------------------------!
   cpatch%krdepth(ico) = dbh2krdepth(cpatch%hite(ico),cpatch%dbh(ico),cpatch%pft(ico),lsl)
   
   return
end subroutine update_derived_cohort_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute the growth and mortality rates.                          !
!------------------------------------------------------------------------------------------!
subroutine update_vital_rates(cpatch,ico,ilu,dbh_in,bdead_in,balive_in,hite_in,bstorage_in &
                             ,nplant_in,agb_in,ba_in,mort_litter,area,basal_area,agb       &
                             ,basal_area_growth,agb_growth,basal_area_mort,agb_mort)
   
   use ed_state_vars , only : patchtype    ! ! structure
   use ed_max_dims   , only : n_pft        & ! intent(in)
                            , n_dbh        & ! intent(in)
                            , n_dist_types ! ! intent(in)
   use ed_misc_coms  , only : ddbhi        ! ! intent(in)
   use consts_coms   , only : pio4         ! ! intent(in)
   use pft_coms      , only : agf_bs       & ! intent(in)
                            , q            & ! intent(in)
                            , qsw          ! ! intent(in)
   use allometry     , only : ed_biomass   ! ! function
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(patchtype)              , target        :: cpatch
   real                         , intent(in)    :: dbh_in
   real                         , intent(in)    :: bdead_in
   real                         , intent(in)    :: balive_in
   real                         , intent(in)    :: hite_in
   real                         , intent(in)    :: bstorage_in
   real                         , intent(in)    :: nplant_in
   real                         , intent(in)    :: agb_in
   real                         , intent(in)    :: ba_in
   real                         , intent(in)    :: mort_litter
   real                         , intent(in)    :: area
   integer                      , intent(in)    :: ico
   integer                      , intent(in)    :: ilu
   real, dimension(n_pft, n_dbh), intent(inout) :: basal_area
   real, dimension(n_pft, n_dbh), intent(inout) :: agb
   real, dimension(n_pft, n_dbh), intent(inout) :: basal_area_growth
   real, dimension(n_pft, n_dbh), intent(inout) :: agb_growth
   real, dimension(n_pft, n_dbh), intent(inout) :: basal_area_mort
   real, dimension(n_pft, n_dbh), intent(inout) :: agb_mort
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: ipft
   integer                                      :: idbh
   !---------------------------------------------------------------------------------------!


   !----- Make the alias for PFT type. ----------------------------------------------------!
   ipft = cpatch%pft(ico)

   !----- Find the DBH bin. ---------------------------------------------------------------!
   idbh = max(1,min(n_dbh,ceiling(dbh_in*ddbhi)))

   !----- Find the new basal area and above-ground biomass. -------------------------------!
   cpatch%basarea(ico)    = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)
   cpatch%agb(ico)        = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)                &
                                      ,cpatch%bleaf(ico),cpatch%pft(ico)                   &
                                      ,cpatch%hite(ico) ,cpatch%bstorage(ico)              &
                                      ,cpatch%bsapwood(ico) ) 

   !---------------------------------------------------------------------------------------!
   !     Change the agb growth to kgC/plant/year, basal area to cm2/plant/year, and DBH    !
   ! growth to cm/year.                                                                    !
   !---------------------------------------------------------------------------------------!
   cpatch%dagb_dt(ico)    = (cpatch%agb(ico)     - agb_in ) * 12.0
   cpatch%dba_dt(ico)     = (cpatch%basarea(ico) - ba_in  ) * 12.0
   cpatch%ddbh_dt(ico)    = (cpatch%dbh(ico)     - dbh_in ) * 12.0

   !---------------------------------------------------------------------------------------!
   !     These are polygon-level variable, so they are done in kgC/m2.  Update the current !
   ! basal area and above-ground biomass.                                                  !
   !---------------------------------------------------------------------------------------!
   basal_area(ipft, idbh) = basal_area(ipft, idbh)                                         &
                          + area * cpatch%nplant(ico) * cpatch%basarea(ico)
   agb(ipft, idbh)        = agb(ipft, idbh)                                                &
                          + area * cpatch%nplant(ico) * cpatch%agb(ico)

   !---------------------------------------------------------------------------------------!
   !    The growth and mortality census are applied only on those cohorts present on the   !
   ! first census.                                                                         !
   !---------------------------------------------------------------------------------------!
   if (cpatch%first_census(ico) /= 1) return

   !---------------------------------------------------------------------------------------!
   !   Computed for plants alive both at past census and current census.  These will be    !
   ! given in cm2/m2/yr and kgC/m2/yr, respectively.                                       !
   !---------------------------------------------------------------------------------------!
   basal_area_growth(ipft,idbh) = basal_area_growth(ipft,idbh)                             &
                                + area * cpatch%nplant(ico) * pio4                         &
                                * (cpatch%dbh(ico) * cpatch%dbh(ico) - dbh_in * dbh_in)    &
                                * 12.0
   agb_growth(ipft,idbh)        = agb_growth(ipft,idbh)                                    &
                                + area * cpatch%nplant(ico)                                &
                                * (cpatch%agb(ico) - agb_in)                               &
                                * 12.0 

   !---------------------------------------------------------------------------------------!
   !    Computed for plants alive at past census but dead at current census.  These        !
   ! variables are also given in cm2/m2/yr and kgC/m2/yr, respectively.                    !
   !---------------------------------------------------------------------------------------!
   basal_area_mort(ipft,idbh) = basal_area_mort(ipft,idbh)                                 &
                              + area * (nplant_in - cpatch%nplant(ico)) * ba_in * 12.0

   !----- Calculation based on mort_litter includes TOTAL biomass, not AGB [[mcd]]. -------!
   agb_mort(ipft,idbh)        = agb_mort(ipft,idbh)                                        &
                              + area * (nplant_in - cpatch%nplant(ico)) * agb_in * 12.0

   return
end subroutine update_vital_rates
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will display the carbon and nitrogen budgets on screen.               !
!------------------------------------------------------------------------------------------!
subroutine print_C_and_N_budgets(cgrid)
   use ed_state_vars , only : edtype       ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)                 , target    :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   integer                                  :: ipy
   real                                     :: soil_C
   real                                     :: soil_N
   real                                     :: veg_C
   real                                     :: veg_N
   !----- Local constants -----------------------------------------------------------------!
   logical                      , parameter :: print_on = .false.
   !---------------------------------------------------------------------------------------!
   


   do ipy = 1,cgrid%npolygons

      call compute_C_and_N_storage(cgrid,ipy, soil_C, soil_N, veg_C, veg_N)

      if (print_on) then
         
         write (unit=*,fmt='(a)') '--------------------------------------------------------'
         write (unit=*,fmt='(a,1x,i6)')     'POLYGON           =',ipy
         write (unit=*,fmt='(a,1x,f9.3)')   'LON               =',cgrid%lon(ipy)
         write (unit=*,fmt='(a,1x,f9.3)')   'LAT               =',cgrid%lat(ipy)
         write (unit=*,fmt='(a,1x,es12.5)') 'C_INITIAL_STORAGE ='                          &
                                                        ,cgrid%cbudget_initialstorage(ipy)
         write (unit=*,fmt='(a,1x,es12.5)') 'SOIL_C+VEG_C-NCEP ='                          &
                                                       ,soil_C+veg_C-cgrid%cbudget_nep(ipy)
         write (unit=*,fmt='(a,1x,es12.5)') 'SOIL_C+VEG_C      = ',soil_C + veg_C
         write (unit=*,fmt='(a,1x,es12.5)') 'NEP               = ',cgrid%cbudget_nep(ipy)
         write (unit=*,fmt='(a,1x,es12.5)') 'N_INITIAL_STORAGE ='                          &
                                                        ,cgrid%nbudget_initialstorage(ipy)
         write (unit=*,fmt='(a,1x,es12.5)') 'SOIL_N+VEG_N      = ',soil_N + veg_N
         write (unit=*,fmt='(a)') '--------------------------------------------------------'
         write (unit=*,fmt='(a)') ' '
      end if
   end do
   return
end subroutine print_C_and_N_budgets
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute the carbon and nitrogen pools.                           !
!------------------------------------------------------------------------------------------!
subroutine compute_C_and_N_storage(cgrid,ipy, soil_C, soil_N, veg_C, veg_N)

   use ed_state_vars , only : edtype         & ! structure
                            , polygontype    & ! structure
                            , sitetype       & ! structure
                            , patchtype      ! ! structure
   use ed_max_dims      , only : n_pft          ! ! intent(in)
   use pft_coms      , only : include_pft    & ! intent(in)
                            , c2n_recruit    & ! intent(in)
                            , c2n_stem       & ! intent(in)
                            , c2n_leaf       & ! intent(in)
                            , c2n_storage    & ! intent(in)
                            , c2n_slow       & ! intent(in)
                            , c2n_structural ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target      :: cgrid
   integer           , intent(in)  :: ipy
   real              , intent(out) :: soil_C
   real              , intent(out) :: soil_N
   real              , intent(out) :: veg_C
   real              , intent(out) :: veg_N
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer     :: cpoly
   type(sitetype)    , pointer     :: csite
   type(patchtype)   , pointer     :: cpatch
   integer                         :: isi
   integer                         :: ipa
   integer                         :: ico
   integer                         :: ipft
   real(kind=8)                    :: area_factor
   real(kind=8)                    :: this_carbon
   real(kind=8)                    :: this_nitrogen
   real(kind=8)                    :: soil_C8
   real(kind=8)                    :: soil_N8
   real(kind=8)                    :: veg_C8
   real(kind=8)                    :: veg_N8
   !----- Local constants -----------------------------------------------------------------!
   real(kind=8)      , parameter   :: almostnothing=1.d-30
   !----- External functions. -------------------------------------------------------------!
   real              , external    :: sngloff
   !---------------------------------------------------------------------------------------!

   !----- Initialize C and N pools. -------------------------------------------------------!
   soil_C8 = 0.0d0
   soil_N8 = 0.0d0
   veg_C8  = 0.0d0
   veg_N8  = 0.0d0

   cpoly => cgrid%polygon(ipy)
   siteloop: do isi = 1,cpoly%nsites

      csite => cpoly%site(isi)
      patchloop: do ipa = 1,csite%npatches
         cpatch => csite%patch(ipa)

         !----- Site area times patch area. -----------------------------------------------!
         area_factor   = dble(cpoly%area(isi)) * dble(csite%area(ipa))

         !----- Find carbon and nitrogen soil pools for this patch. -----------------------!
         this_carbon   = dble(csite%fast_soil_C(ipa)) + dble(csite%slow_soil_C(ipa))       &
                       + dble(csite%structural_soil_C(ipa))
         this_nitrogen = dble(csite%fast_soil_N(ipa))                                      &
                       + dble(csite%mineralized_soil_N(ipa))                               &
                       + dble(csite%slow_soil_N(ipa))                                      &!ATT
                       + dble(csite%structural_soil_C(ipa)) / dble(c2n_structural)

         !----- Add to the full counter. --------------------------------------------------!
         soil_C8 = soil_C8 + area_factor * this_carbon
         soil_N8 = soil_N8 + area_factor * this_nitrogen

         !----- Loop over PFT so we account for veg carbon/nitrogen in repro arrays. ------!
         pftloop: do ipft = 1, n_pft
            if (include_pft(ipft)) then
               veg_C8 = veg_C8 + dble(csite%repro(ipft,ipa)) * area_factor
               veg_N8 = veg_N8 + dble(csite%repro(ipft,ipa))                               &
                               / dble(c2n_recruit(ipft)) * area_factor
            end if
         end do pftloop
         
         cohortloop: do ico = 1,cpatch%ncohorts
            
            !----- Get the carbon and nitrogen in vegetation. -----------------------------!
            veg_C8 = veg_C8 + area_factor                                                  &
                            * ( dble(cpatch%balive(ico)) + dble(cpatch%bdead(ico))         &
                              + dble(cpatch%bstorage(ico)) ) * dble(cpatch%nplant(ico))           
!JL!
           veg_N8 = veg_N8 + area_factor                                                   &
                            * ( dble(cpatch%balive(ico)) / dble(c2n_leaf(cpatch%pft(ico))) &
                              + dble(cpatch%bdead(ico)) / dble(c2n_stem(cpatch%pft(ico)))  &
                              + dble(cpatch%nstorage(ico)))                                &
                            * dble(cpatch%nplant(ico))
!JL!
         end do cohortloop
      end do patchloop
   end do siteloop

   soil_C = sngloff(soil_C8,almostnothing)
   soil_N = sngloff(soil_N8,almostnothing)
   veg_C  = sngloff(veg_C8 ,almostnothing)
   veg_N  = sngloff(veg_N8 ,almostnothing)

   return
end subroutine compute_C_and_N_storage
!==========================================================================================!
!==========================================================================================!


