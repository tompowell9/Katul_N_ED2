!==========================================================================================!
!==========================================================================================!
!     This module contains the wrappers for the Runge-Kutta integration scheme.            !
!==========================================================================================!
!==========================================================================================!
module rk4_driver

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      Main driver of short-time scale dynamics of the Runge-Kutta integrator           !
   !      for the land surface model.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_timestep(cgrid,ifm)
      use rk4_coms               , only : integration_vars     & ! structure
                                        , rk4patchtype         & ! structure
                                        , zero_rk4_patch       & ! subroutine
                                        , zero_rk4_cohort      & ! subroutine
                                        , integration_buff     & ! intent(out)
                                        , rk4site              ! ! intent(out)
      use ed_state_vars          , only : edtype               & ! structure
                                        , polygontype          & ! structure
                                        , sitetype             & ! structure
                                        , patchtype            ! ! structure
      use met_driver_coms        , only : met_driv_state       ! ! structure
      use grid_coms              , only : nzg                  & ! intent(in)
                                        , nzs                  ! ! intent(in)
      use canopy_struct_dynamics , only : canopy_turbulence8   ! ! subroutine
      use ed_misc_coms           , only : current_time         ! ! intent(in)
      implicit none

      !----------- Use MPI timing calls, need declarations --------------------------------!
      include 'mpif.h'
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target      :: cgrid
      integer                   , intent (in) :: ifm
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
      type(patchtype)           , pointer     :: cpatch
      type(met_driv_state)      , pointer     :: cmet
      integer                                 :: ipy
      integer                                 :: isi
      integer                                 :: ipa
      integer                                 :: iun
      integer                                 :: nsteps
      real                                    :: wcurr_loss2atm
      real                                    :: ecurr_loss2atm
      real                                    :: co2curr_loss2atm
      real                                    :: wcurr_loss2drainage
      real                                    :: ecurr_loss2drainage
      real                                    :: wcurr_loss2runoff
      real                                    :: ecurr_loss2runoff
      real                                    :: old_can_theiv
      real                                    :: old_can_shv
      real                                    :: old_can_co2
      real                                    :: old_can_rhos
      real                                    :: old_can_temp
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!

      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               
               !----- Reset all buffers to zero, as a safety measure. ---------------------!
               call zero_rk4_patch(integration_buff%initp)
               call zero_rk4_patch(integration_buff%yscal)
               call zero_rk4_patch(integration_buff%y)
               call zero_rk4_patch(integration_buff%dydx)
               call zero_rk4_patch(integration_buff%yerr)
               call zero_rk4_patch(integration_buff%ytemp)
               call zero_rk4_patch(integration_buff%ak2)
               call zero_rk4_patch(integration_buff%ak3)
               call zero_rk4_patch(integration_buff%ak4)
               call zero_rk4_patch(integration_buff%ak5)
               call zero_rk4_patch(integration_buff%ak6)
               call zero_rk4_patch(integration_buff%ak7)
               call zero_rk4_cohort(integration_buff%initp)
               call zero_rk4_cohort(integration_buff%yscal)
               call zero_rk4_cohort(integration_buff%y)
               call zero_rk4_cohort(integration_buff%dydx)
               call zero_rk4_cohort(integration_buff%yerr)
               call zero_rk4_cohort(integration_buff%ytemp)
               call zero_rk4_cohort(integration_buff%ak2)
               call zero_rk4_cohort(integration_buff%ak3)
               call zero_rk4_cohort(integration_buff%ak4)
               call zero_rk4_cohort(integration_buff%ak5)
               call zero_rk4_cohort(integration_buff%ak6)
               call zero_rk4_cohort(integration_buff%ak7)

               !----- Get velocity for aerodynamic resistance. ----------------------------!
               if (csite%can_theta(ipa) < cmet%atm_theta) then
                  cmet%vels = cmet%vels_stab
               else
                  cmet%vels = cmet%vels_unstab
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Update roughness and canopy depth.                                     !
               !---------------------------------------------------------------------------!
               call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs)
               call update_patch_derived_props(csite,cpoly%lsl(isi),cmet%prss,ipa)
               !---------------------------------------------------------------------------!


               !----- Save the previous thermodynamic state. ------------------------------!
               old_can_theiv    = csite%can_theiv(ipa)
               old_can_shv      = csite%can_shv(ipa)
               old_can_co2      = csite%can_co2(ipa)
               old_can_rhos     = csite%can_rhos(ipa)
               old_can_temp     = csite%can_temp(ipa)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Copy the meteorological variables to the rk4site structure.            !
               !---------------------------------------------------------------------------!
               call copy_met_2_rk4site(nzg,cmet%vels,cmet%atm_theiv,cmet%atm_theta         &
                                      ,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2,cmet%geoht   &
                                      ,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg          &
                                      ,cmet%prss,cmet%rshort,cmet%rlong,cmet%geoht         &
                                      ,cpoly%lsl(isi),csite%ntext_soil(:,ipa)              &!ATT
                                      ,csite%green_leaf_factor(:,ipa)                      &
                                      ,cgrid%lon(ipy),cgrid%lat(ipy),cgrid%cosz(ipy))

               !----- Compute current storage terms. --------------------------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)




               !----- Get photosynthesis, stomatal conductance, and transpiration. --------!
               call canopy_photosynthesis(csite,cmet,nzg,ipa,cpoly%lsl(isi)                &
                                         ,csite%ntext_soil(:,ipa)                          &!ATT
                                         ,csite%leaf_aging_factor(:,ipa)                   &
                                         ,csite%green_leaf_factor(:,ipa))

               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init(csite,ipa,integration_buff%initp)

               !----- Compute root and heterotrophic respiration. -------------------------!
               call soil_respiration(csite,ipa,nzg,csite%ntext_soil(:,ipa))                !ATT

               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init_carbon(csite,ipa,integration_buff%initp)

               !---------------------------------------------------------------------------!
               !    This is the driver for the integration process...                      !
               !---------------------------------------------------------------------------!
               call integrate_patch_rk4(csite,integration_buff%initp,ipa,wcurr_loss2atm    &
                                       ,ecurr_loss2atm,co2curr_loss2atm                    &
                                       ,wcurr_loss2drainage,ecurr_loss2drainage            &
                                       ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps,cgrid%lon(ipy))

               !----- Add the number of steps into the step counter. ----------------------!
               cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)

               !---------------------------------------------------------------------------!
               !    Update the minimum monthly temperature, based on canopy temperature.   !
               !---------------------------------------------------------------------------!
               if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                  cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
               end if
               
               !---------------------------------------------------------------------------!
               !     Compute the residuals.                                                !
               !---------------------------------------------------------------------------!
               call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa           &
                                  ,wcurr_loss2atm,ecurr_loss2atm,co2curr_loss2atm          &
                                  ,wcurr_loss2drainage,ecurr_loss2drainage                 &
                                  ,wcurr_loss2runoff,ecurr_loss2runoff,cpoly%area(isi)     &
                                  ,cgrid%cbudget_nep(ipy),old_can_theiv,old_can_shv        &
                                  ,old_can_co2,old_can_rhos,old_can_temp)

            end do patchloop
         end do siteloop

      end do polygonloop

      return
   end subroutine rk4_timestep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_rk4(csite,initp,ipa,wcurr_loss2atm,ecurr_loss2atm            &
                                 ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage &
                                 ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps,lon)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use ed_misc_coms    , only : dtlsm                ! ! intent(in)
      use soil_coms       , only : soil_rough           & ! intent(in)
                                 , snow_rough           ! ! intent(in)
      use canopy_air_coms , only : exar8                ! ! intent(in)
      use consts_coms     , only : vonk8                & ! intent(in)
                                 , cp8                  & ! intent(in)
                                 , cpi8                 ! ! intent(in)
      use rk4_coms        , only : integration_vars     & ! structure
                                 , rk4patchtype         & ! structure
                                 , rk4site              & ! intent(inout)
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4                & ! intent(inout)
                                 , dtrk4i               ! ! intent(inout)
	  use physiology_coms , only : h2o_plant_lim
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      type(rk4patchtype)    , target      :: initp
      integer               , intent(in)  :: ipa
      real                  , intent(out) :: wcurr_loss2atm
      real                  , intent(out) :: ecurr_loss2atm
      real                  , intent(out) :: co2curr_loss2atm
      real                  , intent(out) :: wcurr_loss2drainage
      real                  , intent(out) :: ecurr_loss2drainage
      real                  , intent(out) :: wcurr_loss2runoff
      real                  , intent(out) :: ecurr_loss2runoff
      integer               , intent(out) :: nsteps
	  real					, intent(in)  :: lon
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                          :: hbeg
      !----- Locally saved variable -------------------------------------------------------!
      logical                  , save       :: first_time=.true.
      !------------------------------------------------------------------------------------!

      !----- Assign some constants which will remain the same throughout the run. ---------!
      if (first_time) then
         first_time = .false.
         tbeg   = 0.d0
         tend   = dble(dtlsm)
         dtrk4  = tend - tbeg
         dtrk4i = 1.d0/dtrk4
      end if

      !------------------------------------------------------------------------------------!
      !      Initial step size.  Experience has shown that giving this too large a value   !
      ! causes the integrator to fail (e.g., soil layers become supersaturated).           !
      !------------------------------------------------------------------------------------!
      hbeg = dble(csite%htry(ipa))


      !------------------------------------------------------------------------------------!
      !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm, !
      ! so they must be zeroed at each call.                                               !
      !------------------------------------------------------------------------------------!
      initp%upwp = 0.d0
      initp%tpwp = 0.d0
      initp%qpwp = 0.d0
      initp%cpwp = 0.d0
      initp%wpwp = 0.d0

	  initp%avg_nutrient_layer_drainage = 0.0d0
      !----- Go into the ODE integrator. --------------------------------------------------!
      call odeint(hbeg,csite,ipa,nsteps)

      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
      initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
      initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
      initp%cpwp = initp%can_rhos * initp%cpwp * dtrk4i
      initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i
      
	  if(h2o_plant_lim == 3) then
		! New water limitaiotn scheme
		! update plant and soil water conditions
		  call update_water_potential(csite,initp,ipa,lon)
	  endif
	  
      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp(tend-tbeg,initp,csite,ipa,wcurr_loss2atm,ecurr_loss2atm            &
                       ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage           &
                       ,wcurr_loss2runoff,ecurr_loss2runoff)

      return
   end subroutine integrate_patch_rk4
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp(hdid,initp,csite,ipa,wbudget_loss2atm,ebudget_loss2atm          &
                          ,co2budget_loss2atm,wbudget_loss2drainage,ebudget_loss2drainage  &
                          ,wbudget_loss2runoff,ebudget_loss2runoff)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4site              & ! intent(in)
                                      , rk4min_veg_temp      & ! intent(in)
                                      , rk4max_veg_temp      & ! intent(in)
                                      , tiny_offset          & ! intent(in) 
                                      , checkbudget          ! ! intent(in)
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            & ! structure
                                      , edgrid_g             ! ! structure
      use consts_coms          , only : day_sec              & ! intent(in)
                                      , t3ple8               ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     & ! intent(in)
                                      , dtlsm                ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , slz8                 ! ! intent(in)
      use grid_coms            , only : nzg                  & ! intent(in)
                                      , nzs                  ! ! intent(in)
      use therm_lib            , only : qwtk                 & ! subroutine
                                      , rslif                ! ! function
      use canopy_air_coms      , only : i_blyr_condct        ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target      :: initp
      type(sitetype)    , target      :: csite
      real(kind=8)      , intent(in)  :: hdid
      integer           , intent(in)  :: ipa
      real              , intent(out) :: wbudget_loss2atm
      real              , intent(out) :: ebudget_loss2atm
      real              , intent(out) :: co2budget_loss2atm
      real              , intent(out) :: wbudget_loss2drainage
      real              , intent(out) :: ebudget_loss2drainage
      real              , intent(out) :: wbudget_loss2runoff
      real              , intent(out) :: ebudget_loss2runoff
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      integer                         :: mould
      integer                         :: ico
      integer                         :: k
      integer                         :: kclosest
      integer                         :: ksn
      integer                         :: nsoil
      integer                         :: nlsw1
      real(kind=8)                    :: tmp_energy
!      real(kind=8), dimension(nzg)                    :: available_water
      !----- Local contants ---------------------------------------------------------------!
      real        , parameter         :: tendays_sec=10.*day_sec
      !----- External function ------------------------------------------------------------!
      real        , external          :: sngloff
      !------------------------------------------------------------------------------------!
      real(kind=8), dimension(nzg) :: relative_porosity, weighted_relative_porosity,  &
           cumul_wei_rel_porosity, avg_rel_porosity



      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.  All floating point variables are converted back  !
      ! to single precision.                                                               !
      !------------------------------------------------------------------------------------!
      csite%can_theiv(ipa)        = sngloff(initp%can_theiv       ,tiny_offset)
      csite%can_theta(ipa)        = sngloff(initp%can_theta       ,tiny_offset)
      csite%can_prss(ipa)         = sngloff(initp%can_prss        ,tiny_offset)
      csite%can_temp(ipa)         = sngloff(initp%can_temp        ,tiny_offset)
      csite%can_shv(ipa)          = sngloff(initp%can_shv         ,tiny_offset)
      csite%can_co2(ipa)          = sngloff(initp%can_co2         ,tiny_offset)
      csite%can_rhos(ipa)         = sngloff(initp%can_rhos        ,tiny_offset)
      csite%can_depth(ipa)        = sngloff(initp%can_depth       ,tiny_offset)
      csite%veg_displace(ipa)     = sngloff(initp%veg_displace    ,tiny_offset)
      csite%rough(ipa)            = sngloff(initp%rough           ,tiny_offset)
      csite%snowfac(ipa)          = sngloff(initp%snowfac         ,tiny_offset)
      csite%total_sfcw_depth(ipa) = sngloff(initp%total_sfcw_depth,tiny_offset)

      csite%ggbare(ipa)           = sngloff(initp%ggbare          ,tiny_offset)
      csite%ggveg (ipa)           = sngloff(initp%ggveg           ,tiny_offset)
      csite%ggnet (ipa)           = sngloff(initp%ggnet           ,tiny_offset)

      csite%ustar (ipa)           = sngloff(initp%ustar           ,tiny_offset)
      csite%tstar (ipa)           = sngloff(initp%tstar           ,tiny_offset)
      csite%qstar (ipa)           = sngloff(initp%qstar           ,tiny_offset)
      csite%cstar (ipa)           = sngloff(initp%cstar           ,tiny_offset)

      csite%zeta  (ipa)           = sngloff(initp%zeta            ,tiny_offset)
      csite%ribulk(ipa)           = sngloff(initp%ribulk          ,tiny_offset)

      csite%upwp  (ipa)           = sngloff(initp%upwp            ,tiny_offset)
      csite%wpwp  (ipa)           = sngloff(initp%wpwp            ,tiny_offset)
      csite%tpwp  (ipa)           = sngloff(initp%tpwp            ,tiny_offset)
      csite%qpwp  (ipa)           = sngloff(initp%qpwp            ,tiny_offset)
      csite%cpwp  (ipa)           = sngloff(initp%cpwp            ,tiny_offset)
         

	  csite%avg_nutrient_layer_drainage        (ipa) = sngloff(initp%avg_nutrient_layer_drainage       ,tiny_offset)
	  csite%dmean_nutrient_layer_drainage        (ipa) = csite%dmean_nutrient_layer_drainage(ipa) + &
	  													csite%avg_nutrient_layer_drainage(ipa)
	  ! drainage should be updated anyway
      !------------------------------------------------------------------------------------!
      !    These variables are fast scale fluxes, and they may not be allocated, so just   !
      ! check this before copying.                                                         !
      !------------------------------------------------------------------------------------!
      if (fast_diagnostics) then
         csite%avg_vapor_lc        (ipa) = sngloff(initp%avg_vapor_lc       ,tiny_offset)
         csite%avg_vapor_wc        (ipa) = sngloff(initp%avg_vapor_wc       ,tiny_offset)
         csite%avg_dew_cg          (ipa) = sngloff(initp%avg_dew_cg         ,tiny_offset)
         csite%avg_vapor_gc        (ipa) = sngloff(initp%avg_vapor_gc       ,tiny_offset)
         csite%avg_wshed_vg        (ipa) = sngloff(initp%avg_wshed_vg       ,tiny_offset)
         csite%avg_intercepted     (ipa) = sngloff(initp%avg_intercepted    ,tiny_offset)
         csite%avg_throughfall     (ipa) = sngloff(initp%avg_throughfall    ,tiny_offset)
         csite%avg_vapor_ac        (ipa) = sngloff(initp%avg_vapor_ac       ,tiny_offset)
         csite%avg_transp          (ipa) = sngloff(initp%avg_transp         ,tiny_offset)
         csite%avg_evap            (ipa) = sngloff(initp%avg_evap           ,tiny_offset)
         csite%avg_drainage        (ipa) = sngloff(initp%avg_drainage       ,tiny_offset)
         csite%avg_drainage_heat   (ipa) = sngloff(initp%avg_drainage_heat  ,tiny_offset)
         csite%avg_rshort_gnd      (ipa) = sngloff(initp%avg_rshort_gnd     ,tiny_offset)
         csite%avg_rlong_gnd       (ipa) = sngloff(initp%avg_rlong_gnd      ,tiny_offset)
         csite%avg_sensible_lc     (ipa) = sngloff(initp%avg_sensible_lc    ,tiny_offset)
         csite%avg_sensible_wc     (ipa) = sngloff(initp%avg_sensible_wc    ,tiny_offset)
         csite%avg_qwshed_vg       (ipa) = sngloff(initp%avg_qwshed_vg      ,tiny_offset)
         csite%avg_qintercepted    (ipa) = sngloff(initp%avg_qintercepted   ,tiny_offset)
         csite%avg_qthroughfall    (ipa) = sngloff(initp%avg_qthroughfall   ,tiny_offset)
         csite%avg_sensible_gc     (ipa) = sngloff(initp%avg_sensible_gc    ,tiny_offset)
         csite%avg_sensible_ac     (ipa) = sngloff(initp%avg_sensible_ac    ,tiny_offset)
         csite%avg_carbon_ac       (ipa) = sngloff(initp%avg_carbon_ac      ,tiny_offset)
         do k = rk4site%lsl, nzg
            csite%avg_sensible_gg(k,ipa) = sngloff(initp%avg_sensible_gg(k) ,tiny_offset)
            csite%avg_smoist_gg  (k,ipa) = sngloff(initp%avg_smoist_gg  (k) ,tiny_offset)
            csite%avg_transloss  (k,ipa) = sngloff(initp%avg_transloss  (k) ,tiny_offset)
         end do
         
         !---------------------------------------------------------------------------------!
         !     These variables are integrated here, since they don't change with time.     !
         !---------------------------------------------------------------------------------!
         csite%avg_rlongup        (ipa) = csite%avg_rlongup        (ipa)                   &
                                        + csite%rlongup            (ipa) * dtlsm
         csite%avg_albedo         (ipa) = csite%avg_albedo         (ipa)                   &
                                        + csite%albedo             (ipa) * dtlsm
         csite%avg_albedo_beam    (ipa) = csite%avg_albedo_beam    (ipa)                   &
                                        + csite%albedo_beam        (ipa) * dtlsm
         csite%avg_albedo_diffuse (ipa) = csite%avg_albedo_diffuse (ipa)                   &
                                        + csite%albedo_diffuse     (ipa) * dtlsm
         csite%avg_rlong_albedo   (ipa) = csite%avg_rlong_albedo   (ipa)                   &
                                        + csite%rlong_albedo       (ipa) * dtlsm
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      if(checkbudget) then
         co2budget_loss2atm    = sngloff(initp%co2budget_loss2atm   ,tiny_offset)
         ebudget_loss2atm      = sngloff(initp%ebudget_loss2atm     ,tiny_offset)
         ebudget_loss2drainage = sngloff(initp%ebudget_loss2drainage,tiny_offset)
         ebudget_loss2runoff   = sngloff(initp%ebudget_loss2runoff  ,tiny_offset)
         wbudget_loss2atm      = sngloff(initp%wbudget_loss2atm     ,tiny_offset)
         wbudget_loss2drainage = sngloff(initp%wbudget_loss2drainage,tiny_offset)
         wbudget_loss2runoff   = sngloff(initp%wbudget_loss2runoff  ,tiny_offset)
      else
         co2budget_loss2atm             = 0.
         ebudget_loss2atm               = 0.
         ebudget_loss2drainage          = 0.
         ebudget_loss2runoff            = 0.
         wbudget_loss2atm               = 0.
         wbudget_loss2drainage          = 0.
         wbudget_loss2runoff            = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following is not a pure diagnostic, it is used for phenology and mortality !
      ! functions, preserve this variable and its dependencies in all contexts.            !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! paw_avg - 10-day average of plant available water.                                 !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
!      do ico = 1,cpatch%ncohorts
!         available_water = 0.d0
!         do k = cpatch%krdepth(ico), nzg - 1
!            nsoil = rk4site%ntext_soil(k)
!            available_water = available_water                                              &
!                            + max(0.d0,(initp%soil_water(k) - soil8(nsoil)%soilwp))        &
!                            * (slz8(k+1)-slz8(k))                                          &
!                            / (soil8(nsoil)%slmsts - soil8(nsoil)%soilwp)
!         end do
!         nsoil = rk4site%ntext_soil(nzg)
!         available_water = available_water                                                 &
!                         + max(0.d0,(initp%soil_water(nzg) - soil8(nsoil)%soilwp))         &
!                         * (-1.d0*slz8(nzg))                                               &
!                         / (soil8(nsoil)%slmsts -soil8(nsoil)%soilwp) 
!         available_water = available_water / (-1.d0*slz8(cpatch%krdepth(ico)))
!
!
!         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)            &
!                             + sngl(available_water)*sngl(hdid)/tendays_sec
!      end do

      ! Compute relative porosity, and then weight by thickness of soil layer.
      do k = 1, nzg
         !nsoil = rk4site%ntext_soil(k) !ATT
	 nsoil = csite%ntext_soil(k,ipa)!ATT
         relative_porosity(k) = max(0.d0,(initp%soil_water(k) - soil8(nsoil)%soilwp))         &
              / (soil8(nsoil)%slmsts -soil8(nsoil)%soilwp) 
         if(k < nzg)then
            weighted_relative_porosity(k) = relative_porosity(k) * (slz8(k+1)-slz8(k))
         else
            weighted_relative_porosity(k) = relative_porosity(k) * (-slz8(k))
         endif
      enddo

      ! Cumulative weighted relative porosity, starting from the top layer
      cumul_wei_rel_porosity(nzg) = weighted_relative_porosity(nzg)
      do k = nzg-1,1,-1
         cumul_wei_rel_porosity(k) = cumul_wei_rel_porosity(k+1) + weighted_relative_porosity(k)
      enddo
         
      ! Average relative porosity down to layer k.
      do k = 1, nzg
         avg_rel_porosity(k) = cumul_wei_rel_porosity(k) / (-slz8(k))
      enddo

      ! Update running average of plant available water.
      do k = 1, nzg
         csite%current_paw(k,ipa) = csite%current_paw(k,ipa) + avg_rel_porosity(k) * hdid
      enddo
      
      do k = rk4site%lsl, nzg
         csite%soil_water(k,ipa)   = sngloff(initp%soil_water(k)  ,tiny_offset)
         csite%soil_energy(k,ipa)  = sngloff(initp%soil_energy(k) ,tiny_offset)
         csite%soil_tempk(k,ipa)   = sngloff(initp%soil_tempk(k)  ,tiny_offset)
         csite%soil_fracliq(k,ipa) = sngloff(initp%soil_fracliq(k),tiny_offset)
      end do
      

      !------------------------------------------------------------------------------------!
      !    Surface water energy is computed in J/m² inside the integrator. Convert it back !
      ! to J/kg in the layers that surface water/snow still exists.                        !
      !------------------------------------------------------------------------------------!
      csite%nlev_sfcwater(ipa) = initp%nlev_sfcwater
      do k = 1, csite%nlev_sfcwater(ipa)
         csite%sfcwater_depth(k,ipa)   = sngloff(initp%sfcwater_depth(k)   ,tiny_offset)
         csite%sfcwater_mass(k,ipa)    = sngloff(initp%sfcwater_mass(k)    ,tiny_offset)
         csite%sfcwater_tempk(k,ipa)   = sngloff(initp%sfcwater_tempk(k)   ,tiny_offset)
         csite%sfcwater_fracliq(k,ipa) = sngloff(initp%sfcwater_fracliq(k) ,tiny_offset)
         tmp_energy                    = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
         csite%sfcwater_energy(k,ipa)  = sngloff(tmp_energy                ,tiny_offset)
      end do
      !------------------------------------------------------------------------------------!
      !    For the layers that no longer exist, assign zeroes for prognostic variables,    !
      ! and something for temperature and liquid fraction (just to avoid singularities,    !
      ! and funny numbers in the output, but these values are meaningless and should never !
      ! be used).                                                                          !
      !------------------------------------------------------------------------------------!
      do k = csite%nlev_sfcwater(ipa)+1,nzs
         csite%sfcwater_energy(k,ipa)  = 0.
         csite%sfcwater_mass(k,ipa)    = 0.
         csite%sfcwater_depth(k,ipa)   = 0.
         if (k == 1) then
            csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
            csite%sfcwater_tempk(k,ipa)   = csite%soil_tempk(nzg,ipa)
         else
            csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
            csite%sfcwater_tempk(k,ipa)   = csite%sfcwater_tempk(k-1,ipa)
         end if
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Cohort variables.  Here we must check whether the cohort was really solved or  !
      ! it was skipped after being flagged as "unsafe".  Here the reason why it was flag-  !
      ! ged as such matters.                                                               !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !  LEAVES                                                                         !
         !---------------------------------------------------------------------------------!
         if (initp%leaf_resolvable(ico)) then
            !------------------------------------------------------------------------------!
            !     Leaves were solved, update water and internal energy, and recalculate    !
            ! the temperature and leaf intercellular specific humidity.  The vegetation    !
            ! dry heat capacity is constant within one time step, so it doesn't need to be !
            ! updated.                                                                     !
            !------------------------------------------------------------------------------!
            cpatch%leaf_water(ico)  = sngloff(initp%leaf_water(ico) , tiny_offset)
            cpatch%leaf_energy(ico) = sngloff(initp%leaf_energy(ico), tiny_offset)
            call qwtk(cpatch%leaf_energy(ico),cpatch%leaf_water(ico),cpatch%leaf_hcap(ico) &
                     ,cpatch%leaf_temp(ico),cpatch%leaf_fliq(ico))

            !------------------------------------------------------------------------------!
            !     The intercellular specific humidity is always assumed to be at           !
            ! saturation for a given temperature.  Find the saturation mixing ratio, then  !
            ! convert it to specific humidity.                                             !
            !------------------------------------------------------------------------------!
            cpatch%lint_shv(ico) = rslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
            cpatch%lint_shv(ico) = cpatch%lint_shv(ico) / (1. + cpatch%lint_shv(ico))
            !----- Convert the wind. ------------------------------------------------------!
            cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Copy the conductances.                                                   !
            !------------------------------------------------------------------------------!
            cpatch%leaf_gbh       (ico) = sngloff(initp%leaf_gbh       (ico), tiny_offset)
            cpatch%leaf_gbw       (ico) = sngloff(initp%leaf_gbw       (ico), tiny_offset)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Divide the values of water demand by the time step to obtain the average !
            ! value over the past DTLSM period.                                            !
            !------------------------------------------------------------------------------!
            cpatch%psi_open  (ico) = sngloff(initp%psi_open  (ico),tiny_offset) / hdid
            cpatch%psi_closed(ico) = sngloff(initp%psi_closed(ico),tiny_offset) / hdid

         elseif (cpatch%hite(ico) <=  csite%total_sfcw_depth(ipa)) then
            !------------------------------------------------------------------------------!
            !    For plants buried in snow, fix the leaf temperature to the snow temper-   !
            ! ature of the layer that is the closest to the leaves.                        !
            !------------------------------------------------------------------------------!
            kclosest = 1
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
            end do
            cpatch%leaf_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
            cpatch%leaf_fliq(ico)   = 0.
            cpatch%leaf_water(ico)  = 0.
            cpatch%leaf_energy(ico) = cpatch%leaf_hcap(ico) * cpatch%leaf_temp(ico)
            !------------------------------------------------------------------------------!
            !     The intercellular specific humidity is always assumed to be at           !
            ! saturation for a given temperature.  Find the saturation mixing ratio, then  !
            ! convert it to specific humidity.                                             !
            !------------------------------------------------------------------------------!
            cpatch%lint_shv(ico) = rslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
            cpatch%lint_shv(ico) = cpatch%lint_shv(ico) / (1. + cpatch%lint_shv(ico))
            !----- Copy the meteorological wind to here. ----------------------------------!
            cpatch%veg_wind(ico) = sngloff(rk4site%vels, tiny_offset)
            !----- Make water demand 0. ---------------------------------------------------!
            cpatch%psi_open  (ico) = 0.0
            cpatch%psi_closed(ico) = 0.0

         else
            !------------------------------------------------------------------------------!
            !     For plants with minimal foliage or very sparse patches, fix the leaf     !
            ! temperature to the canopy air space and force leaf_water to be zero.         !
            !------------------------------------------------------------------------------!
            cpatch%leaf_temp(ico)   = csite%can_temp(ipa)
            cpatch%leaf_fliq(ico)   = 0.
            cpatch%leaf_water(ico)  = 0. 
            cpatch%leaf_energy(ico) = cpatch%leaf_hcap(ico) * cpatch%leaf_temp(ico)
            !------------------------------------------------------------------------------!
            !     The intercellular specific humidity is always assumed to be at           !
            ! saturation for a given temperature.  Find the saturation mixing ratio, then  !
            ! convert it to specific humidity.                                             !
            !------------------------------------------------------------------------------!
            cpatch%lint_shv(ico) = rslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
            cpatch%lint_shv(ico) = cpatch%lint_shv(ico) / (1. + cpatch%lint_shv(ico))
            !----- Copy the meteorological wind to here. ----------------------------------!
            cpatch%veg_wind(ico) = sngloff(rk4site%vels, tiny_offset)
            !----- Make water demand 0. ---------------------------------------------------!
            cpatch%psi_open  (ico) = 0.0
            cpatch%psi_closed(ico) = 0.0
         end if
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !  WOOD                                                                           !
         !---------------------------------------------------------------------------------!
         if (initp%wood_resolvable(ico)) then
            !------------------------------------------------------------------------------!
            !     Wood was solved, update water and internal energy, and recalculate       !
            ! the temperature.  The wood dry heat capacity is constant within one time     !
            ! step, so it doesn't need to be updated.                                      !
            !------------------------------------------------------------------------------!
            cpatch%wood_water(ico)  = sngloff(initp%wood_water(ico) , tiny_offset)
            cpatch%wood_energy(ico) = sngloff(initp%wood_energy(ico), tiny_offset)
            call qwtk(cpatch%wood_energy(ico),cpatch%wood_water(ico),cpatch%wood_hcap(ico) &
                     ,cpatch%wood_temp(ico),cpatch%wood_fliq(ico))

            !----- Convert the wind. ------------------------------------------------------!
            cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Copy the conductances.                                                   !
            !------------------------------------------------------------------------------!
            cpatch%wood_gbh       (ico) = sngloff(initp%wood_gbh       (ico), tiny_offset)
            cpatch%wood_gbw       (ico) = sngloff(initp%wood_gbw       (ico), tiny_offset)
            !------------------------------------------------------------------------------!

         elseif (cpatch%hite(ico) <=  csite%total_sfcw_depth(ipa)) then
            !------------------------------------------------------------------------------!
            !    For plants buried in snow, fix the wood temperature to the snow temper-   !
            ! ature of the layer that is the closest to the branches.                      !
            !------------------------------------------------------------------------------!
            kclosest = 1
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
            end do
            cpatch%wood_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
            cpatch%wood_fliq(ico)   = 0.
            cpatch%wood_water(ico)  = 0.
            cpatch%wood_energy(ico) = cpatch%wood_hcap(ico) * cpatch%wood_temp(ico)

            !----- Copy the meteorological wind to here. ----------------------------------!
            cpatch%veg_wind(ico) = sngloff(rk4site%vels, tiny_offset)

         else
            !------------------------------------------------------------------------------!
            !     For very sparse patches of for when wood thermodynamics is off, fix the  !
            ! wood temperature to the canopy air space and force wood_water to be zero.    !
            !------------------------------------------------------------------------------!
            cpatch%wood_temp(ico)   = csite%can_temp(ipa)
            cpatch%wood_fliq(ico)   = 0.
            cpatch%wood_water(ico)  = 0. 
            cpatch%wood_energy(ico) = cpatch%wood_hcap(ico) * cpatch%wood_temp(ico)


            !----- Copy the meteorological wind to here. ----------------------------------!
            if (.not. cpatch%leaf_resolvable(ico)) then
               cpatch%veg_wind(ico) = sngloff(rk4site%vels, tiny_offset)
            end if
         end if

         !---------------------------------------------------------------------------------!
         !     Final sanity check.  This should be removed soon, since it should never     !
         ! happen (well, if this still happens, then it's a bug, and we should remove the  !
         ! bug first...).                                                                  !
         !---------------------------------------------------------------------------------!
         if (cpatch%leaf_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%leaf_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL LEAF_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4site%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',initp%leaf_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',initp%leaf_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',initp%leaf_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',initp%leaf_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',initp%leaf_hcap(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',cpatch%leaf_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',cpatch%leaf_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',cpatch%leaf_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',cpatch%leaf_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',cpatch%leaf_hcap(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         if (cpatch%wood_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%wood_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL WOOD_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4site%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',initp%wood_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',initp%wood_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',initp%wood_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',initp%wood_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',initp%wood_hcap(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',cpatch%wood_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',cpatch%wood_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',cpatch%wood_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',cpatch%wood_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',cpatch%wood_hcap(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do

      !------ Copy the ground variables to the output. ------------------------------------!
      csite%ground_shv (ipa) = sngloff(initp%ground_shv , tiny_offset)
      csite%ground_ssh (ipa) = sngloff(initp%ground_ssh , tiny_offset)
      csite%ground_temp(ipa) = sngloff(initp%ground_temp, tiny_offset)
      csite%ground_fliq(ipa) = sngloff(initp%ground_fliq, tiny_offset)
      return
   end subroutine initp2modelp
   !=======================================================================================!
   !=======================================================================================! 
   subroutine update_water_potential(csite,initp,ipa,lon)  !XXT
   use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use ed_misc_coms    , only : dtlsm                &
								 , current_time 		&
								 , ifoutput! ! intent(in)
      use soil_coms       , only : dslzi8				&
								 , dslz8				&
								 , nzg					&
								 , soil8				&
								 , slz					&
								 , dslz					&
								 , soil! ! intent(in)
      use consts_coms     , only : wdns8				&
								 , wdnsi8				&
								 , cliqvlme8			&
								 , tsupercool8			&
								 , cliq8				&
								 , cice8				&
								 , pi1					
      use rk4_coms        , only : rk4patchtype        
	  use therm_lib8	 , only	 : qwtk8,qtk8
	  use pft_coms		 , only  : water_conductance	&
	  							 , psi50				&
								 , leaf_psi50			&
								 , Cap_leaf				&
								 , Cap_stem				&
								 , SRA					&
								 , root_dens			&
								 , root_beta			&
								 , cuticular_cond
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      type(rk4patchtype)    , target      :: initp
      integer               , intent(in)  :: ipa
	  real					, intent(in)  :: lon

      !----- Local Vars  --------------------------------------------------------------------!
	  type(patchtype)		, pointer		:: cpatch
	  real(kind=8)							:: transp
	  real(kind=8)							:: dsw
	  real(kind=8)		,dimension(nzg)		:: tmp_water_supply
	  real(kind=8)							:: org_soil_tempk 
	  integer								:: ico
	  integer								:: k
	  integer								:: nsoil
	  integer								:: ipft
	  real									:: J_sr
	  real									:: J_rl
	  real									:: ap
	  real									:: bp
	  real									:: org_psi_leaf
	  real									:: org_psi_stem
	  real									:: weighted_soil_psi
	  real									:: weighted_conductance
	  real									:: above_layer_depth
	  real									:: current_layer_depth
	  real									:: RAI
	  real									:: wgpfrac
	  real									:: c_leaf
	  real									:: c_stem
	  real									:: stem_cond
	  real									:: soil_water_cond
      real   , dimension(nzg)               :: layer_psi

	  real   , parameter					:: vessel_length_coef = 1.5 
	  ! modify the distance of water transporation
   	  real						   :: predawn_start
	  real						   :: predawn_end
	  real						   :: midday_start
	  real						   :: midday_end

	  real(kind=8)						   :: rl_energy
	  real(kind=8)						   :: transp_energy ! the internal energy taken away by transpiration

	  cpatch => csite%patch(ipa)
	  ! Calculate layer_psi
   	  do k = 1,nzg
     	nsoil = csite%ntext_soil(k,ipa)
     	wgpfrac = min(1.0,real(initp%soil_water(k) * initp%soil_fracliq(k) / soil8(nsoil)%slmsts))
		if (wgpfrac < 1e-3) then
			layer_psi(k) = -1e6
		else
	     	layer_psi(k) = soil(nsoil)%slpots / wgpfrac ** soil(nsoil)%slbs
		endif
      enddo
	  
	  !Circulate all the cohorts

cohortloop:do ico = 1,cpatch%ncohorts
	  ipft = cpatch%pft(ico)
	  !1. get transpiration
	  transp = (initp%fs_open(ico) * initp%psi_open(ico) 							&
	  		+ (1.0d0 - initp%fs_open(ico)) * initp%psi_closed(ico)) / dtlsm
!	  if(ico ==1) print*,'transp',transp,'psi_open',initp%psi_open(ico)
	  !2. update leaf, stem psi and water flows

		!2.1. Calculate the change of psi_leaf, assuming stem_psi as constant
		if (cpatch%lai(ico) > 0.0) then
			c_leaf = Cap_leaf(ipft) * cpatch%lai(ico) / cpatch%nplant(ico) !* &
!				max(0.01,min(1.0,1 / (1 + (cpatch%psi_leaf(ico) / leaf_psi50(ipft)) ** 10.0)))
			
			stem_cond =  water_conductance(ipft) / 102.   & !convert to kg H2O m-1 s-1 m-1
						 * (1 / (1 + ((cpatch%psi_leaf(ico)+cpatch%psi_stem(ico))/2.0 / psi50(ipft)) ** 6.0)) &  ! stem cavitation
						 * (pi1 * (cpatch%dbh(ico) / 200.) ** 2 &
						 * (cpatch%bsapwood(ico) / (cpatch%bsapwood(ico) + cpatch%bdead(ico)))) & ! sapwood area m2
						 / (cpatch%hite(ico) * vessel_length_coef)														 !height of the tree   m
			! the unit of stem_cond is kgH2O s-1 m-1

			ap = - stem_cond								&
				 / c_leaf		 !kg H2O /m
			! the unit of ap is s-1

			bp = (cpatch%psi_stem(ico) * stem_cond - &	! unit is kgH2O	s-1
				 transp / cpatch%nplant(ico)) & ! unit is kgH2O s-1
				 / c_leaf									 ! unit is  kgH2O /m
			! the unit of bp is m s-1

			! calculate new psi_leaf
			org_psi_leaf = cpatch%psi_leaf(ico)
			cpatch%psi_leaf(ico) = ((ap * org_psi_leaf + bp) *	&!	m s-1
									exp(ap * dtlsm) - bp) / ap

			! calculate the actual sapflow from stem to leaf
			J_rl = (cpatch%psi_leaf(ico) - org_psi_leaf) * c_leaf /	dtlsm + &
							transp / cpatch%nplant(ico)! kgH2O s-1

		else  ! no leaves
			c_leaf = Cap_leaf(ipft) * cpatch%lai(ico) / cpatch%nplant(ico) !* &
!				max(0.01,min(1.0,1 / (1 + (cpatch%psi_leaf(ico) / leaf_psi50(ipft)) ** 10.0)))
			org_psi_leaf = cpatch%psi_leaf(ico)
			cpatch%psi_leaf(ico) = cpatch%psi_stem(ico) - cpatch%hite(ico) ! equal to psi_stem
					
			! calculate the actual sapflow from stem to leaf
			J_rl = (cpatch%psi_leaf(ico) - org_psi_leaf) * c_leaf /	dtlsm + &
							transp / cpatch%nplant(ico)! kgH2O s-1
		endif

			if(isnan(J_rl)) then
				print*,'c_leaf',c_leaf,'org_psi_leaf',org_psi_leaf,'psi_stem',cpatch%psi_stem(ico),'transp',transp,'lai',cpatch%lai(ico),&
					'ap',ap,'bp',bp
			endif
		! 2.2. Update stem_psi and water absorption from the soil assuming
		! soil water potential is constant
		weighted_soil_psi = 0.
		weighted_conductance = 0.
		cpatch%water_supply_layer_frac(:,ico) = 0.
		cpatch%water_supply(ico) = 0.

		do k = cpatch%krdepth(ico),nzg
			current_layer_depth = -slz(k)
			if (k+1 .le. nzg) then
				above_layer_depth = -slz(k+1)
			else
				above_layer_depth = 0.0
			endif
                                       
			!  Calculate RAI
			RAI = cpatch%broot(ico) * & !kgC
        			(root_beta(ipft) ** (above_layer_depth / -slz(cpatch%krdepth(ico))) - &
        			root_beta(ipft) ** (current_layer_depth / -slz(cpatch%krdepth(ico))) ) * &
					SRA(ipft) / &  !	m2/kgC
					(2.0 * cpatch%crown_area(ico) / cpatch%nplant(ico))     ! m2
			! The unit of RAI is m2/m2

			!  Calculate soil water conductance
     		nsoil = csite%ntext_soil(k,ipa)
			wgpfrac = min(1.0,initp%soil_water(k) * initp%soil_fracliq(k) / soil(nsoil)%slmsts)
			! disable hydraulic redistribution
!			if (layer_psi(k) <= cpatch%psi_stem(ico)) then
!				soil_water_cond = 0.0
!			else
			if (cpatch%crown_area(ico) < 1e-6 .or. cpatch%nplant(ico) < 1e-6) then
				soil_water_cond = 0
			else
				soil_water_cond = soil(nsoil)%slcons * wgpfrac ** (2.0 * soil(nsoil)%slbs + 3.0) * 1.e3 &	! kgH2O m-2 s-1
								  * sqrt(RAI) / (pi1 * dslz(k))  &  ! m-1
								  * cpatch%crown_area(ico) / cpatch%nplant(ico)	! m2 
			endif
!			endif
			! The unit of soil_water_cond is kgH2O m-1 s-1

			
			! Calculate weighted conductance, weighted psi, and
			! water_supply_layer_frac
			weighted_conductance = weighted_conductance + soil_water_cond
			weighted_soil_psi = weighted_soil_psi + soil_water_cond * layer_psi(k) ! kgH2O s-1

			if(isnan(weighted_conductance)) then
				print*,'wgpfrac',wgpfrac,'RAI',RAI,'crown_area',cpatch%crown_area(ico),'nplant',cpatch%nplant(ico)
				print*,'broot',cpatch%broot(ico),'root_beta',root_beta(ipft),'above_layer_depth',above_layer_depth,&
				'current_layer_depth',current_layer_depth,'slz_krdepth',-slz(cpatch%krdepth(ico))
				print*,'soil_water',initp%soil_water(20),'soil_temp',initp%soil_tempk(20)
			endif

			cpatch%water_supply_layer_frac(k,ico) = &
					soil_water_cond * (layer_psi(k) - cpatch%psi_stem(ico))!kgH2O s-1 
			cpatch%water_supply(ico) = cpatch%water_supply(ico) + &
						cpatch%water_supply_layer_frac(k,ico)  ! kgH2O s-1
			enddo
			! Update psi_stem
			org_psi_stem = cpatch%psi_stem(ico)
			c_stem = Cap_stem(ipft) * & ! kg H2O m-3 m-1
					(cpatch%broot(ico) + cpatch%bsapwood(ico)) / root_dens(ipft) ! * &! m3
!					max(0.01, (1 / (1 + (cpatch%psi_stem(ico) / psi50(ipft)) ** 6.0)))  ! stem cavitation
			! the unit of c_stem is kgH2O m-1
			! treat c_stem and c_leaf as constant (Zhang et al. 2014; Sack et
					! al. 2003)

			ap = - weighted_conductance  & !kgH2O m-1 s-1
				/ c_stem  ! kgH2O m-1
			!the unit of ap is s-1

			bp = (weighted_soil_psi - J_rl) & ! kgH2O s-1
				/ c_stem					! kgH2O m-1
			!the unit of bp is m s-1
			
			if (ap == 0.) then
				cpatch%psi_stem(ico) = org_psi_stem - J_rl * dtlsm /c_stem
				J_sr = 0.0
			else
				cpatch%psi_stem(ico) = ((ap * org_psi_stem + bp) * exp(ap * dtlsm)- bp) &
									/ ap
				J_sr = (cpatch%psi_stem(ico) - org_psi_stem) * c_stem  / dtlsm + &! kgH2O s-1
						J_rl
			endif

	 	 ! update several water-related vars
	  		cpatch%m_coef(ico) = max(1e-6,min(1.0,1 / (1 + (cpatch%psi_leaf(ico) / leaf_psi50(ipft)) ** 10.0)))
			cpatch%rl_flow(ico) = J_rl
			cpatch%sr_flow(ico) = J_sr
if(isnan(J_sr))print*,'org_psi_stem',org_psi_stem,'ap',ap,'bp',bp,'c_stem',c_stem,'weighted_cond',weighted_conductance
		 ! update averaged and daily values
		 if(ifoutput > 0) then
		 	cpatch%avg_psi_stem(ico) = cpatch%avg_psi_stem(ico) + cpatch%psi_stem(ico)
		 	cpatch%avg_psi_leaf(ico) = cpatch%avg_psi_leaf(ico) + cpatch%psi_leaf(ico)
		 	cpatch%avg_rl_flow(ico) = cpatch%avg_rl_flow(ico) + cpatch%rl_flow(ico)
		 	cpatch%avg_sr_flow(ico) = cpatch%avg_sr_flow(ico) + cpatch%sr_flow(ico)
		 	cpatch%avg_m_coef(ico) = cpatch%avg_m_coef(ico) + cpatch%m_coef(ico)
		 endif

		 
		 ! Integrate Predawn and Midday Values

				! determine the timing of predawn and midday given the lon of
				! the grid cell. Just consider the simplest case here, every
				! 30deg equals to 2hrs....                            XXT

		predawn_start = 5. - nint(lon/30. * 2.)
		predawn_end = 6. - nint(lon/30. * 2)

		midday_start = 13. - nint(lon/30. * 2.)
		midday_end = 14. - nint(lon/30. * 2.)

		if (current_time%hour >= predawn_start .and. current_time%hour < predawn_end ) then
				cpatch%predawn_psi_leaf(ico) = cpatch%predawn_psi_leaf(ico)	+ &
												cpatch%psi_leaf(ico)
		endif

		if (current_time%hour >= midday_start .and. current_time%hour < midday_end ) then
				cpatch%midday_psi_leaf(ico) = cpatch%midday_psi_leaf(ico)	+ &
											cpatch%psi_leaf(ico)
		endif

			! recalculate water_supply_layer_frac and water_supply according
			! to J_sr, the unit should be kgH2O m-2 s-1

			cpatch%water_supply(ico) = J_sr * cpatch%nplant(ico)  !kgH2O m-2 s-1

			if (sum(cpatch%water_supply_layer_frac(:,ico)) /= 0.) then
				cpatch%water_supply_layer_frac(:,ico) = &
					cpatch%water_supply_layer_frac(:,ico) / &
					sum(cpatch%water_supply_layer_frac(:,ico)) &
					* cpatch%water_supply(ico) * dtlsm ! kgH2O m-2
			endif

!				if(ico == 1)print*,'sum',sum(cpatch%water_supply_layer_frac(:,ico)) / cpatch%water_supply(ico) /dtlsm
		tmp_water_supply = dble(cpatch%water_supply_layer_frac(:,ico))

	  	do k = nzg, 1,-1
			nsoil = csite%ntext_soil(k,ipa)
			dsw = max(min(initp%soil_water(k) - soil8(nsoil)%soilcp, &
				tmp_water_supply(k) * wdnsi8 * dslzi8(k)),&
				initp%soil_water(k) - soil8(nsoil)%slmsts)
			if (k - 1 > 0) then
				tmp_water_supply(k-1) = &
					tmp_water_supply(k-1) +  			&
					tmp_water_supply(k) - dsw * &
					dslz8(k) * wdns8
			endif

			tmp_water_supply(k) = dsw * wdns8 * dslz8(k)

		initp%soil_water(k) = initp%soil_water(k) - dsw
		initp%soil_energy(k) = initp%soil_energy(k) - &
								dsw * cliqvlme8 * &
								(initp%soil_tempk(k) - tsupercool8)

		enddo
if (.false.) then
	  print*,'ico',ico
	  print*,'hour',current_time%hour,'psi_stem',cpatch%psi_stem(ico)
	  print*,'psi_leaf',cpatch%psi_leaf(ico)
	  print*,'J_sr',J_sr,'J_rl',J_rl,'T',transp / cpatch%nplant(ico)
	  print*,'tmp_water_supply',tmp_water_supply
	  print*,'m_coef',cpatch%m_coef(ico),'psi_open',cpatch%psi_open(ico)
print*,'water_supply sum norm',sum(real(tmp_water_supply))/cpatch%water_supply(ico) / dtlsm
	  print*,'top soil water',initp%soil_water(18:20)
endif

!----------------------------
! we also need to consider the flow of energy due to rl_flow
!---------------------------
!if (initp%leaf_hcap(ico) /= 0.) then
!	rl_energy = dble(cpatch%rl_flow(ico) * dtlsm * cpatch%nplant(ico)) & !kgH2O/m2
!						* cliq8 * initp%soil_tempk(nzg) ! J/kg/K * K
!						!final unit is J/m2
!
!	transp_energy = dble(transp * cpatch%nplant(ico) * dtlsm) &! kgH2o/m2
!						* cliq8 * initp%leaf_temp(ico) ! J/kg/K * K
!						!final unit is J/m2
!	! update leaf temperature
!    call qwtk8(initp%leaf_energy(ico) + rl_energy - transp_energy, &
!			initp%leaf_water(ico),&
!			initp%leaf_hcap(ico) &
!			+ dble((cpatch%rl_flow(ico) - transp)*dtlsm*cpatch%nplant(ico)) * &
!			cliq8, &
!			initp%leaf_temp(ico),&
!			initp%leaf_fliq(ico))
	! use the updated temperature to back calculate leaf_energy
!	initp%leaf_energy(ico) = initp%leaf_hcap(ico) * initp%leaf_temp(ico)		&
!						   + initp%leaf_water(ico)								&
!						   * (cliq8 * initp%leaf_fliq(ico)						&
!							* (initp%leaf_temp(ico) - tsupercool8)				&
!							+ cice8 * (1.0d0 - initp%leaf_fliq(ico))				&
!							* initp%leaf_temp(ico)				&
!							)
!endif
	  enddo cohortloop
				
	
	  !Finally, update soil temp and fracliq in initp just in case
	  do k = nzg,1,-1
		nsoil = csite%ntext_soil(k,ipa)
		call qwtk8(initp%soil_energy(k),initp%soil_water(k) * wdns8, soil8(nsoil)%slcpd &
				  ,initp%soil_tempk(k),initp%soil_fracliq(k))
	  enddo

!	  print*,'soil_water',initp%soil_water(18:20),'soil_tempk',initp%soil_tempk(18:20)
   end subroutine update_water_potential
end module rk4_driver

!==========================================================================================!
!==========================================================================================!
