Module katul_stomatal

Contains

  subroutine katul_lphys(par, ca, Tleaf, gbw, ea, Vm0, lambda0,  &
        air_density, cuticular_cond, leaf_psi,leaf_psi50, &
       accepted_fc, accepted_gsw, A_cl, gsw_cl, leaf_resp, &
	   m_coef,hite_coef,ipft,ico,former_gV,former_gJ)

    use consts_coms, only: mmdry1000,mmh2o,mmdry, &   ! molecular mass of dry air
						   lnexp_min8,lnexp_max8
    use therm_lib, only: eslf
    use physiology_coms, only: resp_inhib
	use farq_leuning, only: arrhenius
	use pft_coms,  only: dark_respiration_factor, &
						 vm_hor,&
						 vm_decay_e,&
						 vm_low_temp,&
						 vm_high_temp
	use ed_misc_coms,  only: current_time

    implicit none

    real, parameter :: co=210.  ! O2 concentration, mmol/mol
    real, parameter :: rd_light_threshold = 10. * 4.6  ! 10 W/m2, converted to uEin.  That is, 46 uEin.

    real, intent(in) :: par ! PAR received by the cohort umol/m2/s
    real, intent(in) :: ca ! atmosphere CO2 mixing ratio umol/mol
    real, intent(in) :: Tleaf ! leaf temperature in degC
    real, intent(in) :: gbw ! Areodynamic resistance for CO2 m2 * s / mol
    real, intent(in) :: ea ! Atmosphere vapor pressure kPa
    real, intent(in) :: Vm0 ! umol/m2/s
    real, intent(in) :: air_density ! kg/m3
    real, intent(in) :: cuticular_cond ! umol/m2/s

    real, intent(in) :: lambda0 ! umol/mol/kPa
	real, intent(in) :: leaf_psi
	real, intent(in) :: leaf_psi50
!    real, intent(in) :: growth_resp ! umol/m2/s
    real, intent(in) :: m_coef  ! correcting coefficient due to leaf water potential
	real, intent(in) :: hite_coef ! correcting coefficient due to vertical position of the cohort

    real, intent(out) :: accepted_gsw  ! final stomatal conductance for water
    real, intent(out) :: A_cl
    real, intent(out) :: gsw_cl
    real, intent(out) :: accepted_fc   ! final carbon fluxes
    real, intent(out) :: leaf_resp     ! leaf respiration
	real, intent(inout) :: former_gV
	real, intent(inout) :: former_gJ
	

	real :: raero
    real :: accepted_g
    real :: ei ! kPa
    real :: Vcmax25 ! umol/m2/s
    real :: Jmax25 ! umol/m2/s
	real :: Vcmax20
	real :: Jmax20
    real :: Jmax  ! umol/m2/s
    real :: Vcmax ! umol/m2/s
    real :: Jrate ! umol/m2/s
    real :: cp  ! umol/mol
    real :: kc  ! umol/mol
    real :: ko  ! mmol/mol
    real :: Rdark ! umol/m2/s
    real :: Rlight ! umol/m2/s
    real :: Rleaf ! umol/m2/s
	real :: cuticular_gsc ! mol/m2/s
	real :: lambda ! umol/mol/kPa
	real :: d2fcdg2,d2fedg2,delta_g
    integer :: ipft
	integer :: ico

    real :: a1gk, a2gk
    real :: k1ci, k2ci, k3ci, k4ci, resid, testg
    real :: myg_V, myg_J
    integer :: iter
    real :: myres, lastresid, testfc1, testfc2, testci, dfcdg
    real :: testderiv, dfedg, myfc_V, myfc_J, myfc

	real :: lnexplow,tlowfun,lnexphigh,thighfun
	real :: t_coef
	real :: Rd0

	integer :: temp_scheme = 3
	integer :: integration_scheme = 2 ! Newton's method

	logical :: is_resolvable
    
    ! Set temperature dependence and other derived variables
    call eslf_katul(Tleaf,ei)
    ei = ei * 0.001 ! kPa
	if (temp_scheme == 1) then
	    Vcmax25 = Vm0 / exp(26.35-65330/(8.314*(273.15+15.)))
	    Vcmax = Vcmax25 * exp(26.35-65330./(8.314*(273.15+Tleaf)))
		Vcmax20 = Vcmax25 * exp(26.35-65330./(8.314*(273.15+20)))
	elseif (temp_scheme == 2) then
	! use original temperature function
		Vcmax = arrhenius(dble(Tleaf+273.15),dble(Vm0),dble(vm_hor(ipft)))
    	lnexplow = vm_decay_e(ipft) * (vm_low_temp(ipft) - (Tleaf))
		lnexplow = min(lnexp_max8,max(lnexp_min8,lnexplow))
		tlowfun = 1 + exp(lnexplow)

		lnexphigh = vm_decay_e(ipft) * ((Tleaf) - vm_high_temp(ipft))
		lnexphigh = min(lnexp_max8,max(lnexp_min8,lnexphigh))
		thighfun = 1 + exp(lnexphigh)
		Vcmax = Vcmax / (tlowfun * thighfun)
		
		Vcmax20 = arrhenius(dble(20+273.15),dble(Vm0),dble(vm_hor(ipft)))
    	lnexplow = vm_decay_e(ipft) * (vm_low_temp(ipft) - (20))
		lnexplow = min(lnexp_max8,max(lnexp_min8,lnexplow))
		tlowfun = 1 + exp(lnexplow)

		lnexphigh = vm_decay_e(ipft) * ((20) - vm_high_temp(ipft))
		lnexphigh = min(lnexp_max8,max(lnexp_min8,lnexphigh))
		thighfun = 1 + exp(lnexphigh)
		Vcmax20 = Vcmax20 / (tlowfun * thighfun)
	elseif (temp_scheme == 3) then
		call harley_temp_fun(Tleaf+273.15,288.15,&
							116.3,  & ! Hv
							0.65,   & ! Sv
							202.9,  & ! Hd
							t_coef)
		Vcmax = Vm0 * t_coef
!		if(ico == 1) print*,'Vm02Vcmax',t_coef

		call harley_temp_fun(20.+273.15,288.15,&
							116.3,  & ! Hv
							0.65,   & ! Sv
							202.9,  & ! Hd
							t_coef)
		Vcmax20 = Vm0 * t_coef
!		if(ico == 1) print*,'Vm02Vcmax20',t_coef


	endif

    Jmax20 = 2.68 * Vcmax20 ! Leuning 1997 Journal of Experimental Botany
    
	call harley_temp_fun(Tleaf+273.15,293.15,&
							79.5,  & ! Hv
							0.65,   & ! Sv
							201.0,  & ! Hd
							t_coef)
	Jmax = Jmax20 * t_coef

    cp = exp(19.02-37830./(8.314*(273.15+Tleaf)))
    kc = exp(38.05-79430./(8.314*(273.15+Tleaf)))
    ko = exp(20.30-36380./(8.314*(273.15+Tleaf)))
    ! Calculate respiration in the dark at 25C.  Respiration rates
    ! per unit leaf area are higher when the leaves are expanding
    ! (phenophase == 1).
!    if(phenophase == 1) then  
	! in TDF version, phenology = 1 indicates the
	! trees are growing leaves    !if(phenophase == 2)then
!       Rdark = 11.73 * exp(-0.5 * green_leaf_factor) - 6.55
!	    Rdark = !11.73 * exp(-0.5 * 0.9) - 6.55
!    else
!       Rdark = !11.73 * exp(-0.5) - 6.55
!    endif

    ! Temperature dependence for respiration in the dark.
	Rd0 = Vm0 * dark_respiration_factor(ipft)

	Rdark = Rd0 * exp(46.39/(8.314e-3*288.15) * (1 - 288.15/(Tleaf+273.15))) 
!	Temp dependence from Harley et al. 1992

    ! Calculate actual leaf respiration depending on whether or not
    ! the leaf is shaded.
	Rleaf = Rdark

!    if(par > rd_light_threshold)then
 !      Rleaf = Rdark * (1.0 - resp_inhib)
  !  endif

	! correcting for leaf water stress and canopy position
	Vcmax = Vcmax * hite_coef * m_coef
	Jmax = Jmax * hite_coef * m_coef
	Rleaf = Rleaf * hite_coef
	cuticular_gsc = cuticular_cond * m_coef * 1.0e-6
	lambda = lambda0 * leaf_psi / leaf_psi50

    ! Solve the quadratic
    Jrate = ((Jmax+0.385*par) -   &
         sqrt((Jmax+0.385*par)**2-4.*0.7*0.385*Jmax*par))/1.4

	if (gbw > 0.) then
		raero = mmdry / gbw * 1.4
	else
		raero = 1e10
	endif

	is_resolvable = (Jmax /= 0.) .and. (Vcmax /= 0.) .and. &
			(raero < 1e8)

	if (is_resolvable) then
	    ! Set up the calculation
    	a1gk = Vcmax
	    a2gk = kc * (1. + co/ko)
    	k1ci = a1gk / ca - Rleaf / ca
	    k2ci = a1gk * raero / ca - Rleaf * raero / ca - 1. + a2gk / ca
    	k3ci = -a1gk*cp/ca/ca-Rleaf*a2gk/ca/ca
	    k4ci = -a1gk*cp/ca*raero/ca-Rleaf*a2gk/ca*raero/ca-a2gk/ca
    	
	    if(ei > ea)then
			select case (integration_scheme)
				case (1)
	 	    	  myres = 0.001
	    		  testg = cuticular_gsc !cuticular_cond * 1.0e-6 / 1.6
    		   	  myg_V = cuticular_gsc !cuticular_cond * 1.0e-6 / 1.6
      		   	  myfc_V = -Rleaf
       
    	   		lastresid = 0.
	   		    do iter = 1, 1000
    		      testg = testg + myres
    		      call fluxsolver(testg, raero, ca, k1ci, k2ci, k3ci, k4ci, testci,testfc1)
    	    	  call fcderiv(testg, raero, k1ci, k2ci, k3ci, k4ci, ca, testci, dfcdg)
    	     	 	dfedg = lambda * 1.6 * (ei - ea) / (1. + 1.6/1.4*testg*raero)**2
    	      		if(lastresid * (dfcdg-dfedg) < 0.)then
    	         		myg_V = testg
		    	         myfc_V = testfc1
    			         exit
		    	      endif
	     		    lastresid = dfcdg - dfedg
		        enddo
		   case (2)

			   ! Use newton's method to find the zero point of 
				  ! start with gsw from last time
				  if ((.not. isnan(former_gV)) .and. former_gV > 1e-10) then
		    		  testg = former_gV
	 			  else
			  	 	  testg = cuticular_gsc
				  endif
		!		  if (Vcmax > 0.) then
				  	  do iter = 1, 500
						! calculate dfcdg - dfedg
						call fluxsolver(testg, raero, ca, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
						call fcderiv(testg, raero, k1ci, k2ci, k3ci, k4ci, ca, testci, dfcdg)
					
			            dfedg = lambda * 1.6 * (ei - ea) / (1. + 1.6/1.4*testg*raero)**2
				   		

		  				call dfcdgderiv(testg, raero, k1ci, k2ci, k3ci, k4ci,ca,testci,dfcdg,d2fcdg2)
							d2fedg2 = - 2. * 1.6 /1.4 * raero / (1.6/1.4 * testg * raero + 1) * &
								dfedg
	
						! calculate the derivative of dfcdg - dfedg
						if (d2fcdg2 - d2fedg2 == 0) then
							delta_g = 0.
						else
							delta_g = - (dfcdg - dfedg) / (d2fcdg2 - d2fedg2)
						endif

						! control exit
						if(abs(dfcdg - dfedg) < 1e-4 .or.                     &			! converge
							(testg < cuticular_gsc .and. dfcdg-dfedg < 0.) .or.   &		! close stomatal
							(testg < cuticular_gsc .and. isnan(dfcdg-dfedg)) .or.   &	! close stomatal
							(testg + delta_g < 0.)					.or.			& 	! unrealistic values
							(delta_g == 0.)					.or.					& 	! trapped
							(testg > 0.4 .and. dfcdg-dfedg > 0.)			  &			! fullyopen stomatal
						) then
							exit
						endif
						testg = testg + delta_g

						!print*,'newton'
						!print*,'iter',iter,'tesg',testg,'delta_g',delta_g
						!print*,'dfcdg - dfedg',dfcdg-dfedg,'d2fcfe',d2fcfe

					enddo
		
!	if(ico == 2)print*,'V testg',testg,'dfcdg-dfedg',dfcdg-dfedg,'iter',iter
!if(ico==45) print*,'V iter',iter,'V testg', testg
			   ! check the case that a negative or no optimal value is found
				   if (testg < cuticular_gsc) then ! .or. &!cuticular_cond * 1.0e-6 / 1.6 .or. &
!					   (iter > 500 .and. dfcdg - dfedg < 0.) .or. &
!					   isnan(testg)) then
	   					testg = cuticular_gsc ! cuticular_cond * 1.0e-6 / 1.6
		  	    	endif
	   		
				   if (testg > 0.4) then 
				   		if (par > 50. .and. dfcdg-dfedg > 0.) then ! light
					   		testg = 0.4
			 			  else ! dark
							testg = cuticular_gsc 
					  	endif
			
			 	   endif

				former_gV = testg
				call fluxsolver(testg, raero, ca, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
				myg_V = testg
				myfc_V = testfc1
			end select
	    else
    	   ! No flux from atmosphere into stomata.  So just set stomatal 
       	! conductance to be really high.
	       myg_V = 1000. / raero
    	   call fluxsolver(myg_V, raero, ca, k1ci, k2ci, k3ci, k4ci, testci,myfc_V)
	    endif

    	! Repeat calculation for light-limited case
	    ! To set this up in the same way as the wc equation, divide by 4.
    	Jrate = Jrate * 0.25
	    a1gk = Jrate
    	a2gk = 2. * cp
	    k1ci = a1gk / ca - Rleaf / ca
    	k2ci = a1gk * raero / ca - Rleaf * raero / ca - 1. + a2gk / ca
	    k3ci = -a1gk*cp/ca/ca-Rleaf*a2gk/ca/ca
    	k4ci = -a1gk*cp/ca*raero/ca-Rleaf*a2gk/ca*raero/ca-a2gk/ca
    
	    if(ei > ea)then
			select case (integration_scheme)
				case (1)
      			 myres = 0.001
	      		 testg = cuticular_gsc !cuticular_cond * 1.0e-6 / 1.6
    	  		 myg_J = cuticular_gsc !cuticular_cond * 1.0e-6 / 1.6
      			 myfc_J = -Rleaf
      			 lastresid = 0.
	       		do iter = 1, 1000
    	    	  testg = testg + myres
        		  call fluxsolver(testg, raero, ca, k1ci, k2ci, k3ci, k4ci, testci,testfc1)
        		  call fcderiv(testg, raero, k1ci, k2ci, k3ci, k4ci, ca, testci, dfcdg)
        	  		dfedg = lambda * 1.6 * (ei - ea) / (1. + 1.6/1.4*testg*raero)**2
        	  		if(lastresid * (dfcdg-dfedg) < 0.)then
        	     		myg_J = testg
        	    	 	myfc_J = testfc1
        	    		exit
        	  		endif
		        	  lastresid = dfcdg - dfedg
       			enddo
				case (2)
	   		! Use newton's method to find the zero point of 
	  		! start with gsw from last time
			!	  if(ico == 1)print*,'Jmax accepted_gsw',accepted_gsw / 1.6 / mmdry

				  if ((.not. isnan(former_gJ)) .and. former_gJ > 1e-10) then
	        			testg = former_gJ
		 		  else
			    	  testg = cuticular_gsc 
				  endif

!			 if(ico == 2) print*,'initial value',testg
				   do iter = 1, 500
					! calculate dfcdg - dfedg
						call fluxsolver(testg, raero, ca, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
						call fcderiv(testg, raero, k1ci, k2ci, k3ci, k4ci, ca, testci, dfcdg)
			
		            	dfedg = lambda * 1.6 * (ei - ea) / (1. + 1.6/1.4*testg*raero)**2
				

  					call dfcdgderiv(testg, raero, k1ci, k2ci, k3ci, k4ci,ca,testci,dfcdg,d2fcdg2)
						d2fedg2 = - 2. * 1.6 /1.4 * raero / (1.6/1.4 * testg * raero + 1) * &
						dfedg
	
					! calculate the derivative of dfcdg - dfedg
						if (d2fcdg2 - d2fedg2 == 0.) then
							delta_g = 0.
						else
							delta_g = - (dfcdg - dfedg) / (d2fcdg2 - d2fedg2)
						endif

						! control exit
						if(abs(dfcdg - dfedg) < 1e-4 .or.                     &			! converge
							(testg < cuticular_gsc .and. dfcdg-dfedg < 0.) .or.   &		! close stomatal
							(testg < cuticular_gsc .and. isnan(dfcdg-dfedg)) .or.   &	! close stomatal
							(testg + delta_g < 0.)					.or.			& 	! unrealistic values
							(delta_g == 0.)					.or.					& 	! trapped
							(testg > 0.4 .and. dfcdg-dfedg > 0.)			  &			! fullyopen stomatal
						) then
							exit
						endif

					testg = testg + delta_g

!if(ico==18)	print*,'iter',iter,'dfcdg-dfedg',dfcdg-dfedg,'d2fcdg2-d2fedg2',d2fcdg2-d2fedg2
				   enddo

!	if(ico == 2)print*,'J testg',testg,'dfcdg-dfedg',dfcdg-dfedg,'iter',iter
			   if (testg < cuticular_gsc ) then!.or. &!cuticular_cond * 1.0e-6 / 1.6 .or. &
		   			testg = cuticular_gsc ! cuticular_cond * 1.0e-6 / 1.6
	  		    endif
	   		
		   		if (testg > 0.4) then 
		   			if (par > 50. .and. dfcdg-dfedg>0.) then ! light
			   			testg = 0.4
		 		  	else ! dark
						testg = cuticular_gsc 
			   		endif
		
	 	   		endif

				former_gJ = testg
				call fluxsolver(testg, raero, ca, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
				myg_J = testg
				myfc_J = testfc1
			end select

    	else
	    	   myg_J = 1000. / raero
    	   	call fluxsolver(myg_J, raero, ca, k1ci, k2ci, k3ci, k4ci, testci,myfc_J)
	    endif


		if(myfc_V < myfc_J)then
    		accepted_fc = myfc_V
	    	accepted_g = myg_V
	    else
    	   accepted_fc = myfc_J
	       accepted_g = myg_J
    	endif
	else  ! not resolvable
		accepted_g = cuticular_gsc
		accepted_fc = -Rleaf

	endif


    A_cl = -Rleaf
    leaf_resp = Rleaf

	accepted_gsw = accepted_g * 1.6 * mmdry
	gsw_cl = cuticular_gsc * 1.6 * mmdry!cuticular_cond * 1.0e-6 * mmdry

!	if(ico == 2)print*,'par',par,'Anet',accepted_fc,'Vcmax',Vcmax,'Jrate',Jrate*4,'m_coef',m_coef,&
!					'Vm0',Vm0,'Rleaf',Rleaf,'VPD',ei-ea,'gsc',accepted_g,'leaf_temp',Tleaf,&
!					'myg_V',myg_V,'myg_J',myg_J,'myfc_V',myfc_V,'myfc_J',myfc_J,'raero',raero,'testci',testci,'lambda',lambda,&
!					'leaf_psi',leaf_psi
  end subroutine katul_lphys

!==========================================================

  
  subroutine fluxsolver(g, ra, ca, k1, k2, k3, k4, ci,fc)
    implicit none
    
    real, intent(in) :: g, k1,k2, k3,k4, ra, ca
    real, intent(out) :: ci, fc
    real :: cip, cim, rad
    
    rad = sqrt((k1/g+k2)**2 - 4. * (k3/g + k4))
    
    cip = ca * (-(k1/g+k2) + rad)/2.
    cim = ca * (-(k1/g+k2) - rad)/2.
    
    ci = cip
    
    fc = (ca - ci) / (1./g + ra)
    
    return
  end subroutine fluxsolver
  
  subroutine fcderiv(g, ra, k1, k2, k3, k4, ca, ci, dfcdg)
    implicit none
    
    real, intent(in) :: g, ra, k1, k2, k3, k4, ca, ci
    real, intent(out) :: dfcdg
    real :: dcidg, myroot
    
    myroot = sqrt((k1/g+k2)**2-4.*(k3/g+k4))
    
    dcidg = ca * (0.5*k1/g**2 +   &
         0.25/myroot*(-2.*k1**2/g**3-2.*k1*k2/g**2+4.*k3/g**2))
    
    dfcdg = ((1./g+ra)*(-dcidg) - (ca-ci)*(-1./g**2))/(1./g+ra)**2
    
    return
  end subroutine fcderiv

  subroutine dfcdgderiv(g, ra, k1, k2, k3, k4,ca,ci,dfcdg,d2fcdg2)
	implicit none
	real, intent(in) :: g,ra,k1,k2,k3,k4,ca,ci,dfcdg
	real, intent(out) :: d2fcdg2
	real :: t1, t2

	t1 = 4. * k3 / g - 2. * k1 * (k2 + k1/g) / g
	t2 = sqrt((k1/g + k2) ** 2 - 4. * (k3/g + k4))

	d2fcdg2 = 1. / ((ra + 1./g) * g ** 2) * ( &
				ca / 2. * (t1 ** 2 / (4. * t2 ** 3) - &
				(k1 ** 2 / g ** 2 - t1) / t2 + 2. * k1 * g) + &
				2. * dfcdg - 2. * (ca - ci) / (g * (ra + 1./g)))

  end subroutine dfcdgderiv
  !======================================================

  subroutine eslf_katul(t, eslf)
    ! This function calculates the saturation vapor pressure over liquid water
    ! as a function of Celsius temperature following Flatau et al. 
    ! (JAM, 1992).

    implicit none
    real, intent(in) :: t
    real             :: x
    real, parameter  :: c0 = .6105851e+03, c1 = .4440316e+02, c2 = .1430341e+01
    real, parameter  :: c3 = .2641412e-01, c4 = .2995057e-03, c5 = .2031998e-05
    real, parameter  :: c6 = .6936113e-08, c7 = .2564861e-11, c8 =-.3704404e-13
    real, intent(out) :: eslf

    x = max(-80.,t)
    eslf = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

    return
  end subroutine eslf_katul
  
  subroutine harley_temp_fun(Tleaf,Tref,Hv,Sv,Hd,T_coef) ! based on Harley et al. 1991
  ! does not consider low temperature cut-off
  implicit none
  real, intent(in) :: Tleaf  ! K
  real, intent(in) :: Tref   ! K
  real, intent(in) :: Hv     ! kJ/mol
  real, intent(in) :: Sv     ! kJ/mol/K
  real, intent(in) :: Hd	 ! kJ/mol
  real, intent(out) :: T_coef  ! unitless

  real,  parameter :: R = 8.314e-3 !kJ/mol/K

		T_coef = exp(Hv/(R * Tref) * (1 - Tref/Tleaf)) / &
				(1 + exp((Sv * Tleaf - Hd)/ (R * Tleaf)))

  end subroutine harley_temp_fun
end Module katul_stomatal
