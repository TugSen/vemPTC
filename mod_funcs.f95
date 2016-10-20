module mod_funcs

	use, intrinsic :: iso_fortran_env
	use cpinterface
	use mod_var, only : dp, pi, sig, g, HTF, Gas_as, n_theta, n_z, dz
	implicit none

contains
	
	! mass transfer
	function enthalpy(Tf,P,mdot,rho)
		real(dp) enthalpy
		real(dp) Tf,P,mdot,rho
		real(dp) cp
		
		cp = PropsSI("C"//CHAR(0), "T"//CHAR(0), Tf, "P"//CHAR(0), P, HTF)
		enthalpy = mdot*(cp*Tf+P/rho)
	
		return
	end function enthalpy
	
	
	! convection from the absorber to the HTF
	function Qconv_f(Ta,Tf,P,U_f)
		use mod_var, only: Dlam,rr
		real(dp) Qconv_f
		real(dp) Nu,rho,k,Pr_f,Pr,Re_D,C_f,mu
		real(dp) Ta,Tf,P,U_f		
	    mu = PropsSI("V"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),Tf,HTF)	! dynamic viscosity of fluid (Pa-s)
	    rho = PropsSI("D"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),Tf,HTF)	! density of fluid (kg/m^3)
	    k = PropsSI("CONDUCTIVITY"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),Tf,HTF)	! thermal conductivity
	    Pr_f = PropsSI("PRANDTL"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),Tf,HTF)	! Prandtl number for the fluid
	    Pr =  PropsSI("PRANDTL"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),(Tf+Ta)/2,HTF)	! Prandtl number average
	    Re_D = U_f*Dlam*rho/mu ! Reynolds number
	    C_f = (1.58*log(Re_D)-3.28)**(-2) ! skin friction coefficient
	    Nu = (1-0.14*rr**(0.6))*C_f/2*(Re_D-1000)*Pr_f*(Pr_f/Pr)**(0.11)/(1+12.7*sqrt(C_f/2)*(Pr_f**(0.6666)-1)) ! Nusselt number 
	    Qconv_f = pi*Nu*k*dz*(Tf-Ta)/n_theta ! heat transfer rate (W)
		
		return
	end function Qconv_f
	
	
	! convection between the absorber and the cover through the annular space
	function Qconv_ac(Ta,Tc,A,P_as)
		use mod_var, only: alpha_as,gamma_r,r_ao,r_ci,d_mol,k_cov
		real(dp) Qconv_ac
		real(dp) A
		real(dp) Ta,Tc,Tf
		real(dp) b,lambda,h_as,P_as
		real(dp) rho,mu,c,k,nu,alpha,Pr,beta
		real(dp) RaD_abs,Nu_as_conv,Nu_as_cond,Nu_as
		
		if(P_as < 0.0133322368_dp) then
			! heat transfer occurs in free molecular regime
			! annular space conduction coefficient
		    b = (2_dp-alpha_as)/alpha_as*((9_dp*gamma_r-5_dp)/(2_dp*(gamma_r+1))) 
			! mean free path (m)
		    lambda = 2.331e-20*(Ta+Tc)/2/(P_as*d_mol**2)
			! compute the heat transfer coefficient
		    h_as = k_cov/(r_ao*log(r_ci/r_ao)+b*lambda*(r_ci/r_ao+1))
			! heat transfer rate per unit length (W/m)
		    Qconv_ac = h_as*pi*r_ao*2*(Tc-Ta)/n_theta*dz
		else
			! heat transfer via natural convection
		    Tf = (Ta+Tc)/2
		    mu = PropsSI("V"//CHAR(0),"T"//CHAR(0),Tf,"P"//CHAR(0),P_as,Gas_as )
		    rho = PropsSI("D"//CHAR(0),"T"//CHAR(0),Tf,"P"//CHAR(0),P_as,Gas_as )
		    c = PropsSI("C"//CHAR(0),"T"//CHAR(0),Tf,"P"//CHAR(0),P_as,Gas_as )
		    k = PropsSI("L"//CHAR(0),"T"//CHAR(0),Tf,"P"//CHAR(0),P_as,Gas_as )
		    nu = mu / rho
		    alpha = k / (c*rho)
		    Pr = nu / alpha
		    beta = 1 / Tf
		    RaD_abs = (g*beta*abs(Ta-Tc)*(r_ao*2)**3)/(alpha*nu)
		    Nu_as_conv = 2/log(1+2/((0.518d0*RaD_abs**(0.25)*(1+(0.559/Pr)**(0.6))**(-0.41667))**15+ &
						(0.1*RaD_abs**(0.33333))**15)**(0.066667)) ! from Kuehn and Goldstein, 1976
		    Nu_as_cond = 2/(cosh(((r_ao*2)**2+(r_ci*2)**2)/(2*2*r_ao*2*r_ci)))
		    Nu_as = (Nu_as_conv**15+Nu_as_cond**15)**(0.066667)
		    Qconv_ac = pi*Nu_as*k*(Tc-Ta)/n_theta*dz
		endif
		return
	end function Qconv_ac
	
	! radiation from the absorber to the cover through the annular space
	function Qrad_ac(alpha,Ta,Tc)
		use mod_var, only: r_ao,r_ci,dtheta,epsilon_abs,epsilon_cov
		real(dp) Qrad_ac
		real(dp) Ta,Tc
		real(dp) alpha
		real(dp) dx1,dx2,dnx1,dnx2
		real(dp) Fac
! 		dx1=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha-(3*dtheta/2)) )
! 		dx2=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha+dtheta/2) )
! 		dnx1=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha-dtheta/2) )
! 		dnx2=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha-dtheta/2) )
! 		Fac = (dx1+dx2-dnx1-dnx2)/(2*r_ao*dtheta)

		Qrad_ac = sig*pi*2*r_ao*dz*(Tc**4-Ta**4)/( 1/epsilon_abs+(1-epsilon_cov)/epsilon_cov*(r_ao/r_ci) )/n_theta
		
		return
	end function Qrad_ac
	
! 	! radiation from the cover to the absorber and the cover through the annular space
! 	function Qrad_ca(alpha)
! 		use mod_var, only: r_ao,r_ci,dtheta
! 		real(dp) Qrad_ca
! 		real(dp) alpha
! 		real(dp) dx1,dx2,dnx1,dnx2
! 		real(dp) Fcc,Fca,Fac,Fcc2
!
! ! 		dx1=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha-(3*dtheta/2)) )
! ! 		dx2=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha+dtheta/2) )
! ! 		dnx1=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha-dtheta/2) )
! ! 		dnx2=sqrt( r_ao**2+r_ci**2-2*r_ao*r_ci*cos(alpha-dtheta/2) )
! ! 		Fac = (dx1+dx2-dnx1-dnx2)/(2*r_ao*dtheta)
! ! 		Fca = r_ao/r_ci*Fac
! ! 		Fcc = 1-sqrt(2*r_ci**2*(1-cos(dtheta)))/(r_ci*dtheta)
! !
! ! 		dx1=sqrt(2*r_ci**2*(1-cos(alpha-dtheta/2)))
! ! 		dx2=sqrt(2*r_ci**2-2*r_ci**2*cos(alpha-dtheta/2))
! ! 		dnx1=sqrt(2*r_ci**2-2*r_ci**2*cos(alpha+dtheta/2))
! ! 		dnx2=sqrt(2*r_ci**2*(1-cos(alpha-3*dtheta/2)))
! ! 		Fcc2 = (dx1+dx2-dnx1-dnx2)/(2*r_ci*dtheta)
!
! 		return
! 	end function Qrad_ca
	
	! conduction between elements
	function Qcond(A,L,T1,T2)
		real(dp) Qcond
		real(dp) k
		real(dp) A,L
		real(dp) T1,T2
		
		k = 0.013_dp*T1 + 15.2_dp
		Qcond = A*k*(T2-T1)/L
	
		return
	end function Qcond
	
	! solar absorption
	function Qsun(tau,cf,alpha,Ib,A)
		use mod_var, only: ref_cl,e_ov,A_ap,r_ao,L
		real(dp) Qsun
		real(dp) tau,cf,alpha
		real(dp) A,Ib
		real(dp) eta_opt
		
		eta_opt = ref_cl*e_ov*cf*tau*alpha
		Qsun = eta_opt*1*1*Ib*A_ap/n_theta/n_z
		
		return
	end function Qsun

	! radiation from the cover to the ambient
	function Qrad_cext(Tc,Tsky)
		use mod_var, only: epsilon_cov, r_co
		real(dp) Qrad_cext
		real(dp) Tc,Tsky
		
	    Qrad_cext =  sig*epsilon_cov*2*r_co*pi*dz*(Tsky**4-Tc**4)/n_theta
		
		return
	end function Qrad_cext
	
	! convection from the cover to the ambient
	function Qconv_cext(Uext,P,Tc,Text)
		use mod_var, only: r_co,k_cov
		real(dp) Qconv_cext
		real(dp) Uext,P,Tc,Text
		real(dp) beta,mu,rho,Pr_c,Pr_ext
		real(dp) Ra_D,Re_D,Nu_w,c,m,n,p2
		real(dp) h_ext
		
	    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),(Tc+Text)/2,"AIR"//CHAR(0))
	    mu = PropsSI("V"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),(Tc+Text)/2,"AIR"//CHAR(0))
	    rho = PropsSI("D"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),(Tc+Text)/2,"AIR"//CHAR(0))
	    Pr_c = PropsSI("PRANDTL"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),Tc,"AIR"//CHAR(0))
	    Pr_ext = PropsSI("PRANDTL"//CHAR(0),"P"//CHAR(0),P,"T"//CHAR(0),Text,"AIR"//CHAR(0))

	    if(Uext <= 0.001) then ! if wind speed ~= 0 m/s
	        Ra_D = g*beta*(Tc-Text)*r_co**3/(mu/rho)**2*Pr_ext
	        Nu_w=(0.6+0.387*(Ra_D/(1+(0.559/Pr_ext)**(9/16))**(16/9))**(1/6))**2
	    else ! if wind speed > 0 m/s
	        Re_D = Uext * r_co*2 / (mu/rho)
	        ! define constants for Zukauskas equation
	        if ((Re_D >= 1_dp) .and. (Re_D < 40_dp)) then
	            c = 0.75
				m = 0.4
	        elseif ((Re_D >= 40_dp) .and. (Re_D < 1000_dp)) then
	            c = 0.51
				m = 0.5
	        elseif ((Re_D >= 1000_dp) .and. (Re_D < 200000_dp)) then
	            c = 0.26
				m = 0.6
	        elseif ((Re_D >= 200000_dp) .and. (Re_D < 1000000_dp)) then
	            c = 0.076
				m = 0.7
			endif
	        if (Pr_ext <= 10_dp) then
	            n = 0.37
	        else
	            n = 0.36
			endif
	        ! Zukauskas Equation
	        p2=0.25	! for fluid heating
	        Nu_w=c*Re_D**m*Pr_ext**n*(Pr_ext/Pr_c)**p2
		endif
	    h_ext = Nu_w*k_cov/(r_co*2)
	    Qconv_cext = h_ext*pi*r_co*2*(Text-Tc)*dz/n_theta
		
		return
	end function Qconv_cext
	
	
	! conduction through the bracket
	function Qcond_brack(T,Text,P,Uext)
		use mod_var, only: D_brck,k_brck,L_f,L_hce,Wct
		real(dp) Qcond_brack
		real(dp) P,Uext
		real(dp) T,Text,Tb,T_avg
		real(dp) beta,mu,rho,c,k,nu,alpha,Pr
		real(dp) Ct,m,C2,G2,C1
		real(dp) Ra_L,Nu_T,Nu_T2,Nu_l,Nu_tot,Re_L,hb
		real(dp) z1,z2,m1,theta_r,theta
		
		Tb = T-10_dp
		beta = 1/((Tb + Text)/2)
		T_avg = (Tb + Text)/3 ! recommended by Forrestall et al.
		mu = PropsSI("V"//CHAR(0),"T"//CHAR(0),T_avg,"P"//CHAR(0),P,"AIR"//CHAR(0)) ! dynamic viscosity
		rho = PropsSI("D"//CHAR(0),"T"//CHAR(0),T_avg,"P"//CHAR(0),P,"AIR"//CHAR(0)) ! density
		c = PropsSI("C"//CHAR(0),"T"//CHAR(0),T_avg,"P"//CHAR(0),P,"AIR"//CHAR(0)) ! specific heat
		k = PropsSI("L"//CHAR(0),"T"//CHAR(0),T_avg,"P"//CHAR(0),P,"AIR"//CHAR(0)) ! thermal conductivity
		nu = mu / rho  ! kinematic viscosity
		alpha = k / (c*rho)	! thermal diffusivity
		Pr = nu / alpha ! Prandtl number
		if(Uext == 0) then
			! natural convection
		   Ct = 0.087 ! constant given in Padilla et al.
		   m = 4.5	! constant given in Padilla et al.
		   C2 = 1.3	! constant given in Padilla et al.
		   G2 = 0.735	! constant given in Padilla et al.
		   Ra_L = (g*beta*(abs(Tb-Text))*D_brck**3)/(alpha*nu) ! Rayleigh number
		   Nu_T = G2*Ct*Ra_L**0.25
		   Nu_T2 = Ct*Ra_L**(0.3333333)
		   Nu_l = C2/(log(1+C2/Nu_T))
		   Nu_tot = (Nu_l**m+Nu_T2**m)**(1/m);
		else
			! forced convection
		    C1 = 0.14 ! constant from Padilla et al.
		    m = 0.666 ! constant from Padilla et al.
		    Re_L = Uext*D_brck/nu ! Reynolds number
		    Nu_tot = C1*Re_L**m	! Neusselt number
		endif
		hb = Nu_tot*k/D_brck ! convective heat transfer coefficient
		theta = T - Text
		z1 = sqrt(hb*k_brck*(2*Wct+2*L_hce)*Wct*L_hce)
		z2 = sqrt(hb*k_brck*4*L_hce*L_hce**2)
		m1 = sqrt(hb*(2*Wct+2*L_hce)/(k_brck*Wct*L_hce))
		theta_r = 1/(cosh(m1*L_hce)+z1/z2*sinh(m1*L_hce)) 
		Qcond_brack = 2*z1*theta/L_f*(cosh(m1*L_hce)-theta_r)/sinh(m1*L_hce)*dz
		return
	end function Qcond_brack

end module mod_funcs