subroutine fcn( n, Tp2, dTpdt, iflag )
! ------------------------------------------------------------------
! DERIVES THE SYSTEM OF EQUATIONS
! ------------------------------------------------------------------
	use cpinterface
! 	use omp_lib
	use mod_var, only: dp, degrad, dr_abs, dr_as, dr_cov, dz, HTF, Gas_as, &
						n_z, n_abs, n_as, n_theta, n_cov, nel, &
						Elgeo, Eltyp, Neigh, Tp, Qloss, &
						k_cov, cp_a, cp_c, rho_a, rho_c, &
						r_ci,r_co,r_ai,r_ao,dtheta,P_fi,T_fi, &
						P_as, U_ext, P_ext, T_ext, T_sky, DNI, A_ap, &
						tau_cov, alpha_abs, alpha_cov, V_fi, Dlam, &
						mdot, P_f, U_fi,U_f
	use mod_funcs 
	implicit none
	integer*4 iflag
	integer*4 i,j,k,n
	integer*4 iw,ie,itn,its,ib,it
	integer*4 itype
	real(dp), dimension(nel) :: Tp2, dTpdt
	real(dp), dimension(n_z) :: deltaP
	
	real(dp) Qe,Qw,Qn,Qs,Qt,Qb,Qconv,Qsolar,Qrad,Qbrack,Qf,Qfsum,Qdrop
	real(dp) Awe,Ans,At,Ab,vol
	real(dp) rho,cp,C_f,rho1,rho2
	
	
! 	print *, 'FUNCTION VALUES'

!	pass the values to a new vector to neglect the annular space
	j=1
	do i=1,nel
		if(Eltyp(i)/=2) then	! ignore the annular space
			Tp(i) = Tp2(j)
			j=j+1
		else
		endif
	enddo
	j=1
	
	
	do i=1,n_z
		P_f(i)=P_fi
	enddo

! 	!$OMP DO
! 	do i=1,n_abs*n_theta*n_z
! 		Tp(i)=Tp2(i)
! 	enddo
! 	!$OMP END DO
! 	!$OMP DO
! 	do i=n_theta*n_z*(n_abs+n_as)+1,nel
! 		Tp(i)=Tp2(i-n_theta*n_z*(n_abs+n_as))
! 	enddo
! 	!$OMP END DO
	
	Qfsum = 0_dp
	Qloss = 0_dp
	
	! outermost element loop
	do i=1,nel
		! identify the elements
        iw = Neigh(i,1)	! west
        ie = Neigh(i,2)	! east
        its = Neigh(i,3)	! south
        itn = Neigh(i,4)	! north
        ib = Neigh(i,5)	! bottom
        it = Neigh(i,6)	! top
		! obtain the areas
		Awe = Elgeo(i,1)	! west, east
		Ans = Elgeo(i,2)	! north, south
		At = Elgeo(i,4)	! top
		Ab = Elgeo(i,3)	! bottom
		vol = Elgeo(i,5)	! volume
		! element type
		itype = Eltyp(i)
		! initialiiy set variables to "0"
		Qn=0_dp
		Qs=0_dp
		Qe=0_dp
		Qw=0_dp
		Qt=0_dp
		Qb=0_dp
		Qsolar=0_dp
		Qbrack=0_dp
		Qf=0_dp	! heat from HTF to the absorber
				
		! ------------------------------------------------------------------
		! ABSORBER
		! ------------------------------------------------------------------
		if(itype==1) then	
			rho = rho_a
			cp = cp_a
			
			! recall that east and west faces will never be boundaries since the mesh is closed; thus,
			! there will always be conduction between elements.
			Qw = Qcond(Awe,(r_ai+r_ao)/2*dtheta*degrad,Tp(i),Tp(iw))
			Qe = Qcond(Awe,(r_ai+r_ao)/2*dtheta*degrad,Tp(i),Tp(ie))
			
			! for north and south, the possible interactions are boundary OR
			! conduction to other absorber elements.
			if(itn==-1) then	! boundary
				Qn = 0_dp	! insulated
			else	! other elements
				Qn = Qcond(Ans,dz,Tp(i),Tp(itn))
			endif
			if(its==-1) then	! boundary
				Qs = 0_dp	! insulated
			else	! other elements
				Qs = Qcond(Ans,dz,Tp(i),Tp(its))
			endif
			
			! for top, the element may be below the cover (or annular space) OR
			! touching other absorber elements
			if(Eltyp(it)==2) then	! cover (or annular space)
! 				print *, i, i+n_theta*n_z+n_theta*n_z*n_as, Eltyp(i+n_theta*n_z+n_theta*n_z*n_as)
				Qrad = Qrad_ac(0.d0,Tp(i),Tp(i+n_theta*n_z+n_theta*n_z*n_as))
				Qconv = Qconv_ac(Tp(i),Tp(i+n_theta*n_z+n_theta*n_z*n_as),At,P_as)
				Qsolar = Qsun(tau_cov,1.01_dp,alpha_abs,DNI,At)
				Qbrack = 0 !-Qcond_brack(Tp(i),T_ext,P_ext,U_ext)	! thermal loss through the bracket
				Qt = Qconv + Qrad
				Qloss = Qloss + abs(Qt) + abs(Qbrack)
			else	! adjacent to other elements
				Qt = Qcond(At,dr_abs,Tp(i),Tp(it))
			endif
			! for bottom, the element may be above the HTF OR
			! touching other absorber elements.			
			if(Eltyp(ib)==0) then	! HTF
				Qb = Qconv_f(Tp(i),Tp(ib),P_f(i),U_f(i))
				Qfsum = Qfsum + Qb				
			else	! other elements
				Qb = Qcond(Ab,dr_abs,Tp(i),Tp(ib))
			endif
			
		! ------------------------------------------------------------------
		! ANNULAR SPACE
		! ------------------------------------------------------------------				
		elseif(itype==2) then	
			
		! ------------------------------------------------------------------
		! COVER
		! ------------------------------------------------------------------			
		elseif(itype==3) then	
			rho = rho_c
			cp = cp_c
			
			! recall that east and west faces will never be boundaries since the mesh is closed; thus,
			! there will always be conduction between elements.
			Qw = Qcond(Awe,(r_ci+r_co)/2*dtheta*degrad,Tp(i),Tp(iw))
			Qe = Qcond(Awe,(r_ci+r_co)/2*dtheta*degrad,Tp(i),Tp(ie))
			
			! for north and south, the possible interactions are boundary OR
			! conduction to other cover elements.
			if(itn==-1) then	! boundary
				Qn = 0_dp	! insulated
			else	! other elements
				Qn = Qcond(Ans,dz,Tp(i),Tp(itn))
			endif
			if(its==-1) then	! boundary
				Qs = 0_dp	! insulated
			else	! other elements
				Qs = Qcond(Ans,dz,Tp(i),Tp(its))
			endif
			
			! for top, the element may be at the boundary with the ambient OR
			! touching other cover elements
			if(it==-1) then	! boundary
				Qrad = Qrad_cext(Tp(i),T_sky)
				Qconv = Qconv_cext(U_ext,P_ext,Tp(i),T_ext)
				Qsolar = Qsun(1.0_dp,1.0_dp,alpha_cov,DNI,At)
				Qt = Qconv + Qrad
			else
				Qt = Qcond(At,dr_cov,Tp(i),Tp(it))
			endif

			! for bottom, the element may be above the annular space OR
			! touching other absorber elements.
			if(Eltyp(ib)==2) then	! boundary with the annular space
				Qrad = Qrad_ac(0.d0,Tp(i),Tp(i-(n_theta*n_z+n_theta*n_z*n_as)))
				Qconv = Qconv_ac(Tp(i),Tp(i-(n_theta*n_z+n_theta*n_z*n_as)),Ab,P_as)
				Qb = Qrad + Qconv
			else	! other elements
				Qb = Qcond(Ab,dr_cov,Tp(i),Tp(ib))
			endif
			
	
		! ------------------------------------------------------------------
		! HTF
		! ------------------------------------------------------------------					
		else	
			rho1 = PropsSI("D"//CHAR(0), "T"//CHAR(0), Tp(i), "P"//CHAR(0), P_f(i), HTF)
			rho2 = PropsSI("D"//CHAR(0), "T"//CHAR(0), Tp(i-1), "P"//CHAR(0), P_f(i-1), HTF)!previous density
			rho=abs(rho1-rho2)/2
			cp = PropsSI("C"//CHAR(0), "T"//CHAR(0), Tp(i), "P"//CHAR(0), P_f(i), HTF)

			! the HTF may have inlet, outlet, or other HTF elements as their boundaries,
			! in addition to the absorber. 
			if(n_z==1) then
				U_f(i:n_theta) = (mdot/rho1)/(pi*(Dlam/2)**2)
				U_f(n_theta*n_z*i+1:n_theta*n_z*i+n_theta) = (mdot/rho2)/(pi*(Dlam/2)**2)
				deltaP(i)= 2*C_f*(rho1+rho2)/2*((U_f(i)+U_fi)/2)**2*dz/Dlam !2*C_f*rho*U_f**2*n_z*Dlam
				P_f(i)=P_fi-deltaP(i)
				Qs = enthalpy(T_fi,P_fi,mdot,rho)
				Qn = -enthalpy(Tp(i),P_f(i),mdot,rho)
				
			else	
				if(its==-1) then	! at the inlet
					U_f(i) = (mdot/rho1)/(pi*(Dlam/2)**2)
					U_f(i-1) = (mdot/rho2)/(pi*(Dlam/2)**2)
					deltaP(i)=2*C_f*(rho1+rho2)/2*((U_f(i)+U_fi)/2)**2*dz/Dlam
					P_f(i)=P_fi-deltaP(i)
			
					Qs = enthalpy(T_fi,P_fi,mdot,rho)
					Qn = -enthalpy(Tp(i),P_f(i),mdot,rho)
					
				elseif(itn==-1)	then ! at the outlet
					U_f(i) = (mdot/rho1)/(pi*(Dlam/2)**2)
					U_f(i-1) = (mdot/rho2)/(pi*(Dlam/2)**2)
					deltaP(i)=2*C_f*(rho1+rho2)/2*((U_f(i)+U_f(i-1))/2)**2*dz/Dlam
					P_f(i)=P_f(i-1)-deltaP(i)
					Qs = enthalpy(Tp(i-1),P_f(i-1),mdot,rho)
					Qn = -enthalpy(Tp(i),P_f(i),mdot,rho)
					
				else
					U_f(i) = (mdot/rho1)/(pi*(Dlam/2)**2)
					U_f(i-1) = (mdot/rho2)/(pi*(Dlam/2)**2)
					deltaP(i)=2*C_f*(rho1+rho2)/2*((U_f(i)+U_f(i-1))/2)**2*dz/Dlam
					P_f(i)=P_f(i-1)-deltaP(i)
					Qs = enthalpy(Tp(i-1),P_f(i-1),mdot,rho)
					Qn = -enthalpy(Tp(i),P_f(i),mdot,rho)
					
				endif
			endif
			Qf = -Qfsum/n_z	
		endif	
		
		
		if(itype/=2) then
			dTpdt(j) = 1_dp/( rho*cp*vol )*( Qe+Qw+Qn+Qs+Qt+Qb+Qf+Qsolar+Qbrack )
			j=j+1
	! 		print *, Qe,Qw,Qn,Qs,Qt,Qb,Qsolar,Qf
		else
		endif
	enddo
	
	
	return
	
end subroutine fcn