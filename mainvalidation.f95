program vemPTC
! ------------------------------------------------------------------
!    vemPTR IS A VOLUME ELEMENT MODEL OF A PARABOLIC TROUGH SOLAR RECEIVER
!    (PTR) FORMULATED TO COMPUTE ITS TEMPERATURE GRADIENT WITH RESPECT TO
!    TIME. THE SPATIAL DOMAIN IS DISCRETIZED USING LUMPED CONTROL VOLUMES
!    TO WHICH THE ENERGY BALANCE IS APPLIED AND SOLVED IN TIME.
!
!    CREATORS:     	SAM YANG
!    DATE:          SEPTEMBER 2016
!    AFFILIATION:   FLORIDA STATE UNIVERSITY
!    VERSION:       1.1 (FORTRAN)
!
!    NON-STANDARD LIBRARIES NEEDED:
!       CoolProp - http://www.coolprop.org/ for Fortran
!		LC++ - for c++ interfacing with CoolProp
!		OPENMP - for parallel processing
!	 PLEASE REFER TO Makefile FOR ADDITIONAL LIBRARIES NEEDED BY vemPTC.
!
!    PLEASE REFER TO ReadMe FILE FOR THE ALGORITHM STRUCTURE AND UPDATES
!    THAT HAVE BEEN MADE.
! ------------------------------------------------------------------

	use, intrinsic :: iso_fortran_env
	use cpinterface
	use mod_var, only: Ne,Np,dp,read_input,nel,pi, &
						dTpdt,Tp2,Eltyp,n_as,n_theta,n_z,&
						atol,rtol,Qloss,T_fi,T_sky,U_ext,&
						V_fi,A_ap,DNI,rho,mdot,r_ai,P_fi,HTF,&
						T_ext,U_fi,Dlam
					
	implicit none	
	!include 'nlopt.f'
	character(len=10) time
	character(len=5) user
	integer*4 INFO
	integer*4 i,j,k
	integer*4 nthreads
	real(dp) cpu_ini,cpu_fin
	real(dp), dimension(7) :: Qloss2,DeltaT
	external fcn
	
	real, dimension(7) :: DNI_vec,T_ext_vec,T_fi_vec,DeltaT_loss_vac,Vdot_f_vec,Loss_Vac_Dudley,Loss_Vac_Dudley_Error
	real, dimension(7) :: U_ext_vec
	
	call date_and_time(TIME=time)
	write(*,*) "----------------------------------------------------"
	write(*,*) " vemPTC"
	write(*,*) " CREATOR: SAM YANG"
	write(*,*) " E-MAIL: smyng91@gmail.com"
	write(*,*) " PROGRAM UPDATED: OCTOBER 2016 "
!  	write(*,'(a,2x,a,2x,a)') time
	write(*,*) "----------------------------------------------------"
	
    call cpu_time(cpu_ini)
	
	! READ INPUT DATA
	call read_input
	
	write(*,*) " 1. Generating the mesh..."
	! GENERATE THE MESH AND COORDINATE MATRIX
	call MeshGen
	
	! SEARCH ELEMENT NEIGHBORS
	call NeighSearch		
		
	! EXPORT THE MESH FOR DEBUGGING PURPOSES
	call Visit_VTK
	
	! EXPORT A PROPERTY MATRIX CONTAINING INFORMATION ON:
	!	ELEMENT TYPE, EACH FACE AREA, VOLUME, ARC LENGTH
	call EleProp

		
    write(*,*) "####################################"
    write(*,*) " GENERATED MESH INFO."
    write(*,*) "####################################"
    write(*,'(a11,i10)') " Elements =",nel
    write(*,'(a9,i10)') " Points =",Np
10	write(*,*) "---------------------------------------------------------"
	write(*,*) " Do you want to solve for steady-state solution? (Y/N)"
	write(*,*) "---------------------------------------------------------"
	read *, user
	if((user=='Y').or.(user=='y')) then
		! initialize the steady-state loop
		write(*,*) " 2. Solving the equations for steady-state solution..."

! STEADY-STATE SOLVER

!define vectors to match experimental parameters for thermal loss in vacuum
DNI_vec(1)=889.7 !direct normal insolation (W/m^2)
DNI_vec(2)=889.7 !direct normal insolation (W/m^2)
DNI_vec(3)=889.7 !direct normal insolation (W/m^2)
DNI_vec(4)=889.7 !direct normal insolation (W/m^2)
DNI_vec(5)=872   !direct normal insolation (W/m^2)
DNI_vec(6)=872   !direct normal insolation (W/m^2)
DNI_vec(7)=872   !direct normal insolation (W/m^2)
U_ext_vec(1) = 3.2  !wind speed (m/s)
U_ext_vec(2) = 2.9  !wind speed (m/s)
U_ext_vec(3) = 0.1  !wind speed (m/s)
U_ext_vec(4) = 2    !wind speed (m/s)
U_ext_vec(5) = 1.1  !wind speed (m/s)
U_ext_vec(6) = 1.5  !wind speed (m/s)
U_ext_vec(7) = 0.6  !wind speed (m/s)
T_ext_vec(1) = 26.3 + 273.15 	!ambient temperature (m/s)
T_ext_vec(2) = 25.4 + 273.15 	!ambient temperature (m/s)
T_ext_vec(3) = 22.5 + 273.15 	!ambient temperature (m/s)
T_ext_vec(4) = 26.7 + 273.15 	!ambient temperature (m/s)
T_ext_vec(5) = 19.9 + 273.15 	!ambient temperature (m/s)
T_ext_vec(6) = 24.2 + 273.15 	!ambient temperature (m/s)
T_ext_vec(7) = 27.6 + 273.15 	!ambient temperature (m/s)
T_fi_vec(1) = 99.55 + 273.15 	!inlet temperature (K)
T_fi_vec(2) = 100.02 + 273.15 	!inlet temperature (K)
T_fi_vec(3) = 199.4 + 273.15 	!inlet temperature (K)
T_fi_vec(4) = 299 + 273.15 		!inlet temperature (K)
T_fi_vec(5) = 153.4 + 273.15 	!inlet temperature (K)
T_fi_vec(6) = 253.8 + 273.15 	!inlet temperature (K)
T_fi_vec(7) = 348.3 + 273.15 	!inlet temperature (K)
DeltaT_loss_vac(1) = 74.2		!temp above ambient (K)
DeltaT_loss_vac(2) = 74.6		!temp above ambient (K)
DeltaT_loss_vac(3) = 176.3		!temp above ambient (K)
DeltaT_loss_vac(4) = 271.9		!temp above ambient (K)
DeltaT_loss_vac(5) = 133.1		!temp above ambient (K)
DeltaT_loss_vac(6) = 229.2 		!temp above ambient (K)
DeltaT_loss_vac(7) = 319.9 		!temp above ambient (K)
Vdot_f_vec(1) = 27.4 * 1.66667e-5 !volumetric flow rate (L/min)
Vdot_f_vec(2) = 27.4 * 1.66667e-5 !volumetric flow rate (L/min)
Vdot_f_vec(3) = 54.7 * 1.66667e-5 !volumetric flow rate (L/min)
Vdot_f_vec(4) = 56 * 1.66667e-5   !volumetric flow rate (L/min)
Vdot_f_vec(5) = 53.6 * 1.66667e-5 !volumetric flow rate (L/min)
Vdot_f_vec(6) = 55.6 * 1.66667e-5 !volumetric flow rate (L/min)
Vdot_f_vec(7) = 56.8 * 1.66667e-5 !volumetric flow rate (L/min)
Loss_Vac_Dudley(1) = 0.3		!measured efficiency
Loss_Vac_Dudley(2) = 0.85		!measured efficiency
Loss_Vac_Dudley(3) = 14.04		!measured efficiency
Loss_Vac_Dudley(4) = 36.7		!measured efficiency
Loss_Vac_Dudley(5) = 5.3		!measured efficiency
Loss_Vac_Dudley(6) = 23.4       !measured efficiency
Loss_Vac_Dudley(7) = 55.8       !measured efficiency
Loss_Vac_Dudley_Error(1) = 3.7	!associated error
Loss_Vac_Dudley_Error(2) = 4	!associated error
Loss_Vac_Dudley_Error(3) = 8.5	!associated error
Loss_Vac_Dudley_Error(4) = 8	!associated error
Loss_Vac_Dudley_Error(5) = 7.6	!associated error
Loss_Vac_Dudley_Error(6) = 8.5   !associated error
Loss_Vac_Dudley_Error(7) = 7.3   !associated error


do i = 1,7
    DNI = 0
    T_ext = T_ext_vec(i)
    U_ext = U_ext_vec(i)
    V_fi = Vdot_f_vec(i)
    T_sky = 0.0553*T_ext**(1.5)
    T_fi = T_fi_vec(i)
    rho = PropsSI("D"//CHAR(0), "T"//CHAR(0), T_fi, "P"//CHAR(0), P_fi, HTF)
    mdot = rho*V_fi*10
    U_fi=V_fi/(pi*(Dlam/2)**2) !U_f = mdot/(rho*pi*(r_ai)**2)
    do j=1,nel-(n_as*n_theta*n_z)
		Tp2(j)= T_fi
 	enddo
    call hybrd1 ( fcn, nel-n_z*n_theta*n_as, Tp2, dTpdt, rtol, info )
    Qloss2(i) = Qloss/A_ap   

	if(info==1) then
		write(*,*) " 3. The solution converged successfully!"
	else
		write(*,*) " 3. The solution did not converge for the following reason = ",info
		write(*,*) " 	0, improper input parameters."
		write(*,*) " 	1, algorithm estimates that the relative error between X and the solution is at most TOL.."
		write(*,*) " 	2, number of calls to FCN has reached or exceeded 200*(N+1)."
		write(*,*) " 	3, TOL is too small.  No further improvement in the approximate solution X is possible."
		write(*,*) " 	4, the iteration is not making good progress."		
	endif

!     efficiency(i) =  OPTICAL.ref_cl*OPTICAL.e_ov*1.01*OPTICAL.tau_cover*OPTICAL.alpha_absorber - ...
!         Qloss(i)/(AMBIENT.DNI);
    DeltaT(i) = T_fi - T_ext
        	
enddo

!collecting data for figure(1) Qloss/DeltaT
open (unit=1,file='Qloss2.txt')
do i=1,7
	write (1,*) Qloss2(i),DeltaT(i)
enddo
	endif	! initialize the steady-state loop
    call cpu_time(cpu_fin)
    write(*,*) "####################################"
	print '(" CPU Time = ",f10.3," s.")',cpu_fin-cpu_ini
    write(*,*) "####################################"

end program