module mod_var

	use, intrinsic :: iso_fortran_env
	use cpinterface
	implicit none
	integer*4, parameter :: dp=kind(0.0d0) 
	real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
	real(dp), parameter :: g = 9.807_dp
	real(dp), parameter :: sig = 5.669e-08_dp
	real(dp), parameter :: degrad = pi/180_dp
	
	character*120 trash
	character(len=32) Gas_as, HTF
	
! 	integer*4, ALLOCATABLE :: M(:,:)
	integer*4, ALLOCATABLE :: E(:,:)
	integer*4, ALLOCATABLE :: Neigh(:,:)
	integer*4, ALLOCATABLE :: Eltyp(:)
	integer*4 nel	! number of elements
	integer*4 Ne,Np
	integer*4 n_z, n_theta, n_abs, n_as, n_cov
	
	real(dp), ALLOCATABLE :: P(:,:)
	real(dp), ALLOCATABLE :: cxl(:,:)
	real(dp), ALLOCATABLE :: xl(:,:,:)
	real(dp), ALLOCATABLE :: step_z(:)
	real(dp), ALLOCATABLE :: Tp(:)
	real(dp), ALLOCATABLE :: Tp2(:)
	real(dp), ALLOCATABLE :: Elgeo(:,:)
	real(dp), ALLOCATABLE :: dTpdt(:)
	real(dp), ALLOCATABLE :: P_f(:)
	real(dp), ALLOCATABLE :: U_f(:)
	real(dp) L, r_co, r_ci, r_ao, r_ai, t_cov, t_abs, A_ap
	real(dp) rr,Dlam
	real(dp) dtheta,dz,dr_cov,dr_abs,dr_as
	real(dp) d_co,d_ci,d_ao,d_ai
	real(dp) T_ext, P_ext, U_ext, T_sky, DNI
	real(dp) P_as,gamma_r,d_mol,alpha_as
	real(dp) P_fi, V_fi, T_fi,U_fi
	real(dp) rho,mdot
	real(dp) L_hce, k_brck, D_brck, L_f, Wct
	real(dp) cp_a,cp_c,rho_a,rho_c
	real(dp) epsilon_cov,epsilon_abs,alpha_cov,alpha_abs,tau_cov,k_cov
	real(dp) rim_theta,e_sh,e_tr,e_ge,e_dm,e_un,e_ov,ref_cl
	real(dp) T0,tsim,rtol,atol,sstol
	real(dp) time
	real(dp) Qloss


contains
	
	subroutine read_input
		
		Gas_as = "AIR"//CHAR(0)
		HTF = "INCOMP::S800"//CHAR(0)
		
		open(unit=10,file='input/input.txt')
		
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)L
		read(10,*)d_co
		r_co = d_co/2
		read(10,*)d_ci
		r_ci = d_ci/2
		read(10,*)d_ao
		r_ao = d_ao/2
		read(10,*)d_ai
		r_ai = d_ai/2
		read(10,*)A_ap
		rr = d_ai/d_ao
		Dlam = d_ai*(1+rr**2+(1-rr**2)/log(rr))/(1-rr)**2
					
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)T_ext
		read(10,*)P_ext
		read(10,*)U_ext
		T_sky=0.0553_dp*T_ext**(1.5_dp)
		read(10,*)DNI
		
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)P_as
		read(10,*)gamma_r
		read(10,*)d_mol
		read(10,*)alpha_as
			
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)T_fi
		read(10,*)P_fi
		read(10,*)V_fi
		
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)L_hce
		read(10,*)L_f
		read(10,*)k_brck
		read(10,*)D_brck
		read(10,*)Wct
		
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)epsilon_cov
		read(10,*)epsilon_abs
		read(10,*)alpha_cov
		read(10,*)alpha_abs
		read(10,*)tau_cov
		read(10,*)k_cov
		read(10,*)cp_a
		read(10,*)cp_c
		read(10,*)rho_a
		read(10,*)rho_c
		
		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)rim_theta
		read(10,*)e_sh
		read(10,*)e_tr
		read(10,*)e_ge
		read(10,*)e_dm
		read(10,*)e_un
		e_ov = e_sh*e_tr*e_ge*e_dm*e_un    ! overall error
		read(10,*)ref_cl

		read(10,*)trash
		read(10,*)trash
		read(10,*)trash
		read(10,*)n_z
		read(10,*)n_theta
		read(10,*)n_abs
		read(10,*)n_as	! note that the annular space is ignored for now
		read(10,*)n_cov
		read(10,*)T0
		read(10,*)tsim
		read(10,*)rtol
		read(10,*)atol
		read(10,*)sstol

		! compute the total number of elements and allocate array sizes
		nel=n_z*n_theta*(n_abs+n_as+n_cov)+n_z
		allocate(E((nel),8))
		allocate(Neigh(nel,6))
		allocate(P((nel)*8,3))
		allocate(cxl(nel,3))
		allocate(xl(nel,8,3))
		allocate(step_z(n_z+1))
		allocate(P_f(n_z*n_theta))
		allocate(U_f(nel))
		allocate(Tp(nel))
		allocate(Tp2(nel-n_z*n_theta*n_as))
		allocate(dTpdt(nel-n_z*n_theta*n_as))
		allocate(Elgeo(nel,5))
		allocate(Eltyp(nel))
		
		! compute the element size
		dtheta = 360_dp/n_theta
		dz = L/n_z
		dr_cov = (r_co-r_ci)/n_cov
		dr_abs = (r_ao-r_ai)/n_abs
		dr_as = (r_ci-r_ao)/n_as
		
		close(10)
		
		! compute the initial mean fluid velocity
		U_fi=V_fi/(pi*(Dlam/2)**2)
		rho = PropsSI("D"//CHAR(0), "T"//CHAR(0), T_fi, "P"//CHAR(0), P_fi, HTF)
		mdot = V_fi*rho*10
		
		return
		
	end subroutine read_input

end module mod_var