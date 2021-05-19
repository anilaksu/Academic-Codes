!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                         !
!            Written by Anil A. Aksu on 04.07.2020 update on 16.05.2021    			      !
!  				Numerical subroutines of Thermoelastic beam solver					      !
!                                                                                         !
!                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getBeamBetaCoefficients(Betas, DiffX, Nx)
	! this function computes the external forcing
	integer i, j, k 
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(inout),dimension(Nx,Nx) :: DiffX	                ! Differentiaon matrix on the Mother interval
	real*8, intent(inout),dimension(Nx-4,2) :: Betas                ! Beta Coefficients
	real*8, dimension(2,2) :: M, M_inv								! Relation matrix between coefficients and its inverse
    real*8, dimension(2) :: V										! Corresponding coefficients
	
	M(1,1) = DiffX(1,2)
	M(1,2) = DiffX(1,Nx-1)
	M(2,1) = DiffX(Nx,2)
	M(2,2) = DiffX(Nx,Nx-1)
	
	call invertSVD(2,M,M_inv)
	
	k = 2 ! to call DGEMM function 
	j = 1 ! to call DGEMM function
	
	
	do i=3,Nx-2
		 V(1) =  -DiffX(1,i)
		 V(2) =  -DiffX(Nx,i)
		 call DGEMM("N","N",k,j,k,1.d0,M_inv,k,V,k,0.d0,Betas(i-2,:),k)
	end do

end subroutine getBeamBetaCoefficients

subroutine getBeamStiffnessMatrix(K_stiff, Betas, DiffXX, w_x, Nx)
	! this function computes the external forcing
	integer i, j, k 
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(in),dimension(Nx) :: w_x			                ! integration weights
	real*8, intent(in),dimension(Nx,Nx) :: DiffXX	                ! Second Differentiaon matrix on the Mother interval
	real*8, intent(in),dimension(Nx-4,2) :: Betas                   ! Beta Coefficients
	real*8, intent(inout),dimension(Nx-4,Nx-4) :: K_stiff           ! Stiffness matrix
	real*8, dimension(Nx,Nx-4) :: RightDiffXX					    ! Right differentiation matrix
	real*8, dimension(Nx-4,Nx) :: LeftDiffXX					    ! Left differentiation matrix

	
	
	do i=3,Nx-2
		RightDiffXX(:,i-2) = DiffXX(:,i) + Betas(i-2,1)*DiffXX(:,2) + Betas(i-2,2)*DiffXX(:,Nx-1)
	end do
	
	call getTranspose(RightDiffXX,LeftDiffXX,Nx,Nx-4)               ! Here we generate Left differentiation matrix by taking the transpose
	
	do i=1,Nx
		!print*, "Right Differentiation Matrix"
		!print*, RightDiffXX(i,:), DiffXX(i,2), DiffXX(i,Nx-1)
		!print*, "Betas"
		!print*, Betas(i,1), Betas(i,2) 
		RightDiffXX(i,:) =  w_x(i) * RightDiffXX(i,:) 
	end do
	
	call matmultnon(LeftDiffXX,RightDiffXX,K_stiff,Nx-4,Nx,Nx-4)
	
end subroutine getBeamStiffnessMatrix

subroutine getBeamExternalForcing(F_beam, F_tot, Betas, w_x, Nx)
	! this function computes the external forcing
	integer i, j, l                                                 ! dummy integer
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(in),dimension(Nx) :: w_x			                ! integration weights
	real*8, intent(in),dimension(Nx-4,2) :: Betas                   ! Beta Coefficients
	real*8, intent(in),dimension(Nx) :: F_tot	                ! Total external forcing
	real*8, intent(inout),dimension(Nx-4) :: F_beam	                ! Total external forcing

	do i=3,Nx-2
		 F_beam(i-2) =  F_tot(i)*w_x(i) + F_tot(2)*Betas(i-2,1)*w_x(2)  + F_tot(Nx-1)*Betas(i-2,2)*w_x(Nx-1)  
	end do
	
	
end subroutine getBeamExternalForcing

subroutine getBeamAllCoefficients(w_disp, Betas, Nx)
	! this function calculates the coefficient 3 and Nx-2 in terms of the rest of the coefficients using Betas
	integer i	                                                   ! dummy integer
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(in),dimension(Nx-4,2) :: Betas                   ! Beta Coefficients
	real*8, intent(inout),dimension(Nx) :: w_disp	                ! Domain coefficients which include from 3 to Nx -2
	
    ! Note that due to boundary condition w(-1)=w(1) =0, coefficient 1 and Nx are already zero 
	do i=3,Nx-2
		 w_disp(2) = w_disp(2) + Betas(i-2,1)*w_disp(i) 
		 w_disp(Nx-1) = w_disp(Nx-1) + Betas(i-2,2)*w_disp(i) 
	end do
	
	
end subroutine getBeamAllCoefficients

subroutine getBeamMassMatrix(M, w_x, Betas, Nx)
	! this function generates the Mass matrix for the beam
	integer i, j	                                                ! dummy integer
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(in),dimension(Nx) :: w_x			                ! integration weights
	real*8, intent(in),dimension(Nx-4,2) :: Betas                   ! Beta Coefficients
	real*8, intent(inout),dimension(Nx-4,Nx-4) :: M			        ! Mass matrix
	
    M = 0.
	
	do i=1,Nx-4
		do j=1,Nx-4
		 	M(i,j) = Betas(i,1)*Betas(j,1)*w_x(2) + Betas(i,2)*Betas(j,2)*w_x(Nx-1)
			if (i == j) then
				M(i,j) = M(i,j)+ w_x(i+2)
			end if
		end do
	end do
	
	
end subroutine getBeamMassMatrix

subroutine integrateRungeKutta4(w_next, wp_next ,w_prev, wp_prev, M_inv, K_stiff, F_beam,  Lx, Ly,  E, Inertia, rho, dt, Nx)
	! this function computes the external forcing
	integer i, j, l                                                 ! dummy integer
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(in) :: rho                       		        ! Density
	real*8, intent(in) :: Inertia                      		        ! Inertia
	real*8, intent(in) :: Lx, Ly                   		        	! Height
	real*8, intent(in) :: E                      		            ! Elasticity
	real*8, intent(in) :: dt                      		            ! time step size
	real*8, intent(in),dimension(Nx-4) :: F_beam		     	    ! Forcing term
	real*8, intent(in),dimension(Nx-4,Nx-4) :: M_inv, K_stiff       ! Inverse mass matrix and stiffness matrix
	real*8, intent(inout),dimension(Nx) :: w_prev, wp_prev          ! Displacement at previous time step
	real*8, intent(inout),dimension(Nx) :: w_next, wp_next          ! Displacement at next time step
	real*8, dimension(Nx-4) :: w_int, wp_int       			        ! Displacement at nintermedite time step
	real*8, dimension(2*(Nx-4),2*(Nx-4)) :: K_mod			        ! Modified stiffness matrix
	real*8, dimension(2*(Nx-4)) :: v_all, k_1, k_2, k_3, k_4, k_d        ! Modified stiffness matrix
	real*8, dimension(2*(Nx-4)) :: F_mod		     	    			! Forcing term multiplied with inverse mass matrix

	
	! Initial Preparation of the integration
	!print*, E, Inertia, rho
	K_mod(1:(Nx-4),(Nx-3):2*(Nx-4)) = -16.*(E*Inertia/(rho*Ly))*M_inv*K_stiff/(Lx**4.)             ! Stiffness matrix multiplied with parameters	
	j=1
	do i=(Nx-3),2*(Nx-4)
		K_mod(i,i-Nx+4) = 1.
	end do
	
	do i=1,2*(Nx-4)
	!	print*, K_mod(i,:)
	end do

	call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,M_inv,Nx-4,F_beam,Nx-4,0.d0,F_mod(1:(Nx-4)),Nx-4)
	F_mod = -F_mod/(rho*Ly)
	
	v_all(1:Nx-4) = wp_prev(3:Nx-2)
	v_all(Nx-3:2*(Nx-4)) = w_prev(3:Nx-2)
	
	!print*, "Runge Kutta Forcing"
	!print*, F_mod
	
	call DGEMM("N","N",2*(Nx-4),j,2*(Nx-4),1.d0,K_mod,2*(Nx-4),v_all,2*(Nx-4),0.d0,k_1,2*(Nx-4))
	
	k_1 = dt*(k_1+F_mod)
	
	!print*, "Runge Kutta Forcing k 1"
	!print*, k_1
	
	k_d = v_all + 0.5* k_1 
	
	call DGEMM("N","N",2*(Nx-4),j,2*(Nx-4),1.d0,K_mod,2*(Nx-4),k_d,2*(Nx-4),0.d0,k_2,2*(Nx-4))
	
	k_2 = dt*(k_2+F_mod)
	
	!print*, "Runge Kutta Forcing k 2"
	!print*, k_2
	
	k_d = v_all + 0.5* k_2
	
	call DGEMM("N","N",2*(Nx-4),j,2*(Nx-4),1.d0,K_mod,2*(Nx-4),k_d,2*(Nx-4),0.d0,k_3,2*(Nx-4))
	
	k_3 = dt*(k_3+F_mod)
	
	!print*, "Runge Kutta Forcing k 3"
	!print*, k_3
	
	k_d = v_all +  k_3
	
	call DGEMM("N","N",2*(Nx-4),j,2*(Nx-4),1.d0,K_mod,2*(Nx-4),k_d,2*(Nx-4),0.d0,k_4,2*(Nx-4))
	
	k_4 = dt*(k_4+F_mod)
	
	!print*, "Runge Kutta Forcing k 4"
	!print*, k_4
	
	v_all = v_all + k_1/6. + k_2/3. + k_3/3. + k_4/6.
	
	wp_next(3:Nx-2) = v_all(1:Nx-4)
	w_next(3:Nx-2) = v_all(Nx-3:2*(Nx-4)) 
	
	
end subroutine integrateRungeKutta4

subroutine getBeamNextTime(w_next, wp_next ,w_prev, wp_prev, M_inv, K_stiff, F_beam,  Lx, Ly,  E, Inertia, rho, dt, Nx)
	! this function computes the external forcing
	integer i, j, l                                                 ! dummy integer
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	real*8, intent(in) :: rho                       		        ! Density
	real*8, intent(in) :: Inertia                      		        ! Inertia
	real*8, intent(in) :: Lx, Ly                   		        	! Height
	real*8, intent(in) :: E                      		            ! Elasticity
	real*8, intent(in) :: dt                      		            ! time step size
	real*8, intent(in),dimension(Nx-4) :: F_beam		     	    ! Forcing term
	real*8, intent(in),dimension(Nx-4,Nx-4) :: M_inv, K_stiff       ! Inverse mass matrix and stiffness matrix
	real*8, intent(inout),dimension(Nx) :: w_prev, wp_prev          ! Displacement at previous time step
	real*8, intent(inout),dimension(Nx) :: w_next, wp_next          ! Displacement at next time step
	real*8, dimension(Nx-4) :: w_int, wp_int       			        ! Displacement at nintermedite time step
	real*8, dimension(Nx-4,Nx-4) :: K_mod					        ! Modified stiffness matrix
	real*8, dimension(Nx-4) :: F_mod		     	    			! Forcing term multiplied with inverse mass matrix
	real*8 w_0, w_1, c_1, c_2, c_3, c_4, d_1, d_2, d_3              ! Yoshida time integration parameters
	
	
	w_0 = -(2.**(1./3.))/(2.-2**(1./3.))
	w_1 = 1./(2.-2.**(1./3.))
	
	c_1 = 0.5 * w_1
	c_2 = 0.5 * (w_0 + w_1)
	c_3 = 0.5 * (w_0 + w_1)
	c_4 = 0.5 * w_1
	
	d_1 = w_1
	d_2 = w_0
	d_3 = w_1
	
	! Initial Preparation of the integration
	K_mod = 16.*(E*Inertia/(rho*Ly))*K_stiff/(Lx**4.)             ! Stiffness matrix multiplied with parameters	
	j=1
	call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,M_inv,Nx-4,F_beam,Nx-4,0.d0,F_mod,Nx-4)
	F_mod = F_mod/(rho*Ly)
	
	! First step of the integration
	w_int = w_prev(3:Nx-2) + c_1*wp_prev(3:Nx-2)*dt
	call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,w_int,Nx-4,0.d0,wp_next(3:Nx-2),Nx-4)
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,wp_prev(3:Nx-2),Nx-4,0.d0,w_next(3:Nx-2),Nx-4) ! Damping term
	wp_int = wp_prev(3:Nx-2)- d_1*dt*(wp_next(3:Nx-2) + F_mod + 0.24*wp_prev(3:Nx-2))
	!wp_int = wp_prev(3:Nx-2)- d_1*dt*(wp_next(3:Nx-2) + F_mod)
	! Second step of the integration
	w_int = w_int  + c_2*wp_int*dt
	call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,w_int,Nx-4,0.d0,wp_next(3:Nx-2),Nx-4)
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,wp_int,Nx-4,0.d0,w_next(3:Nx-2),Nx-4) ! Damping term
	wp_int = wp_int  - d_2*dt*(wp_next(3:Nx-2) + F_mod + 0.24*wp_int)	
	!wp_int = wp_int  - d_2*dt*(wp_next(3:Nx-2) + F_mod)	
	! Third step of the integration
	w_int = w_int  + c_3*wp_int*dt
	call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,w_int,Nx-4,0.d0,wp_next(3:Nx-2),Nx-4)
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,wp_int,Nx-4,0.d0,w_next(3:Nx-2),Nx-4) ! Damping term
	wp_int = wp_int - d_3*dt*(wp_next(3:Nx-2) + F_mod + 0.24*wp_int)
	!wp_int = wp_int - d_3*dt*(wp_next(3:Nx-2) + F_mod)
	
	! Forth step of the integration
	w_int = w_int  + c_4*wp_int*dt
	
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,w_prev(3:Nx-2),Nx-4,0.d0,wp_next(3:Nx-2),Nx-4)
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,wp_prev(3:Nx-2),Nx-4,0.d0,w_next(3:Nx-2),Nx-4)
	
	!wp_int = wp_prev(3:Nx-2) + 0.5*dt*(wp_next(3:Nx-2) - F_mod)
	!w_int = w_prev(3:Nx-2)+ dt* wp_prev(3:Nx-2) + 0.5*(dt**2.)*(wp_next(3:Nx-2) - F_mod)
	
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,w_int,Nx-4,0.d0,wp_next(3:Nx-2),Nx-4)
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_mod,Nx-4,wp_int,Nx-4,0.d0,w_next(3:Nx-2),Nx-4)
	
	!wp_int = wp_int + 0.5*dt*(wp_next(3:Nx-2) - F_mod) 
	
	!-0.01*w_next(3:N-2)
	
	wp_next(3:Nx-2) = wp_int
	w_next(3:Nx-2) = w_int
	!print*, "Time integration output"
	!print*, w_int
	
	
	
end subroutine getBeamNextTime

