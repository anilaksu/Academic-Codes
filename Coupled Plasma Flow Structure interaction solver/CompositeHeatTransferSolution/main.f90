
Program main

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                                                         !
   	!                                  Written by Anil A. Aksu                                !
	!     This program is written to solve transient coupled thermoelastic beam equation:     !
	!     At each time step, it first solves 2-D heat equation and generates thermal moment   !
	!    then solves transient Euler beam equation under the effect of the thermal moment     !
	!                                                                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	integer i,j,k,l, jstart,jend
	! the start and the end point of the grid 
	real*8 x_ini, x_end
	! the grid points array and corresponding weight
	real*8, allocatable:: x_grid(:),y_grid(:),w_x(:),w_y(:)
	real*8 arg   	! the dummy argument
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The dimensions of dependent and 			  !
	!             independent variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the length of domain in each direction
	real*8 :: Lx,Ly1, Ly2
	! the number of points in each subgrid in x and y direction
	integer Nx,Ny
	! Boundary condition logical variables to state if it is Dirichlet or Neumann boundary conditions
	logical LeftBC, RightBC, TopBC, BottomBC
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The time integration parameters             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the time integration parameter and the step size
	real*8 theta, dt,time_tot
	! the number of time steps
	integer Ntime
	! the time array 
	real*8, allocatable:: time(:)

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!              Heat Transfer Parameters			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	real*8 kappa1, kappa2, R, rho, cp_al, h, g             ! the heat diffusion coeffients
	real*8, allocatable::uL(:),uR(:),uT(:),uB(:)   		   ! Temperature boundary conditions values 
	real*8, allocatable::Temp_ADI(:,:), Temp(:,:)		   ! Temperature array
	real*8, allocatable::Temp_aux(:,:), Temp_Err(:)		   ! Auxilary Temperature array for ADI time integration and temperature error vector
	real*8, allocatable::M_t(:,:), q(:)			   		   ! Thermal moment and distributed load 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Numerical Operators			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 1D weak diffusion operator 
	real*8, allocatable:: DwDx(:,:), DwDy(:,:) 
	! Boundary condition, Governing Equation and System Matrix and inverse of the system matrix
	real*8, allocatable:: K_right(:,:),K_inv(:,:)
	! ADI matrices in x and y directions
	real*8, allocatable:: Kx1_inv(:,:), Kx1_right(:,:),Kx2_inv(:,:), Kx2_right(:,:)
	real*8, allocatable:: Ky_inv(:,:), Ky_right(:,:)
	!! Differentiation Matrix and the product of the derivatices of lagrange interpolants in x and y direction
	real*8 ,allocatable:: DiffX(:,:),DiffY(:,:)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        2-D Grid and related parameters			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the grid points
	real*8, allocatable:: x_grid2D(:,:),wx(:),wy(:)	
	! Forcing functions
	real*8, allocatable:: F_lbc(:),F_rbc(:),F_bbc(:),F_tbc(:),F_tot(:),F_boundary(:)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                  Useful Constants                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! pi number
	real*8 pi
	! for complex exponential representation
	complex phi
	pi=4.*datan(1.d0)
	phi=CMPLX(0,2.*pi)

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!           2D Transient Heat Equation 	  	          !
	!			   on GLL-GLL Grid     		        	  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! domian properties: number of points in domain and length of domain
	Nx=10
	Ny=5
	! the number of time steps
	Ntime = 100
	! the subdomain length in x direction and y direction
	Lx=100.
	Ly1=20.
	Ly2=10.

	
	rho = 2.8       ! Density of aluminium 
	cp_al = 0.91    ! Specifici heat of aluminium
	kappa1=176.       ! the heat diffusion coefficient in the upper layer
	kappa2=20.        ! the heat diffusion coefficient in the lower layer
	h = 20.;          ! Outside convection coefficient
	!g = 553.;         ! conductance between plates
	g = 200.
	R = 0.5 		  ! the radiation coupling term
	E = 68.           ! Elasticity 
	alp = 2.32*(10.**-5.)! Thermal expansion 
	

	allocate(time(2))                                                            ! Time array 
	
	! Numerical operation variables
	allocate(x_grid(Nx))                                                         ! Lobatto Gauss Legendre point in 1-D in x direction
	allocate(y_grid(2*Ny))                                                       ! Lobatto Gauss Legendre point in 1-D in x direction
	allocate(w_x(Nx))															 ! Corresponding weight in x direction
	allocate(w_y(Ny))															 ! Corresponding weight in y direction
	allocate(DiffX(Nx,Nx))														 ! Differentiation matrix defined on LGL points in x direction
	allocate(DiffY(Ny,Ny))														 ! Differentiation matrix defined on LGL points in y direction
	allocate(DwDx(Nx,Nx))														 ! Weak Diffusion matrix defined on LGL points in x direction
	allocate(DwDy(Ny,Ny))														 ! Weak Diffusion matrix defined on LGL points in y direction
	
	allocate(Kx1_inv(Nx,Nx))											         ! Inverse of the left handside of ADI time integration in x direction for domain 1
	allocate(Kx1_right(Nx,Nx))											         ! The right handside of ADI time integration in x direction for domain 1
	allocate(Kx2_inv(Nx,Nx))											         ! Inverse of the left handside of ADI time integration in x direction for domain 2
	allocate(Kx2_right(Nx,Nx))											         ! The right handside of ADI time integration in x direction for domain 2
	allocate(Ky_inv(2*Ny,2*Ny))											         ! Inverse of the left handside of ADI time integration in y direction
	allocate(Ky_right(2*Ny,2*Ny))											     ! The right handside of ADI time integration in y direction
	
	allocate(K_inv(2*Nx*Ny,2*Nx*Ny))											 ! Inverse of left handside of full heat matrix
	allocate(K_right(2*Nx*Ny,2*Nx*Ny))											 ! Right handside of full heat matrix
	
	allocate(F_lbc(Nx))															 ! Left boundary forcing vector
	allocate(uL(2*Ny))															 ! Left boundary condition
	allocate(F_rbc(Nx))															 ! Right boundary forcing vector
	allocate(uR(2*Ny))															 ! Right boundary condition
	allocate(F_bbc(2*Ny))														 ! Left boundary forcing vector
	allocate(uB(Nx))															 ! Left boundary condition
	allocate(F_tbc(2*Ny))														 ! Right boundary forcing vector
	allocate(uT(Nx))															 ! Right boundary condition
	allocate(Temp_ADI(Ntime,2*Nx*Ny))										     ! Temperature array for ADI
	allocate(Temp(Ntime,2*Nx*Ny))										         ! Temperature array for Full System
	allocate(Temp_aux(2,2*Ny))										         	 ! Auxilarry temperature array for ADI time integration
	allocate(F_tot(2*Nx*Ny))													 ! Total forcing vector
	allocate(F_boundary(2*Nx*Ny))											     ! Boundary forcing vector
		
	call GLLPoints(x_grid,w_x,Nx) 												 ! Grid points and correspondin weights in x direction
	call FirstDiff(DiffX,Nx)													 ! First derivative matrix in x direction
	call GLLPoints(y_grid(1:Ny),w_y,Ny) 												 ! Grid points and correspondin weights in y direction
	call FirstDiff(DiffY,Ny)													 ! First derivative matrix in y direction
	
	x_grid = 0.5*(Lx + Lx*x_grid)                                                ! Mappping form -1 to 1 to 0 to Lx
	DiffX  = (2./Lx)*DiffX														 ! Mappping form -1 to 1 to 0 to Lx
	y_grid(1:Ny) = 0.5*(Ly1+Ly1*y_grid(1:Ny))                                    ! Mappping form -1 to 1 to 0 to Ly1
	y_grid(Ny+1:2*Ny) = 0.001 + Ly1 + (Ly2*y_grid(1:Ny))/Ly1                     ! Mappping form -1 to 1 to Ly1 to Ly1 + Ly2
	DiffY  = (2./Ly1)*DiffY														 ! Mappping form -1 to 1 to 0 to Ly
	 
	
	time(1) = 1500.
	time(2) = 200.
	
	call get1DHeatWeakForm(DwDx, DiffX, w_x, Nx)                                 ! Weak form of the heat equation in x direction
	call get1DHeatWeakForm(DwDy, DiffY, w_y, Ny)								 ! Weak form of the heat equation in x direction
	

	dt =  time(2)/(Ntime-1)        		! time step 
	print*,"Time Step Size", dt
	
	LeftBC  = .TRUE.                   ! It states Neumann boundary condition on the left boundary
	RightBC = .FALSE.				   ! It states Dirichlet boundary condition on the right boundary
	
	TopBC  = .TRUE.                    ! It states Neumann boundary condition on the top boundary
	BottomBC  = .TRUE.                 ! It states Neumann boundary condition on the bottom boundary
	
	uL = 200. 						   ! Heat condition on the left boundary
	uR = 20. 						   ! Heat condition on the right boundary
	uT = 500. 						   ! Heat condition on the top boundary
	uB = 20. 						   ! Heat condition on the bottom boundary
	
	F_lbc = 0.                         ! Left forcing vector
	F_rbc = 0. 						   ! Right forcing vector
	F_tbc = 0.                         ! Top forcing vector
	F_bbc = 0. 						   ! Right forcing vector
	

	
	!Temp_ADI = 20. 						   ! Initial condition 
	Temp = 20.	 						   ! Initial condition 
	
	! Heat matrices in x direction for the first domain 
	! call getHeatSystemXdirection(Kx1_inv, Kx1_right, DwDx, DiffX, LeftBC, RightBC, rho, cp_al, kappa1, Lx, dt, w_x, Nx)  
	! Heat matrices in x direction for the second domain 
	! call getHeatSystemXdirection(Kx2_inv, Kx2_right, DwDx, DiffX, LeftBC, RightBC, rho, cp_al, kappa2, Lx, dt, w_x, Nx)
	! Heat matrix in y directgion
	! call getCompHeatSystemYdirection(Ky_inv,Ky_right,DwDy,DiffY,TopBC,BottomBC,rho,rho,cp_al,cp_al,kappa1,kappa2,g,Ly1,Ly2,dt,w_y,Ny)
	
	! Full heat matrix 
	call getCompHeatSystem(K_inv,K_right,DwDx,DiffX,DwDy,DiffY,rho,rho,cp_al,cp_al,kappa1,kappa2,g,dt,Lx,Ly1,Ly2,w_x,w_y,Nx,Ny)

	
	
	l = 1   ! To DGEMM function
	do i = 1, Ntime -1 
		
		! Full system time integration
		F_tot = 0.
		F_boundary = 0. 
		call DGEMM("N","N",2*Nx*Ny,l,2*Nx*Ny,1.d0,K_right,2*Nx*Ny,Temp(i,:),2*Nx*Ny,0.d0,F_tot,2*Nx*Ny)               ! Forcing due to explicit time part		
		call getTotalCompBoundarForcing(F_boundary,LeftBC, RightBC, TopBC, BottomBC, uL,uR,uT,uB,kappa1,kappa2,Nx,Ny) ! Boundary forcing vector
		F_tot = F_tot + F_boundary
		call DGEMM("N","N",2*Nx*Ny,l,2*Nx*Ny,1.d0,K_inv,2*Nx*Ny,F_tot,2*Nx*Ny,0.d0,Temp(i+1,:),2*Nx*Ny)               ! Implicit part integration	
		
	end do

!	l = 1   ! To DGEMM function
!	do i = 1, Ntime -1 	
!		
!		! Alternating direction time integration
!		
!		! Time integration in y direction
!		do j = 1, Nx
!			! Generate Auxilarry temperature array in y direction
!			do k = 1, 2*Ny
!				Temp_aux(1,k) = Temp_ADI(i,j+(k-1)*Nx)
!			end do
!			! The forcing generated by the previous Temperature distribution in implicit time integration in y direction
!			call DGEMM("N","N",2*Ny,l,2*Ny,1.d0,Ky_right,2*Ny,Temp_aux(1,:),2*Ny,0.d0,F_tbc,2*Ny) 
!			F_tbc(1)  = uB(j)/kappa1     ! Bottom boundary condition
!			F_tbc(2*Ny) = uT(j)/kappa2	 ! Top boundary condition
!			! Temperature distribution in y direction at the next time step
!			call DGEMM("N","N",2*Ny,l,2*Ny,1.d0,Ky_inv,2*Ny,F_tbc,2*Ny,0.d0,Temp_aux(2,:),2*Ny)
!			
!			! Redistribute auxilarry temperature array in y direction to actual temperature array
!			do k = 1, 2*Ny
!				Temp_ADI(i+1,j+(k-1)*Nx) = Temp_aux(2,k) 
!			end do
!			
!		end do
!				
!		! Time integration in x direction
!		do j = 1, 2*Ny
!			jstart = (j-1)*Nx + 1 
!			jend  =   j * Nx
!		
!			if ( j <= Ny) then               ! Heat operator in the first domain
!			
!				! The forcing generated by the previous Temperature distribution in implicit time integration in the first domain
!			    call DGEMM("N","N",Nx,l,Nx,1.d0,Kx1_right,Nx,Temp_ADI(i+1,jstart:jend),Nx,0.d0,F_rbc,Nx) 
!				F_rbc(1)  = uL(j)/kappa1     ! Left boundary condition
!				F_rbc(Nx) = uR(j)			 ! Right boundary condition
!				! Temperature distribution in the first domain at the next time step
!				call DGEMM("N","N",Nx,l,Nx,1.d0,Kx1_inv,Nx,F_rbc,Nx,0.d0,Temp_ADI(i+1,jstart:jend),Nx) 
!				
!			else                             ! Heat operator in the second domain
!				! The forcing generated by the previous Temperature distribution in implicit time integration in the second domain
!			    call DGEMM("N","N",Nx,l,Nx,1.d0,Kx2_right,Nx,Temp_ADI(i+1,jstart:jend),Nx,0.d0,F_rbc,Nx) 
!				F_rbc(1)  = uL(j)/kappa2     ! Left boundary condition
!				F_rbc(Nx) = uR(j)			 ! Right boundary condition
!				! Temperature distribution in the first domain at the next time step
!				call DGEMM("N","N",Nx,l,Nx,1.d0,Kx2_inv,Nx,F_rbc,Nx,0.d0,Temp_ADI(i+1,jstart:jend),Nx) 
!				
!			end if
!			
!		end do
!	end do
	
	open(171,file='Temperature2D.dat',status='unknown')
	open(172,file='FullHeatMatrix.dat',status='unknown')

	
	do j=1, 2*Ny
		do i=1, Nx
			write(171,*) x_grid(i), y_grid(j), Temp(Ntime,i+(j-1)*Nx)
		end do
	end do
	
	do i=1, 2*Nx*Ny
		write(172,*) K_inv(i,:)
		!print*, i, F_boundary(i)
	end do
	
	
	close(171)
	close(172)
		
	deallocate(K_inv)												     ! Resultant system matrices inverse deallaction
	deallocate(K_right)												     ! Resultant system matrices forcing deallaction
	deallocate(F_boundary)												 ! Boundary heating vector
	deallocate(F_tot)													 ! Total forcing vector deallocation
	
	
	allocate(K_right(Nx*Ny,Nx*Ny))										 ! Forcing matrix
	allocate(K_inv(Nx*Ny,Nx*Ny))										 ! Resultant system matrices inverse
	!
	allocate(Temp_Err(Nx*Ny))										     ! Temperature Error vector
	allocate(F_tot(Nx*Ny))												 ! Total forcing vector
	allocate(F_boundary(Nx*Ny))											 ! Boundary heating vector
		
	LeftBC  = .FALSE.                  ! It states Dirichlet boundary condition on the left boundary
	RightBC = .FALSE.				   ! It states Dirichlet boundary condition on the right boundary
	
	TopBC  = .TRUE.                    ! It states Neumann boundary condition on the top boundary
	BottomBC  = .TRUE.                 ! It states Neumann boundary condition on the bottom boundary	
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     Single Domain 2-D Heat Equation Solution        !
	!	  to check accuracy of the order of method        !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	call get2DHeatSystem(K_inv,K_right,DwDx,DiffX,DwDy,DiffY,LeftBC,RightBC,BottomBC,TopBC,rho,cp_al,kappa1,dt,w_x,w_y,Nx,Ny)                ! 2-D Single Domain Heat System for Crank-Nicholson Time Integration

	Temp = 0.                         ! Initial Condition
	
	uL = 200. 						   ! Heat condition on the left boundary
	uR = 2000. 						   ! Heat condition on the right boundary
	
	uT = 0. 						   ! Heat condition on the top boundary
	uB = 0. 						   ! Heat condition on the bottom boundary
	
	! Here, we use analytical solution to get initial temperature distribution
	call getHeat2DAnalyticalSolution(Temp(1,1:Nx*Ny), uL, uR, Lx, Ly1, Nx, Ny, x_grid, y_grid, kappa1/(rho * cp_al), 0. )
	
	l = 1   ! To DGEMM function
	do i = 1, Ntime -1 
		
		! Full system time integration
		F_tot = 0.
		F_boundary = 0. 
		call DGEMM("N","N",Nx*Ny,l,Nx*Ny,1.d0,K_right, Nx*Ny,Temp(i,1:Nx*Ny),Nx*Ny,0.d0,F_tot,Nx*Ny)               ! Forcing due to explicit time part		
		call getTotalBoundarForcing(F_boundary,LeftBC, RightBC, TopBC, BotBC, uL(1:Ny),uR(1:Ny),uT,uB,kappa1,Nx,Ny)   	! Boundary forcing vector
		F_tot = F_tot + F_boundary
		call DGEMM("N","N", Nx*Ny,l, Nx*Ny,1.d0,K_inv, Nx*Ny,F_tot, Nx*Ny,0.d0,Temp(i+1,1:Nx*Ny),Nx*Ny)               ! Implicit part integration	
		
	end do
	
	! Analytical solution
	call getHeat2DAnalyticalSolution(Temp(1,1:Nx*Ny), uL, uR, Lx, Ly1, Nx, Ny, x_grid, y_grid, kappa1/(rho * cp_al), time(2))
	
	open(174,file='SingleHeatMatrix.dat',status='unknown')
	!open(175,file='Temperature2D.dat',status='unknown')
	open(176,file='AbsTempError.dat',status='unknown')
		
	do i=1, Nx*Ny
		write(174,*) K_inv(i,1: Nx * Ny)
		!print*, i, F_boundary(i)
	end do
	
	do j=1, Ny
		do i=1, Nx
			!write(175,*) x_grid(i), y_grid(j), Temp(Ntime,i+(j-1)*Nx)
			Temp_Err(i+(j-1)*Nx) = dabs((Temp(1,i+(j-1)*Nx)-Temp(Ntime-1,i+(j-1)*Nx))/Temp(1,i+(j-1)*Nx))     ! Relative Temperature Error
			!write(176,*) x_grid(i), y_grid(j), Temp_Err(i+(j-1)*Nx)
		end do
	end do
	
	print*, "Maximum Relative Error: ", MAXVAL(Temp_Err)
	write(176,*) Nx, MAXVAL(Temp_Err)
	
	close(174)
	!close(175)
	close(176)
	
End Program main