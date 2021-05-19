 Program main

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                                                         !
   	!                                  Written by Anil A. Aksu                                !
	!  This program is written to solve transient beam equation: It uses Legendre polynomials !
	!  as basis functions and evaluates on Legendre-Gauss-Lobatto points. Here, Euler beam    !
	!  equation is solved under cantilevered conditions at both boundaries which means both   !
	!  funciton and its first derivative are zero at both boundaries. Moreover, Runge Kutta 4 !
	!  is used for numerical integration.                                                     !
	!                                                                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	integer i,j,k,jstart,jend
	! the start and the end point of the grid 
	real*8 x_ini, x_end
	! the grid points array and corresponding weight
	real*8, allocatable:: x_grid(:),y_grid(:),w(:),w_y(:)
	real*8 arg   	! the dummy argument
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The dimensions of dependent and 			  !
	!             independent variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the number of subgrids 
	integer Nsub
	! the length of domain in each direction
	real*8 :: Lx,Ly
	! the number of points in each subgrid in x and y direction
	integer Nx,Ny
	! Number of eigenvalues of the analytical solution
	integer N_eigen
	! the length of subdomain in x direction
	real*8, allocatable::dx(:)
	
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
	!               Dependent Variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! displacement field in x and y direction and theit time derivatives
	real*8, allocatable::w_disp(:,:),wp_disp(:,:)
	real*8, allocatable::Evalues(:), Nfreq(:)				      ! Beam eigenvalues and natural frequencies
	real*8, allocatable::Beta_n(:,:)						      ! Beam eigenrelation coefficients
	real*8, allocatable::M_nm(:,:)						      	  ! Element-wise coefficient matrix
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!              Elastisity Parameters			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	real*8 rho, E, Inertia, Area                               	   ! Density, Elasticity, Inertia and Cross-sectional area
	real*8, allocatable:: q(:)			   		  				   ! Distributed load 

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Numerical Operators			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	real*8, allocatable:: M(:,:),M_inv(:,:)                  	! Mass matrix
	real*8, allocatable:: K_stiff(:,:), K_sinverse(:,:)      	! Beam Stiffness matrix and its inverse
	real*8 ,allocatable:: DiffX(:,:),DiffXX(:,:)  			    ! First and second differentiation Matrices
	real*8, allocatable::Betas(:,:)         	   				! Fourier Mode coefficients 
	real*8, allocatable:: F_beam(:),F_tot(:)   	   				! Beam Forcing functions
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!           2D Transient Heat Equation 	  	          !
	!			   on GLL-GLL Grid     		        	  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the number of subgrids
	Nsub=5 
	! domian properties: number of points in domain and length of domain
	Nx= 15
	Ny= 10

	! the subdomain length in x direction and y direction
	Lx=2.0
	Ly=0.2
	
	! Number of eigenvalues in the analytical solution
	N_eigen = 10
	! the number of time steps
	Ntime = 10000
	
	rho = 2.8       ! Density of aluminium 

	E = 6.9*10.**3.                          ! Elasticity 
	Inertia = (Ly**3.)/12.                   ! Inertia


	allocate(time(2))                                                            ! Time array 
	
	! Numerical operation variables
	allocate(x_grid(Nx))                                                         ! Lobatto Gauss Legendre point in 1-D in x direction
	allocate(y_grid(Ny))                                                         ! Lobatto Gauss Legendre point in 1-D in x direction
	allocate(w(Nx))																 ! Corresponding weight in x direction
	allocate(w_y(Ny))															 ! Corresponding weight in y direction
	allocate(DiffX(Nx,Nx))														 ! Differentiation matrix defined on LGL points in x direction
	
    ! Analytical Solution operators
	allocate(Evalues(N_eigen))													 ! Eigenvalues of the beam equation	
	allocate(Nfreq(N_eigen))													 ! Natural frequencies of the beam equation	
	allocate(Beta_n(N_eigen,2))													 ! Eigenrelation coefficient between basis function of modes
	allocate(M_nm(2,2))															 ! Element-wise coefficient matrix
			
	allocate(F_tot(Nx))											         ! External Forcing of the beam
	allocate(F_beam(Nx-4))												 ! Forcing of the beam due to Crank Nicholson integration scheme
	allocate(q(Nx))										     			 ! Distributed load	
	allocate(w_disp(Ntime,Nx))                                           ! Displacment function 
	allocate(wp_disp(Ntime,Nx))                                           ! Displacment function 
	allocate(Betas(Nx-4,2))                                              ! Beta coefficients for numerical solution of the beam
	
	allocate(M(Nx-4,Nx-4))									     	     ! Resultant Mass matrix
	allocate(M_inv(Nx-4,Nx-4))									     	 ! Inverse of Mass matrix
	allocate(K_stiff(Nx-4,Nx-4))									     ! Resultant stiffness matrix
	allocate(K_sinverse(Nx-4,Nx-4))                                      ! Inverse of he beam stiffness matrix
	allocate(DiffXX(Nx,Nx))                                              ! Second derivative matrix
	
	
	call GLLPoints(x_grid,w,Nx) 												 ! Grid points and correspondin weights in x direction
	call FirstDiff(DiffX,Nx)													 ! First derivative matrix in x direction	
	x_grid = 0.5*(Lx*x_grid)	 	                                             ! Mappping form -1 to 1 to 0 to Lx
	call DGEMM("N","N",Nx,Nx,Nx,1.d0,DiffX,Nx,DiffX,Nx,0.d0,DiffXX,Nx)           ! Second differentiation matrix


	q = 300000.		                         ! distirbuted load
	F_tot = q

	
	!!! Transient solution of the beam
	dt =  10.**-8.       				! time step 
	time(1) = 0.
	time(2) = dt * (Ntime-1) 		    ! Total time
	
	w_disp = 0.                         ! Initial conditions
	wp_disp = 0.

	
	call getBeamBetaCoefficients(Betas, DiffX, Nx)                                      ! Here we get coefficient beta values
	call getBeamStiffnessMatrix(K_stiff, Betas, DiffXX, w, Nx)						 	! Here we compute stiffness matrix
	call getBeamExternalForcing(F_beam, F_tot, Betas, w, Nx)							! Here we compute the corresponding forcing vector
	
	call getBeamMassMatrix(M, w, Betas, Nx) 													! Mass Matrix
	call invertSVD(Nx-4,M,M_inv)																! Inverse of the mass matrix
	call getBeamExternalForcing(F_beam, F_tot, Betas, w, Nx)						        	! Here we compute the corresponding forcing vector
	call DGEMM("N","N",Nx-4,Nx-4,Nx-4,1.d0,M_inv,Nx,K_stiff,Nx,0.d0,K_sinverse,Nx)

	do i=1, Ntime-1	
		call getBeamExternalForcing(F_beam, F_tot, Betas, w, Nx)						        	! Here we compute the corresponding forcing vector
		call integrateRungeKutta4(w_disp(i+1,:),wp_disp(i+1,:),w_disp(i,:),wp_disp(i,:), M_inv, K_stiff,F_beam,Lx,Ly,E,Inertia,rho,dt,Nx)
		!call getBeamNextTime(w_disp(i+1,:), wp_disp(i+1,:),w_disp(i,:),wp_disp(i,:),M_inv, K_sinverse,F_beam,Lx,Ly,E,Inertia,rho,dt,Nx)
		call getBeamAllCoefficients(w_disp(i+1,:), Betas, Nx)
		call getBeamAllCoefficients(wp_disp(i+1,:), Betas, Nx)
	end do
		
	
	!call getBeamEigenvalues(Nfreq,Evalues,N_eigen,Lx, 10.)              		   		! Here we calculate eigenvalues for the Euler beam equation 
	!call getBeamModeCoeffs(Beta_n,Evalues,N_eigen,Lx)					    			! Here we calculate eigenrelation between basis function for each eigenvalue
	!call getCoefficientMatrix(M_nm, Evalues(1),Evalues(2), w, x_grid, Nx,Lx)	   		! Here we generate element-wise coefficient matrix
	!call getBeamModes(M_nm,Evalues, N_eigen, Nx, w, x_grid, Lx)                        ! Let's generate the each mode of the eigen expansion
	
	call getBeamAnalyticalSolution( w_disp(1,:), E*Inertia, rho*Ly, q , time(2), N_eigen, Nx, w, x_grid, Lx)  ! Analytical solution of the beam
	

	!!! Steady-state solution of the beam
	!call invertSVD(Nx-4,K_stiff,K_sinverse)
	!F_beam = -((Lx**4)/(16.*E*Inertia))*F_beam
	!j =1 
	!call DGEMM("N","N",Nx-4,j,Nx-4,1.d0,K_sinverse,Nx-4,F_beam,Nx-4,0.d0,w_disp(2,3:Nx-2),Nx-4)
	
	!call getBeamAllCoefficients(w_disp(2,:), Betas, Nx)                                 ! Let's calculate coefficients at 2 and Nx-1
	
	
	open(176,file='DisplacementAnalytical.dat',status='unknown')
	open(177,file='DisplacementNumerical.dat',status='unknown')
	open(178,file='DisplacementError.dat',status='unknown')
	
	do i=1,  Nx
		write(176,*) x_grid(i), w_disp(1,i)
		write(177,*) x_grid(i), w_disp(Ntime,i)
	end do
	
	! L_inf of the error
	write(178,*) Nx, maxval(dabs(w_disp(1,:)- w_disp(2,:)))

	do i=3,  Nx-2
		print*,x_grid(i), w_disp(Ntime,i), w_disp(1,i)
	end do
	
	close(176)
	close(177)
	close(178)

End Program main