subroutine getHeatEigenvalues(lambda, Eval_x, Eval_y, Lx, Ly, Nx, Ny, alpha)
	! this function calculates the eigenvalues in x and y directions
	integer i, j, l                                                  ! dummy integer
	integer, intent(in):: Nx, Ny    	                             ! the number of points on grid in x and y direction
	real*8, intent(in) :: Lx, Ly                                	 ! Length and width of the domain
	real*8, intent(in) :: alpha	                                	 ! Heat constant = kappa/ rho * cp
	real*8, intent(inout),dimension(Nx) :: Eval_x                    ! Eigenvalues in x direction
	real*8, intent(inout),dimension(Ny) :: Eval_y	                 ! Eigenvalues in y direction
	real*8, intent(inout),dimension(Nx*Ny) :: lambda	             ! Eigenvalues in time
	real*8 pi														 ! Pi number

	pi = 4. * datan(1.d0)
	
	
	! Eigenvalues in x direction
	do i=1,Nx
		Eval_x(i) = 1. * i * pi / Lx 
	end do
	
	! Eigenvalues in y direction
	do i=1,Ny
		!if (i == 1) then
			!Eval_y(i) = 0.
		!else
			Eval_y(i) = 2. * (i) * pi / Ly 
		!end if
	end do
	
	! Eigenvalues in time
	do j=1,Ny
		do i=1, Nx
			lambda(i+(j-1)*Nx) = alpha * (  Eval_x(i)**2. + Eval_y(j)**2.)
		end do
	end do

	
end subroutine  getHeatEigenvalues

subroutine getHeatEigencoefficients(beta, u_L, Lx, Ly, Nx, Ny)
	! this function calculates the eigenvalues in x and y directions
	integer i, j, l                                                  ! dummy integer
	integer, intent(in):: Nx, Ny    	                             ! the number of points on grid in x and y direction
	real*8, intent(in) :: Lx, Ly                                	 ! Length and width of the domain
	real*8, intent(in) :: u_L	                                	 ! Initial temperature distribution constant
	real*8, intent(inout),dimension(Nx*Ny) :: beta	                 ! Eigenvalues in time
	real*8 pi														 ! Pi number

	pi = 4. * datan(1.d0)
	
	! Eigenvalues in time
	do j=1,Ny
		do i=1, Nx
			if ( j == 1) then
				beta(i+(j-1)*Nx) = u_L*4.*(dcos(1. * i * pi )-1.)/((1. * i *pi)**3.)
			else
				!beta(i+(j-1)*Nx) = u_L*12.*((Lx*Ly)**2.)*(dcos(1. * i * pi )-1.)/((1. * i *pi)**3.)
				!beta(i+(j-1)*Nx) = beta(i+(j-1)*Nx) * dcos((j-1) * pi )/ ((pi * (j-1) )** 4.)
				beta(i+(j-1)*Nx) = 0.
			end if 
		end do
	end do

	
end subroutine  getHeatEigencoefficients

subroutine getHeat2DAnalyticalSolution(u, u_0, u_L, Lx, Ly, Nx, Ny, x_grid, y_grid, alpha, t)
	! this function calculates the eigenvalues in x and y directions
	integer i, j, k, l                                               ! dummy integer
	integer, intent(in):: Nx, Ny    	                             ! the number of points on grid in x and y direction
	real*8, intent(in) :: Lx, Ly                                	 ! Length and width of the domain
	real*8, intent(in),dimension(Nx) :: x_grid                     	 ! Grid points in x direction
	real*8, intent(in),dimension(Ny) :: y_grid	                     ! Grid points in y direction
	real*8, intent(in) :: alpha, t	                                 ! Heat constant = kappa/ rho * cp and time 
	real*8, intent(in) :: u_0, u_L	                                 ! Boundary temperatures at Left and Right boundaries
	real*8, intent(inout),dimension(Nx*Ny) :: u	                 	 ! Temperature at time t
	real*8 pi, dummy					 							 ! Pi number
	real*8, dimension(Nx) :: Eval_x                     			 ! Eigenvalues in x direction
	real*8, dimension(Ny) :: Eval_y	                   			     ! Eigenvalues in y direction
	real*8, dimension(Nx*Ny) :: lambda, beta	        	         ! Eigenvalues in time and Eigencoefficients

	pi = 4. * datan(1.d0)
	
	call getHeatEigenvalues(lambda, Eval_x, Eval_y, Lx, Ly, Nx, Ny, alpha)
	
	call getHeatEigencoefficients(beta, u_L, Lx, Ly, Nx, Ny)

	u = 0.  ! Let's reset it
	!print*, "Temperature Distribution"
	do j = 1, Ny           ! Evaluating at grid points in y direction
		do i = 1, Nx       ! Evaluating at grid points in y direction
			do k = 1, Ny   ! Summing eigenmodes in y direction
				do l=1, Nx ! Summing eigenmodes in x direction
					dummy = beta(l+(k-1)*Nx)*dexp(-lambda(l+(k-1)*Nx)*t)*dsin(Eval_x(l)*x_grid(i))*dcos(Eval_y(k)*(y_grid(j)-0.5*Ly))
					u(i+(j-1)*Nx)=u(i+(j-1)*Nx)+dummy
				end do
			end do
			! Non homogeneous part of the solution
			u(i+(j-1)*Nx) = u(i+(j-1)*Nx) + u_0 + (u_L-u_0) * x_grid(i) / Lx
			!u(i+(j-1)*Nx) =  u_0 + (u_L-u_0) * x_grid(i) / Lx
			!print*, u(i+(j-1)*Nx)
		end do
	end do 
	
	
end subroutine  getHeat2DAnalyticalSolution