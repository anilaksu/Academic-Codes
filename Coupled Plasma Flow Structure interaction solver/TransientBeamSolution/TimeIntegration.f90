!! this routine contains all required subroutines associated with time integration
!! it has to be used with MatrixOperations.f90 which includes all matrix operations routines used here
subroutine getImplicitIntegrate(M,K,MK,Nx,Ny,theta,dt)
	!! This function generates the matrix required to perform time integration for 
	!! given diffusice matrix and mass matrix
	!! It also requires the time integration coefficient theta
	integer i,j
	! the number of grid points in and y direction
	integer, intent(in):: Nx, Ny
	! theta coefficient and the time step
	real*8, intent(in):: theta, dt
	!the diffusion and the mass matrix and their combination
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: K,M,MK
	! the quantity to be integrated in time and 
	
	MK=M-dt*(1.-theta)*K
	
end subroutine getImplicitIntegrate

subroutine getRHSMatrix(M,K,MK,Nx,Ny,theta,dt)
	!! This function generates the matrix required to perform time integration for 
	!! given diffusive matrix and mass matrix, it results in the explicit integration term
	!! It also requires the time integration coefficient theta
	integer i,j
	! the number of grid points in and y direction
	integer, intent(in):: Nx, Ny
	! theta coefficient and the time step
	real*8, intent(in):: theta, dt
	!the diffusion and the mass matrix and their combination
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: K,M,MK
	! the quantity to be integrated in time and 
	
	MK=M+dt*theta*K
	
end subroutine getRHSMatrix

subroutine integrateAdamBashforth(f_next,f_previous,f_t,dt)
	!! This function integrates function in time
	integer i,j
	real*8, intent(in):: f_previous, f_t, dt  	                ! theta coefficient and the time step
	real*8, intent(inout) :: f_next  	                        ! the diffusion and the mass matrix and their combination
	
    f_next=f_previous+dt*f_t
	
end subroutine integrateAdamBashforth
