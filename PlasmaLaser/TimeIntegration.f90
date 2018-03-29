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

subroutine integrateConvective(u,ux,ua,Nx,Ny,dt,time_s,time_tot)
	!! This function performs the integration of the convective part
	integer i
	! the number of grid points
	integer, intent(in):: Nx,Ny
	! the time step size, the current time and the total simuation time
	real*8,intent(in)::dt,time_s,time_tot
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny):: u,ux,ua
	
	! let's perform the integration 
	!do i=1,Nx*Ny
	!	ua(i)=u(i)-u(i)*(ux(i))*dt
	!end do
	! it is set to zero at boundaries
	do i=1,Ny
		do j=1,Nx
			if(j==Nx .or. i==1 .or. i==Ny) then
				ua((i-1)*Nx+j)=0.
			else if (j==1) then
					if(time_s/time_tot<0.1) then
						ua((i-1)*Nx+j)=u((i-1)*Nx+j)
					else 
						ua((i-1)*Nx+j)=0.
					end if
				 else
					ua((i-1)*Nx+j)=u((i-1)*Nx+j)-u((i-1)*Nx+j)*(ux((i-1)*Nx+j))*dt
			end if		
		end do
	end do
	
end subroutine integrateConvective