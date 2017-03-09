!! this routine does all RHS treatment of the system 
subroutine OneDRHS(RHS,f,N)
	!! This function generates RHS of the solution matrix
	integer i,j,k 
	! the number of grid points
	integer, intent(in):: N
	!the differentiation matrix
	real*8, intent(inout),dimension(N):: RHS,f
	! 1D GLL points and corresponding weight w
	real*8, allocatable::x_1D(:),w(:)
	! 1D GLL points and corresponding weight
	allocate(x_1D(N))
	allocate(w(N))
	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1D,w,N)
	
	do i=1,N
		! at first point is not multiplied with weight
		if(i==1 .or. i==N) then
			RHS(i)=f(i)
		else 
			RHS(i)=f(i)*w(i)
		end if 
	end do
	
end subroutine OneDRHS

subroutine TwoDRHS(RHS,f,Nx,Ny)
	!! This function generates RHS of the solution matrix
	integer i,j,k 
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny):: RHS,f
	! 1D GLL points and corresponding weight w in x and y direction
	real*8, allocatable::x_1Dx(:),wx(:),x_1Dy(:),wy(:)
	! 1D GLL points and corresponding weight in x and y direction
	allocate(x_1Dx(Nx))
	allocate(wx(Nx))
	allocate(x_1Dy(Ny))
	allocate(wy(Ny))
	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dx,wx,Nx)
	call GLLPoints(x_1Dy,wy,Ny)
	
	do i=1,Ny
		do j=1,Nx
			! at first point is not multiplied with weight
			if(i==1 .or. i==Ny .or. j==1 .or. j==Nx) then
				RHS((i-1)*Nx+j)=f((i-1)*Nx+j)
			else 
				RHS((i-1)*Nx+j)=f((i-1)*Nx+j)*wy(i)*wx(j)
			end if 
		end do
	end do
	
end subroutine TwoDRHS