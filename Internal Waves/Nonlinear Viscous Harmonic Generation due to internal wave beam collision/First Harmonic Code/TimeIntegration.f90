!! this routine contains all required subroutines associated with time integration
subroutine TwoDIDMatrix(Id,Nx,Ny)
	!! This function generates Identity matrix to integrate the system in time
	integer i,j,k 
	! the number of grid points in and y direction
	integer, intent(in):: Nx, Ny
	!the identity matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: Id
	
	Id=0.
	
	do i=2,Ny-1
		do j=2,Nx-1
			! Note that the rows spared for boundary conditions are intentionally left empty
				Id((i-1)*Nx+j,(i-1)*Nx+j)=1.
		end do
	end do
	
end subroutine TwoDIDMatrix
