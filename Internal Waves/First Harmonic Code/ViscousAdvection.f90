!! this routine generates 1-D and 2-D Viscous Dissipation Operators with Spectral Methods
!! on GLL-GLL grid and finite differences
subroutine ViscousAdvectionSpectral(VisAdvMatrix,Nx,Ny,Lx,Ly,Cgx,Cgy,alpha)
	!! This function returns Viscous Advection operator on 2D GLL-GLL grid
	integer i,j,k 
	!! the start and end of indices
	integer istart,iend,jstart,jend
	! the number of grid points
	integer, intent(in):: Nx,Ny
	! the domain length in x and y direction
	real*8, intent(in):: Lx,Ly
	! the group velocities and the viscous dissipation coefficients
	real*8, intent(in):: Cgx,Cgy,alpha
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: VisAdvMatrix
	! 1D GLL points and corresponding weight w in x and y direction
	real*8, allocatable::x_1Dx(:),wx(:),x_1Dy(:),wy(:)
	!! Differentiation Matrix and the product of the derivatices of lagrange interpolants in x and y direction
	real*8,allocatable:: D1x(:,:),ldpnx(:,:), D1y(:,:),ldpny(:,:)
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dx(Nx))
	allocate(wx(Nx))
	allocate(x_1Dy(Ny))
	allocate(wy(Ny))
	! the allocation of first derivative matrix and the product of the derivatices of lagrange interpolants
	allocate(D1x(Nx,Nx))
	allocate(ldpnx(Nx,Nx))
	allocate(D1y(Ny,Ny))
	allocate(ldpny(Ny,Ny))
	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dx,wx,Nx)
	call GLLPoints(x_1Dy,wy,Ny)
	
	! let's first generate differentiation matrices
	call FirstDiff(D1x,Nx)
	call FirstDiff(D1y,Ny)
	
	! the derivative in x direction multiplied with group velocity in x direction
	do i=1,Ny
		!do j=1,N
		istart=Nx*(i-1)
		iend=Nx*i
		! x-derivative 
		do j=1,Nx
			do k=1,Nx
				VisAdvMatrix(istart+j,istart+k)=2.*Cgx*D1x(j,k)/Lx
				!*wx(j)*wy(i)*Ly
				!print*,LaplaceMatrix(istart+j,istart+k)
			end do 
		end do
	end do
	
	! the derivative in y direction multiplied with group velocity in y direction
	do i=1,Nx
		! y-derivative of laplacian
		do j=1,Ny
			do k=1,Ny
				VisAdvMatrix(i+Nx*(j-1),i+Nx*(k-1))=VisAdvMatrix(i+Nx*(j-1),i+Nx*(k-1))+2.*Cgy*D1y(j,k)/Ly
				!*wy(j)*wx(i)*Lx
				!print*,LaplaceMatrix(istart+j,istart+k)
			end do 
		end do
	end do
	
	! the viscous dissipation term
	do i=1,Ny		
		do j=1,Nx
			! the counter 
			istart=(i-1)*Nx+j
			
			VisAdvMatrix(istart,istart)=VisAdvMatrix(istart,istart)+alpha
			!*wx(j)*wy(i)*Lx*Ly/4.
		end do
	end do
	
end subroutine ViscousAdvectionSpectral


subroutine ViscousAdvectionFinite(VisAdvMatrix,Nx,Ny,Lx,Ly,Cgx,Cgy,alpha)
	!! This function returns Viscous Advection operator on uniform 2D grid 
	integer i,j,k 
	!! the start and end of indices
	integer istart,iend,jstart,jend
	! the number of grid points
	integer, intent(in):: Nx,Ny
	! the domain length in x and y direction
	real*8, intent(in):: Lx,Ly
	! the group velocities and the viscous dissipation coefficients
	real*8, intent(in):: Cgx,Cgy,alpha
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: VisAdvMatrix
	! the finite step sizes
	real*8 :: dx,dy
	
	! the step sizes 
	dx=Lx/(Nx-1)
	dy=Ly/(Ny-1)
	
	! the derivative in x direction multiplied with group velocity in x direction
	do i=1,Ny
		istart=Nx*(i-1)
		! x-derivative 
		do j=1,Nx		
			if(j==1) then
				VisAdvMatrix(istart+j,istart+j)=-Cgx/dx
				VisAdvMatrix(istart+j,istart+j+1)=Cgx/dx			
			else
				VisAdvMatrix(istart+j,istart+j-1)=-Cgx/dx
				VisAdvMatrix(istart+j,istart+j)=Cgx/dx
			end if
			
		end do
	end do
	
	! the derivative in y direction multiplied with group velocity in y direction
	do i=1,Nx
		! y-derivative of laplacian
		do j=1,Ny			
				if(j==Ny) then
					VisAdvMatrix(i+Nx*(j-1),i+Nx*(j-2))=VisAdvMatrix(i+Nx*(j-1),i+Nx*(j-2))-Cgy/dy	
					VisAdvMatrix(i+Nx*(j-1),i+Nx*(j-1))=VisAdvMatrix(i+Nx*(j-1),i+Nx*(j-1))+Cgy/dy			
				else
					VisAdvMatrix(i+Nx*(j-1),i+Nx*(j-1))=VisAdvMatrix(i+Nx*(j-1),i+Nx*(j-1))-Cgy/dy	
					VisAdvMatrix(i+Nx*(j-1),i+Nx*j)=VisAdvMatrix(i+Nx*(j-1),i+Nx*j)+Cgy/dy	
				end if
		end do
	end do
	
	! the viscous dissipation term
	do i=1,Ny		
		do j=1,Nx
			! the counter 
			istart=(i-1)*Nx+j
			
			VisAdvMatrix(istart,istart)=VisAdvMatrix(istart,istart)+alpha
		end do
	end do
	
end subroutine ViscousAdvectionFinite

