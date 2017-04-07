!! this routine generates system matrices

subroutine getStiffMatrix(StiffMatrix,Nx,Ny,lambda,mu)
	!! This function returns differentian matrix in 1D on mother interval -1 to 1
	integer i,j,k 
	!! the start and end of indices
	integer istart,iend,jstart,jend
	!! lame coefficients
	real*8,intent(inout):: lambda,mu
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: StiffMatrix
	! 1D GLL points and corresponding weight w in x and y direction
	real*8, allocatable::x_1Dx(:),wx(:),x_1Dy(:),wy(:)
	!! Differentiation Matrix in x and y direction
	real*8,allocatable:: DiffX(:,:),WeakDiffXX(:,:), DiffY(:,:),WeakDiffYY(:,:)
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dx(Nx))
	allocate(wx(Nx))
	allocate(x_1Dy(Ny))
	allocate(wy(Ny))
	! the allocation of first derivative matrix and the product of the derivatices of lagrange interpolants
	allocate(DiffX(Nx*Ny,Nx*Ny))
	allocate(WeakDiffXX(Nx*Ny,Nx*Ny))
	allocate(DiffY(Nx*Ny,Nx*Ny))
	allocate(WeakDiffYY(Nx*Ny,Nx*Ny))
	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dx,wx,Nx)
	call GLLPoints(x_1Dy,wy,Ny)
	
	! let's generate the weak second derivatives
	call getWeakDiffXXMatrix(WeakDiffXX,Nx,Ny)
	call getWeakDiffYYMatrix(WeakDiffYY,Nx,Ny)
	
	! the stiffness matrix
	StiffMatrix=WeakDiffXX+WeakDiffYY
	
	
end subroutine getStiffMatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                     !
!   All differentiation matrices both direct first    !
!   derivatives and weak second derivatives 		  !
!                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getWeakDiffXXMatrix(WeakDiffXX,Nx,Ny)
	!! This function returns the weak second differentiation matrix in x direction
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::WeakDiffXX
	!! DiffX Matrix and Mass Matrix
	real*8,allocatable:: DiffX(:,:),DiffXt(:,:),DiagX(:,:)
	
	
	! the allocation of first derivative matrix 
	allocate(DiffX(Nx*Ny,Nx*Ny))
	! and its transpose
	allocate(DiffXt(Nx*Ny,Nx*Ny))
	! the allocation of the mass matrix
	allocate(DiagX(Nx*Ny,Nx*Ny))
	
	! let's generate the differentiation matrix
	call getDiffXMatrix(DiffX,Nx,Ny)
	! let's generate the mass matrix
	call getDiagYMatrix(DiagX,Nx,Ny)
	! let' get transpose of the differentiation matrix 
	call getTranspose(DiffX,DiffXt,Nx*Ny,Nx*Ny)
	
	! let's generate the weak differentiation matrix
	call matmultnon(DiffXt,DiagX,WeakDiffXX,Nx*Ny,Nx*Ny,Nx*Ny)
	call matmultnon(WeakDiffXX,DiffX,DiffXt,Nx*Ny,Nx*Ny,Nx*Ny)
	
	WeakDiffXX=DiffXt
	
end subroutine getWeakDiffXXMatrix

subroutine getWeakDiffYYMatrix(WeakDiffYY,Nx,Ny)
	!! This function returns the weak second differentiation matrix in y direction
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::WeakDiffYY
	!! DiffX Matrix and Mass Matrix
	real*8,allocatable:: DiffY(:,:),DiffYt(:,:),DiagY(:,:)
	
	
	! the allocation of first derivative matrix 
	allocate(DiffY(Nx*Ny,Nx*Ny))
	! and its transpose
	allocate(DiffYt(Nx*Ny,Nx*Ny))
	! the allocation of the mass matrix
	allocate(DiagY(Nx*Ny,Nx*Ny))
	
	! let's generate the differentiation matrix
	call getDiffYMatrix(DiffY,Nx,Ny)
	! let's generate the mass matrix
	call getDiagYMatrix(DiagY,Nx,Ny)
	! let' get transpose of the differentiation matrix 
	call getTranspose(DiffY,DiffYt,Nx*Ny,Nx*Ny)
	
	! let's generate the weak differentiation matrix
	call matmultnon(DiffYt,DiagY,WeakDiffYY,Nx*Ny,Nx*Ny,Nx*Ny)
	call matmultnon(WeakDiffYY,DiffY,DiffYt,Nx*Ny,Nx*Ny,Nx*Ny)
	
	WeakDiffYY=DiffYt
	
end subroutine getWeakDiffYYMatrix

subroutine getDiffXMatrix(DiffX,Nx,Ny)
	!! This function returns the inverse of the mass matrix and the identity matrix 
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::DiffX
	!! Differentiation Matrix
	real*8,allocatable:: D1x(:,:)
	!! the upper limit and the lower limit and counter
	integer u_limit,l_limit,counter
	
	! the allocation of first derivative matrix 
	allocate(D1x(Nx,Nx))
	! let's first generate differentiation matrix
	call FirstDiff(D1x,Nx)

	! let's set it to zero first then fill it 
	DiffX=0.
	! the mass matrix part
	do i=1,Ny
		do j=1,Nx
			! the counter
			counter=(i-1)*Nx+j
			! the upper limit
			u_limit=i*Nx
			! the lower limit
			l_limit=(i-1)*Nx+1
			! let'fill the differentiation matrix
			DiffX(counter,l_limit:u_limit)=D1x(j,:)
		end do
	end do

end subroutine getDiffXMatrix

subroutine getDiffYMatrix(DiffY,Nx,Ny)
	!! This function returns the inverse of the mass matrix and the identity matrix 
	integer i,j,k
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::DiffY
	!! Differentiation Matrix
	real*8,allocatable:: D1y(:,:)
	!! the upper limit and the lower limit and counter
	integer counter
	
	! the allocation of first derivative matrix 
	allocate(D1y(Ny,Ny))
	! let's first generate differentiation matrix
	call FirstDiff(D1y,Ny)

	! let's set it to zero first then fill it 
	DiffY=0.
	! the mass matrix part
	do i=1,Ny
		do j=1,Nx
			! the counter
			counter=(i-1)*Nx+j
			do k=1,Ny
				! let'fill the differentiation matrix
				DiffY(counter,(k-1)*Nx+j)=D1y(j,k)
			end do
		end do
	end do

end subroutine getDiffYMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                     !
!   Mass Matrix and Diagonal Matrices with weights    !
!                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getInverseMassMatrix(InvMass,Nx,Ny,Lx,Ly)
	!! This function returns the inverse of the mass matrix and the identity matrix 
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	! the domain length in x and y direction
	real*8, intent(in):: Lx,Ly
	!the inverse mass matrix
	real*8, intent(inout), dimension(2*Nx*Ny,2*Nx*Ny)::InvMass
	!! the GLL points in x and y direction and corresponding mass
	real*8,allocatable:: x_1Dx(:),wx(:),x_1Dy(:),wy(:)
	!! the counter
	integer counter
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dx(Nx))
	allocate(wx(Nx))
	allocate(x_1Dy(Ny))
	allocate(wy(Ny))

	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dx,wx,Nx)
	call GLLPoints(x_1Dy,wy,Ny)

	! let's set it to zero first then fill it 
	InvMass=0.
	! the mass matrix part
	do i=1,Ny
		do j=1,Nx
			counter=(i-1)*Nx+j
			InvMass(counter,counter)=1./(wx(j)*wy(i))
		end do
	end do
	! the identity matrix part
	do i=1,Ny
		do j=1,Nx
			counter=(i-1)*Nx+j
			InvMass(counter+Nx*Ny,counter+Nx*Ny)=1.
		end do
	end do
end subroutine getInverseMassMatrix

subroutine getMassMatrix(MassMatrix,Nx,Ny)
	!! This function returns the mass matrix
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::MassMatrix
	!! the GLL points in x and y direction and corresponding mass
	real*8,allocatable:: x_1Dx(:),wx(:),x_1Dy(:),wy(:)
	!! the counter
	integer counter
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dx(Nx))
	allocate(wx(Nx))
	allocate(x_1Dy(Ny))
	allocate(wy(Ny))

	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dx,wx,Nx)
	call GLLPoints(x_1Dy,wy,Ny)

	! the mass matrix part
	do i=1,Ny
		do j=1,Nx
			counter=(i-1)*Nx+j
			MassMatrix(counter,counter)=wx(j)*wy(i)
		end do
	end do

end subroutine getMassMatrix

subroutine getDiagXMatrix(DiagX,Nx,Ny)
	!! This function returns the mass matrix
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::DiagX
	!! the GLL points in x direction and corresponding mass
	real*8,allocatable:: x_1Dx(:),wx(:)
	!! the counter
	integer counter
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dx(Nx))
	allocate(wx(Nx))


	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dx,wx,Nx)


	! the mass matrix part
	do i=1,Ny
		do j=1,Nx
			counter=(i-1)*Nx+j
			DiagX(counter,counter)=wx(j)
		end do
	end do

end subroutine getDiagXMatrix

subroutine getDiagYMatrix(DiagY,Nx,Ny)
	!! This function returns the mass matrix
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(Nx*Ny,Nx*Ny)::DiagY
	!! the GLL points in y direction and corresponding mass
	real*8,allocatable:: x_1Dy(:),wy(:)
	!! the counter
	integer counter
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dy(Ny))
	allocate(wy(Ny))

	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dy,wy,Ny)

	! the mass matrix part
	do i=1,Ny
		do j=1,Nx
			counter=(i-1)*Nx+j
			DiagY(counter,counter)=wy(i)
		end do
	end do

end subroutine getDiagYMatrix


