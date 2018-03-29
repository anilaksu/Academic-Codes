!! this routine generates stiffness matrices for Non-linear Scrodinger Equation developed by Anil Aksu

!subroutine getStateSpace2ndDiff(StateMatrix,LaplaceMatrix,Nr,Nz)
	!! This function returns differentian matrix in 1D on mother interval -1 to 1
!	integer i,j,k 
	!! the start and end of indices
!	integer istart,iend,jstart,jend
	! the number of grid points
!	integer, intent(in):: Nx,Ny
	!the weak laplacian 
!	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: LaplaceMatrix
	!the state space matrix 
!	real*8, intent(inout),dimension(2*Nx*Ny,2*Nx*Ny):: StateMatrix
	!! the counter
!	integer counter
	
	! the state space matrix 
!	do i=1,Ny
!		do j=1,Nx
!			counter=(i-1)*Nx+j
			! the stiffness part 
!			StateMatrix(counter,counter)=LaplaceMatrix(i,j)
			! the zeroth derivative mass matrix part
!			StateMatrix(counter+Nx*Ny,counter+Nx*Ny)=1.
!		end do
!	end do

!end subroutine getStateSpace2ndDiff

subroutine getStiffMatrix(DervMatrix,Nr,Nz,Lr,Lz)
	!! This function returns differentian matrix in 2D 
	integer i,j,k 
	!! the start and end of indices
	integer istart,iend,jstart,jend
	! the number of grid points
	integer, intent(in):: Nr,Nz
	! the domain length in each direction 
	real*8, intent(in)::Lr,Lz
	!the differentiation matrix
	real*8, intent(inout),dimension(Nr*Nz,Nr*Nz):: DervMatrix
	! 1D GLL points and corresponding weight w in r and z direction
	real*8, allocatable::x_1Dr(:),wr(:),x_1Dz(:),wz(:)
	!! Differentiation Matrix and the product of the derivatices of lagrange interpolants in x and y direction
	real*8,allocatable:: D1r(:,:),ldpnr(:,:), D1z(:,:),ldpnz(:,:),Diffz(:,:)
	
	! 1D GLL points and corresponding weight
	allocate(x_1Dr(Nr))
	allocate(wr(Nz))
	allocate(x_1Dz(Nz))
	allocate(wz(Nz))
	! the allocation of first derivative matrix and the product of the derivatices of lagrange interpolants
	allocate(D1r(Nr,Nr))
	allocate(ldpnr(Nr,Nr))
	allocate(D1z(Nz,Nz))
	allocate(ldpnz(Nz,Nz))
	! the derivative matrix in z direction
	allocate(Diffz(Nr*Nz,Nr*Nz))
	! let's first generate grid points and corresponding weight
	call GLLPoints(x_1Dr,wr,Nr)
	call GLLPoints(x_1Dz,wz,Nz)
	! product of lagrange interpolants
	call Lagder(ldpnr,Nr)
	call Lagder(ldpnz,Nz)
	
	do i=1,Nz
		!do j=1,N
		istart=Nr*(i-1)
		iend=Nr*i
		! r-derivative of laplacian
		do j=1,Nr
			do k=1,Nr
				DervMatrix(istart+j,istart+k)=-1.*wr(i)*ldpnr(j,k)*0.25/(Lr**2.)
				!print*,LaplaceMatrix(istart+j,istart+k)
			end do 
		end do
	end do
	! let's generate the derivative in z direction 
	call getDiffZMatrix(Diffz,Nr,Nz,Lz)
	
	do i=1,Nz*Nr
		do j=1,Nz*Nr
			DervMatrix(i,j)=DervMatrix(i,j)+Diffz(i,j)
		end do
	end do
	
end subroutine getStiffMatrix

subroutine getDiffZMatrix(Diffz,Nr,Nz,Lz)
	!! This function returns differentian matrix in 2D domain with Lz length in z direction 
	integer i,j,k 
	!! the start and end of indices
	integer istart,iend,jstart,jend
	! the number of grid points
	integer, intent(in):: Nr,Nz
	! the domain length in x direction
	real*8, intent(in)::Lz
	!the differentiation matrix
	real*8, intent(inout),dimension(Nr*Nz,Nr*Nz):: DiffZ
	!! Differentiation Matrix in x  directioın
	real*8,allocatable:: D1z(:,:)

	! the allocation of first derivative matrix and the product of the derivatices of lagrange interpolants
	allocate(D1z(Nz,Nz))
	! let's first generate differentiation matrix
	call FirstDiff(D1z,Nz)
	
	do i=1,Nz
		do j=1,Nr
			do k=1,Nz
				istart=Nr*(i-1)+j+Nr*(k-1)
				jstart=Nr*(i-1)+j
				! z-derivative 
				DiffZ(istart,jstart)=D1z(i,k)*2./Lz
			end do
		end do
	end do
	
end subroutine getDiffZMatrix

subroutine Lagder(ldpn,N)
	!! it calculate product matrice of ldplpn
	integer i,j
	! the number of grid points
	integer, intent(in):: N
	!the differentiation matrix
	real*8, intent(inout),dimension(N,N):: ldpn
	!! Differentiation Matrix
	real*8,allocatable:: D1(:,:)
	! 1D GLL points and corresponding weight w
	real*8, allocatable::x_1D(:),w(:)
	! the allocation of first derivative matrix
	allocate(D1(N,N))
	! 1D GLL points and corresponding weight
	allocate(x_1D(N))
	allocate(w(N))
	! let's first generate differentiation matrix
	call FirstDiff(D1,N)
	! let's generate GLL points
	call GLLPoints(x_1D,w,N)
	ldpn=0.
	
	do i=1,N
		do j=1,N	
			ldpn(i,j)=0.
			do k=1,N
				ldpn(i,j)=ldpn(i,j)+D1(k,i)*D1(k,j)*w(k)
			end do
		!	print*,ldpn(i,j)
		end do
	end do
	
end subroutine Lagder


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                     !
!   Mass Matrix and Diagonal Matrices with weights    !
!                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine getMassMatrixWave(MassMatrix,Nx,Ny)
	!! This function returns the mass matrix
	integer i,j
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the inverse mass matrix
	real*8, intent(inout), dimension(2*Nx*Ny,2*Nx*Ny)::MassMatrix
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
			! the first derivative mass matrix part 
			MassMatrix(counter,counter+Nx*Ny)=wx(j)*wy(i)
			! the zeroth derivative mass matrix part
			MassMatrix(counter+Nx*Ny,counter)=1.
		end do
	end do

end subroutine getMassMatrixWave

