!! this routine generates stiffness matrices for discontinous galerkin method
subroutine getDisStiffMatrix(DisStiffMatrix2D,dy,dx,Nsub,Nx,Ny)
	!! this function calculates the discontinous laplace operator in 2D with Nsub subdomain in x direction 
	!! Nx and Ny grid points in each domain
	integer i,j
	!! the start and end indices
	integer jstart, jend
	! the number of grid points
	integer, intent(in):: Nsub,Nx,Ny
	!the laplacian matrix
	real*8, intent(inout),dimension(Nsub*Nx*Ny,Nsub*Nx*Ny):: DisStiffMatrix2D
	!the length of each domain
	real*8, intent(inout),dimension(Nsub):: dx
	!! Single Domain Laplacian in 2-D
	real*8,allocatable:: StiffMatrix(:,:)
	!! Differentiation Matrix in x direction
	real*8,allocatable:: D1(:,:)
	! the allocation of first derivative matrix and the laplacian
	allocate(D1(Nx,Nx))
	allocate(StiffMatrix(Nx*Ny,Nx*Ny))
	! let's generate the single domain laplacian
 	call getStiffMatrix(StiffMatrix,Nx,Ny)
	! let's first generate differentiation matrix
	call FirstDiff(D1,Nx)
	
	LapDis1D=0.
	do i=1,Nsub
		jstart=(i-1)*Nx*Ny+1
		jend=i*Nx*Ny
		!! let's fill it
		DisStiffMatrix2D(jstart:jend,jstart:jend)=0.5*dx(i)*StiffMatrix
	end do
	! let's add interface fluxes
	!do i=1,Nsub-1
	!	jstart=(i-1)*Ngrid+1
	!	jend=i*Ngrid
		!! let's add the interface fluxes it
	!	LapDis1D(jend,jstart:jend)=LapDis1D(jend,jstart:jend)+0.25*dx(i)*D1(Ngrid,:)
	!	LapDis1D(jend+1,jstart:jend)=LapDis1D(jend+1,jstart:jend)+0.25*dx(i)*D1(Ngrid,:)
		!! the next element
	!	LapDis1D(jend,jend+1:jend+Ngrid)=LapDis1D(jend,jend+1:jend+Ngrid)-0.25*dx(i)*D1(1,:)
	!	LapDis1D(jend+1,jend+1:jend+Ngrid)=LapDis1D(jend+1,jend+1:jend+Ngrid)-0.25*dx(i)*D1(1,:)
	!end do
end subroutine getDisStiffMatrix

subroutine getStiffMatrix(LaplaceMatrix,Nx,Ny)
	!! This function returns differentian matrix in 1D on mother interval -1 to 1
	integer i,j,k 
	!! the start and end of indices
	integer istart,iend,jstart,jend
	! the number of grid points
	integer, intent(in):: Nx,Ny
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: LaplaceMatrix
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
	! product of lagrange interpolants
	call Lagder(ldpnx,Nx)
	call Lagder(ldpny,Ny)

	do i=1,Ny
		!do j=1,N
		istart=Nx*(i-1)
		iend=Nx*i
		! x-derivative of laplacian
		do j=1,Nx
			do k=1,Nx
				LaplaceMatrix(istart+j,istart+k)=-1.*wx(i)*ldpnx(j,k)
				!print*,LaplaceMatrix(istart+j,istart+k)
			end do 
		end do
	end do

	do i=1,Nx
		! y-derivative of laplacian
		do j=1,Ny
			do k=1,Ny
				LaplaceMatrix(i+Nx*(j-1),i+Nx*(k-1))=LaplaceMatrix(i+Nx*(j-1),i+Nx*(k-1))-1.*wy(i)*ldpny(j,k)
				!print*,LaplaceMatrix(istart+j,istart+k)
			end do 
		end do
	end do
end subroutine getStiffMatrix



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

