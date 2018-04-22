!! this routine generates 1-D and 2-D Laplacian Operators
subroutine Laplace(LaplaceMatrix,Nx,Ny)
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
end subroutine Laplace

subroutine OneDDisLaplace(LapDis1D,dx,Nsub,Ngrid)
	!! this function calculates the discontinous laplace operator in 1D with Nsub subdomain with 
	!! Ngrid grid points in each domain
	integer i,j
	!! the start and end indices
	integer jstart, jend
	! the number of grid points
	integer, intent(in):: Nsub,Ngrid
	!the laplacian matrix
	real*8, intent(inout),dimension(Nsub*Ngrid,Nsub*Ngrid):: LapDis1D
	!the length of each domain
	real*8, intent(inout),dimension(Nsub):: dx
	!! Single Domain Laplacian in 1-D
	real*8,allocatable:: Lap(:,:)
	!! Differentiation Matrix
	real*8,allocatable:: D1(:,:)
	! the allocation of first derivative matrix and the laplacian
	allocate(D1(Ngrid,Ngrid))
	allocate(Lap(Ngrid,Ngrid))
	! let's generate the single domain laplacian
 	call Lagder(Lap,Ngrid)
	! let's first generate differentiation matrix
	call FirstDiff(D1,Ngrid)
	LapDis1D=0.
	do i=1,Nsub
		jstart=(i-1)*Ngrid+1
		jend=i*Ngrid
		!! let's fill it
		LapDis1D(jstart:jend,jstart:jend)=0.5*dx(i)*Lap
	end do
	! let's add interface fluxes
	do i=1,Nsub-1
		jstart=(i-1)*Ngrid+1
		jend=i*Ngrid
		!! let's add the interface fluxes it
		LapDis1D(jend,jstart:jend)=LapDis1D(jend,jstart:jend)+0.25*dx(i)*D1(Ngrid,:)
		LapDis1D(jend+1,jstart:jend)=LapDis1D(jend+1,jstart:jend)+0.25*dx(i)*D1(Ngrid,:)
		!! the next element
		LapDis1D(jend,jend+1:jend+Ngrid)=LapDis1D(jend,jend+1:jend+Ngrid)-0.25*dx(i)*D1(1,:)
		LapDis1D(jend+1,jend+1:jend+Ngrid)=LapDis1D(jend+1,jend+1:jend+Ngrid)-0.25*dx(i)*D1(1,:)
	end do
end subroutine OneDDisLaplace

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