!!! this file contains all numerical differentiation routines 
!!! including laplacian and gradient operators in 2-D and 3-D
!!! Jacobiam Matrices are also created here
subroutine FirstDiff(D1,N)
	!! This function returns differentian matrix in 1D on mother interval -1 to 1
	integer i,j,k 
	! the number of grid points
	integer, intent(in):: N
	!the differentiation matrix
	real*8, intent(inout),dimension(N,N):: D1
	! 1D GLL points and corresponding weight w
	real*8, allocatable::x_1D(:),w(:)
	!! Legendre Polynomial and Its derivative
	real*8,allocatable:: Ln(:)
	!!  Its derivative
	real*8 Lpn
	!!  q(x) and Its derivative
	real*8 q,qp
	! 1D GLL points and corresponding weight
	allocate(x_1D(N))
	allocate(w(N))
	allocate(Ln((N)))
	! let's generate GLL points
	call GLLPoints(x_1D,w,N)

	! let's first generate legendre polynomial at Nth order at grid points
	do i=1,N
		call LegendrePolynomials(Ln(i),Lpn,q,qp,x_1D(i),(N-1))
	end do
	
	do i=1,N
		do j=1,N
		! the diagonal elements of the differentiation
			if(i==j) then
			! the first element
				if(i==1) then
					D1(i,j)=-0.25*(N-1)*(N)
					! the last element
				else if(i==N) then
						D1(i,j)=0.25*(N-1)*(N)
						else
							D1(i,j)=0.
				end if
			else 
				D1(i,j)=Ln(i)/(Ln(j)*(x_1D(i)-x_1D(j)))
			end if
		end do
		!print*,D1(i,:)
	end do
end subroutine FirstDiff

subroutine LagrangeInterpolant(lij,xi,j,x,N)
	!!this routine calculates lagrange interpolants
	integer i,k 
	! the number of grid points,  xj the location of lagrange interpolant
	integer, intent(in):: N,j
	!the grid coordinate array
	real*8, intent(inout),dimension(N):: x
	! lagrange interpolant
	real*8, intent(inout)::lij
	! xi the location of lagrange interpolant evaluated
	real*8, intent(inout)::xi
	
	! it is initially set to 1
	lij=1.
	
	do i=1,N
		if(i==j) then 
			continue 
		else 
			lij=lij*(xi-x(i))/(x(j)-x(i))
		end if
	end do
		
end subroutine LagrangeInterpolant

