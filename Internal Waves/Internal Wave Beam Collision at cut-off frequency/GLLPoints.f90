!!! This file contains all grid generation operations including:
!!! 1) 1D GLL and Fourier Points
!!! 2) 2D GLL-GLL, GLL-Fourier, Fourier-Fourier points
!!! 3) Coordinate Mapping functions like Numerical Jacobian Matrix calculations
subroutine OneDGrid(x_ini,x_end,x_grid,N)
	integer i
	! the number of grid points
	integer, intent(in):: N
	!the start and end coordinates of mesh
	real*8, intent(in):: x_ini,x_end
	!the grid coordinate array
	real*8, intent(inout),dimension(N):: x_grid
	! the length of domain and the uniform distance between points
	real*8 domLength, dx
	
	! the domain length is calculated as
	domLength=x_end -x_ini
	! let's calculate the uniform distance between points
	dx=domLength/(N-1)
	! first point of grid starts with x_ini
	x_grid(1)=x_ini
	
	do i=1,(N-1)
		x_grid(i+1)=x_grid(i)+dx
		!print*, x_grid(i+1)
	end do
		
end subroutine OneDGrid

subroutine GLLPoints(x_grid,w,N)
	!!! this routine generates Gauss-Legendre-Lobatto Points in mother interval
	integer i,j
	!!! the number points
	integer,intent(inout)::N
	!!! GLL points and corresponding weights
	real*8,intent(inout),dimension(N):: x_grid,w
	!!! The order of GLL function 
	integer Nl
	!! maximum number of iterations
	integer maxiter
	!! tolerance for convergence and the root of the function  
	real*8 tol,x_root
	!! the number pi
	real*8 pi
	!! asymptotical root to serve as initial Guess
	real*8,allocatable:: x_asymptotic(:)
	!! it shows if it converged or not
	logical success
	!! Legendre Polynomial and Its derivative
	real*8 Ln,Lpn
	!!  q(x) and Its derivative
	real*8 q,qp
     
	allocate(x_asymptotic(N-2))
    !! the order of GLL function
	Nl=N-1
	!! convergence criterias 
	tol=10E-14
	maxiter=100
	!pi number
	pi=4.*datan(1.d0)
	
	!! let's start finding GLL points
	if(Nl==1) then 
		x_grid(1)=-1.
		x_grid(N)=1.
		w(1)=1.
		w(N)=1.
	!end if
	else 
		! first and last points
	    x_grid(1)=-1.
		x_grid(N)=1.
		w(1)=2./(Nl*(Nl+1.))
		w(N)=2./(Nl*(Nl+1.))
		!print*,w(1),w(N)
		!! since points and weight are symmetric there is no need to generate positive 
		!! and negative side of the function
		do i=1,ceiling(real((Nl)/2))
			! asymptotic relation for initial guess
			x_asymptotic(i)=-dcos((1.*i+0.25)*pi/Nl-(3./(8.*Nl*pi))/(1.*i+0.25))
			call find_root( x_asymptotic(i),Nl,tol, maxiter, x_root, success )
			x_grid(i+1)=x_root
			call LegendrePolynomials(Ln,Lpn,q,qp,x_grid(i+1),Nl)
			w(i+1)=2/(Nl*(Nl+1)*Ln**2.)
		end do
		!! let's fill the positive side
		do i=1,(Nl/2)
			x_grid(N-i)=-1.*x_grid(i+1)
			w(N-i)=w(i+1)
		end do
	end if
	

end subroutine GLLPoints

subroutine LegendrePolynomials(Ln,Lpn,q,qp,x,N)
	!!! this routine is generated to calculate the Legendre Polynomial
	!! this subroutine generates q(x)=LN+1(x)-LN-1(x) to find GLL points
	!!! and its derivative at Nth order on x coordinate
	integer i,j
	!! the order of Legendre Polynomial
	integer, intent(in)::N
	!! the x coordinate 
	real*8, intent(in)::x
	!! Legendre Polynomial and Its derivative
	real*8, intent(inout)::Ln,Lpn
	!!  q(x) and Its derivative
	real*8, intent(inout)::q,qp
	!! the lower order Legendre polynomials and their derivatives 
	real*8 Ln1,Lpn1,Ln2,Lpn2
	!! one higher order Legendre polynomial and its derivative
	real*8 Lnp,Lpnp
	
	!! the convention used above 
	!! Ln1= L^(n-1)
	!! Lpn1= Lp^(n-1)
	!! Ln2= L^(n-2)
	!! Lpn2= Lp^(n-2)

	!! note that zeroth order Legendre Polynomial is L0=1 and Its Derivative Lp0=0
	!! note that first order Legendre Polynomial is L1=x and Its Derivative Lp1=1
	Ln2=1
	Ln1=x
	Lpn2=0
	Lpn1=1
	
	if(N==0) then 
		Ln=1.
		Lpn=0.
		!! q(x)			
		q=0.
		!! qp(x)
		qp=0.
	else if (N==1) then
			Ln=x
			Lpn=1.
			!! q(x)			
			q=0.
			!! qp(x)
			qp=0.
		 else  

			!! Let's start a loop to calculate Nth order Legendre Polynomial and Its Derivative
			do i=1,(N-1)
				!! the recursion relation between legendre polynomial are given as:
				!!Ln+1(x)=(2n+1)/(n+1)xLn(x)-n/(n+1)Ln-1(x) 
				!! therefore Legendre Polynomial
				Ln=(2*i+1)*x*Ln1/(i+1)-i*Ln2/(i+1)
				!! therefore Legendre Polynomial's Derivative
				Lpn=(2*i+1)*(x*Lpn1+Ln1)/(i+1)-i*Lpn2/(i+1)
				!! let's update lower order polynomials 
				Ln2=Ln1
				Ln1=Ln
				!! and their derivatives
				Lpn2=Lpn1
				Lpn1=Lpn
			end do

			!! therefore Legendre Polynomial
			Lnp=(2*N+1)*x*Ln1/(N+1)-N*Ln2/(N+1)
			!! therefore Legendre Polynomial's Derivative
			Lnpm=(2*N+1)*(x*Lpn1+Ln1)/(N+1)-N*Lpn2/(N+1)
			!! q(x)			
			q=Lnp-Ln2
			!! qp(x)
			qp=Lnpm-Lpn2
	end if
end subroutine LegendrePolynomials

subroutine find_root( xinit,N,tol, maxiter, result, success )
	!! this routine is generated to calculate roots of a given function 
	!! Function is called inside the routine therefore for different function 
	!! it has to be modified inside this function
	
	!! initial guess
    real*8, intent(in)   :: xinit
	!! the order of the legendre polynomial
	integer, intent(in)  :: N
	!! tolerance of the convergence
    real*8, intent(in)   :: tol
	!! maximum number of the iterations
    integer, intent(in)  :: maxiter
	!! found root 
    real*8, intent(out)  :: result
	!! it shows if it converged or not
    logical, intent(out) :: success
    !! the derivative of the function
    real*8               :: fprime
    !! the function itself 
	!! legendre polynomials and associated q
	real*8               :: Ln,Lpn,q,qp
    real*8               :: x,xnew
    integer              :: i

    result  = 0.0
    success = .false.

    x = xinit
    do i = 1,max(1,maxiter)
	    !! let's calculate the function at x
		call LegendrePolynomials(Ln,Lpn,q,qp,x,N)
        !! Newton-Raphson iteration
        xnew   = x - q / qp

        if ( abs(xnew-x) <= tol ) then
            success = .true.
            result  = xnew
            exit
        endif

        x = xnew
       ! print*, i, x
     enddo

end subroutine find_root