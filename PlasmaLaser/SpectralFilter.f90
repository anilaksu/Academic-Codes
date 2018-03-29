
!! this routine generates 1-D and 2-D nonlinear convective terms 
subroutine OneDGLLFilter(UF,U,x_grid,w,N)
	!! This function legendre transforms data, filters and backtransforms
	integer i,j,k 
	! the number of grid points
	integer, intent(in):: N
	!the original and filtered data and the grid and the corresponding weight
	real*8, intent(inout),dimension(N):: U,UF,x_grid,w
	! the order and the deviation of the grid
	real*8  sig,order
	!! Legendre Transform of the data
	real*8,allocatable:: ULeg(:),ULegF(:)
	allocate(ULeg(N+1))
	allocate(ULegF(N+1))
	
	sig=10.
	order=4.
	! let's legendre transform the data
 	call LegendreTransform(U,ULeg,x_grid,w,N)
	! here we filter the transform data which basically multiplication of data
	! with filter exp(-((i-1)/sig)^2p)
	do i=1,N+1
		ULegF(i)=ULeg(i)*dexp(-((i-1)/sig)**(2.*order))
		!print*,ULegF(i),ULeg(i),dexp(-((i-1)/sig)**(2.*order))
	end do
	
	! let's back transform the filtered data
	call LegendreBackTransform(UF,ULegF,x_grid,N)

	!open(139,file='LegendreTrans.dat',status='unknown')
	!do i=1,N
		!write(139,*) i-1,ULeg(i)
	!	write(139,*) x_grid(i),UF(i)
	!end do
end subroutine OneDGLLFilter
