!!! this file contains all numerical integration routines both in 1-D, 2-D and 3-D
!!! It also integrates on complex geometries
subroutine OneDIntegration(int_f,f,w,N)
	integer i
	! the number of grid points
	integer, intent(in):: N
	!the start and end coordinates of mesh
	real*8, intent(inout):: int_f
	!the grid coordinate array
	real*8, intent(inout),dimension(N):: f,w
	
	!! integral must start from 0
	int_f=0.
	
	do i=1,N
		int_f=int_f+w(i)*f(i)
	end do
		
end subroutine OneDIntegration