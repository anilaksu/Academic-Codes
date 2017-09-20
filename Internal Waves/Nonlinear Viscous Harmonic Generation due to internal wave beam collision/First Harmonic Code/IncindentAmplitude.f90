!! this routine generates the amplitude of the incident and the reflecting internal wave beam
!! it also performs the required coordinate transformation
subroutine AmplitudeIncidentX(x,z,Cgx,Cgz)
	!! This function generates boundary conditions 
	!! if LB==1, it is neumann boundary condtion 
	!! if LB==0, it is dirichlet boundary condition
	integer i,j,k 
	! the number of grid points and left boundary condition LB and right boundary condtion RB
	integer, intent(in):: N, LB, RB
	!the differentiation matrix
	real*8, intent(inout),dimension(N,N):: BCMatrix
	!! Differentiation Matrix 
	real*8,allocatable:: D1(:,:)

	

	
end subroutine AmplitudeIncidentX



