
Program BrainTumor
	! this program is written to solve 2D transient elastodynamic model with anomalies on the bounries
	! it will be excited on two points seperated by distance d, the hope is to mimic parallax effect to 
	! to detect anomalies
	
	integer i,j,k,jstart,jend
	! the start and the end point of the grid 
	real*8 x_ini, x_end
	! the grid points array and corresponding weight
	real*8, allocatable:: x_grid(:),w(:),x_total(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     Independent Variables  (Coordinate)             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the test grid in coordinate transformation
	real*8, allocatable:: x_2Dreal(:,:),x_2Dmother(:,:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The dimensions of dependent and 			  !
	!             independent variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the number of points in grid and the order of legendre polynomial
	integer N,Nl,Ngrid
	! the number of subgrids 
	integer Nsub
	! the length of domain in each direction
	real*8 :: Lx,Ly
	! the number of points in each subgrid in x and y direction
	integer Nx,Ny
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The time integration parameters             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the time integration parameter and the step size
	real*8 theta, dt
	! the number of time steps
	integer Ntime
	! the explicit part of integration
	real*8, allocatable:: f_time(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!               Dependent Variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! displacement field in x and y direction and theit time derivatives
	real*8, allocatable::u_x(:,:),v_x(:,:),up_x(:,:),vp_x(:,:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!              Elastisity Parameters			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! lame constants and density
	real*8 mu, lambda, rho
	! the excitation wave length 
	real*8 lambdax
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Numerical Operators			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 1-D Discontinuous SEM Laplacian Operator and Boundary Condition Matrix
	real*8, allocatable:: MassMatrix(:,:),InvMassMatrix(:,:)
	! the system matrix and its inverse and the right handside matrix
	real*8, allocatable:: SysMatrix(:,:),InvSysMatrix(:,:),RhsMatrix(:,:)
	! 2-D Laplacian Operator and Boundary Condition Matrix
	real*8, allocatable:: StiffMatrix(:,:),BCMatrix2D(:,:)
	! Boundary condition, Governing Equation and System Matrix and inverse of the system matrix
	real*8, allocatable:: BCMatrix(:,:),GEMatrix(:,:)
	!! Differentiation Matrix and the product of the derivatices of lagrange interpolants in x and y direction
	real*8 ,allocatable:: D1x(:,:),ldpnx(:,:), D1y(:,:),ldpny(:,:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        2-D Grid and related parameters			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the grid points
	real*8, allocatable:: x_grid2D(:,:),wx(:),wy(:)
	
	
	! the number of points on the surface
	integer, allocatable::Nsurf(:)
	! legendre polynomials and associated q
	real*8, allocatable:: Ln(:),Lpn(:),q(:),qp(:)
	! Test Function
	real*8, allocatable:: F(:),dF(:),FB(:),LegenData(:)
	! the right hand side of the equation
	real*8, allocatable:: RHS(:)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                  Useful Constants                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! pi number
	real*8 pi
	! for complex exponential representation
	complex phi
	pi=4.*datan(1.d0)
	phi=CMPLX(0,2.*pi)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!               Elasticity Parameters                 !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! lame parameters
	mu=10 
	lambda=20
	! the excitation wave length in x direction
	lambdax=0.875
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!      2D Transient Elastodynamic Equation 		      !
	!			Solution on GLL-GLL Grid     			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! domian properties: number of points in domain and length of domain
	Nx=40
	Ny=40
	! the number of time steps
	Ntime = 100
	! the domain length in x direction and y direction
	Lx=10.*lambdax
	Ly=5.*lambdax
	
	! let's allocate the grid
	allocate(x_2Dmother(Nx*Ny,2))
	allocate(x_2Dreal(Nx*Ny,2))
	allocate(wx(Nx))
	allocate(wy(Ny))
	

	! let's generate the mother grid
	call TwoDMapXY(x_2Dmother,wx,wy,Nx,Ny)
	! let's generate the real grid
	do i=1,Nx*Ny
		x_2Dreal(i,1)=(x_2Dmother(i,1)+1)*Lx/2.
		x_2Dreal(i,2)=(x_2Dmother(i,2)+1)*Ly/2.
		print*,x_2Dmother(i,:)
	end do
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!  System Matrices: Mass Matrix and Stiffness Matrix  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	! let's allocate all the related matrices
	allocate(BCMatrix2D(Nx*Ny,Nx*Ny))
	allocate(StiffMatrix(Nx*Ny,Nx*Ny))
	allocate(SysMatrix(Nx*Ny,Nx*Ny))
	allocate(InvSysMatrix(Nx*Ny,Nx*Ny))
	allocate(RhsMatrix(Nx*Ny,Nx*Ny))
	allocate(MassMatrix(Nx*Ny,Nx*Ny))
	allocate(InvMassMatrix(Nx*Ny,Nx*Ny))
	! let's allocate all relavent vector
	allocate(u_x(Nx*Ny,Ntime))
	allocate(f_time(Nx*Ny))
	
	! let's generate the mass matrix
	!call getInverseMassMatrix(InvMassMatrix,Nx,Ny,Lx,Ly)
	
	! it is set to zero at boundaries
	do i=1,Ny
		do j=1,Nx
			if(j==1) then
				!F((i-1)*Nx+j)=dsin(2.*pi*x_2Dreal((i-1)*Nx+j,2)/lambdax)
				u_x((i-1)*Nx+j,1)=2.
			else 
				u_x((i-1)*Nx+j,1)=0.
			end if		
		end do
	end do
	print*,"okey"
	
	! let's generate 2D stiffness matrix 
	call getStiffMatrix(StiffMatrix,Nx,Ny,wx,wy,ldpnx,ldpny)
	! let's generate the boundary conditions
	call TwoDBC(BCMatrix2D,0,0,0,0,Nx,Ny)
	! let's apply boundary conditions
	call TwoDBCApplyFirst(BCMatrix2D,StiffMatrix,InvSysMatrix,Nx,Ny)
	! let's generate the required mass matrix
	call getMassMatrix(InvMassMatrix,Nx,Ny)
	! the boundary conditions are also applied to mass matrix in this formulation but this time all of them has to zero
	call TwoDBC(BCMatrix2D,0,0,0,0,Nx,Ny)
	! let's apply them 
	call TwoDBCApplyFirst(BCMatrix2D,InvMassMatrix,MassMatrix,Nx,Ny)
	
	
	! the time integration Crank-Nicholson method is employed to get more stable and accurate 
	! solution, it is commonly referred semi implicit time integration as it has explicit and 
	! implicit part. The calculation of explicit part is straightforward, however implicit part
	! requires matrix inversion
	
	! theta integration parameter
	theta=0.5
	
	! the step size
	dt=1.

	! the RHS term integration matrix
	call getRHSMatrix(MassMatrix,InvSysMatrix,RhsMatrix,Nx,Ny,theta,dt)
	! the integration of RHS
	call getImplicitIntegrate(MassMatrix,InvSysMatrix,SysMatrix,Nx,Ny,theta,dt)
	! let's invert it
	call invertSVD(Nx*Ny,SysMatrix,InvSysMatrix)
	
	do i=2,Ntime
		! the calculation of RHS of time integration, therefore the external forcing of time integration 
		call matvect(RhsMatrix,u_x(:,(i-1)),f_time,Nx*Ny)
		! the calculation of implicit part of integration 
		call matvect(InvSysMatrix,f_time,u_x(:,i),Nx*Ny)		
	end do
	
	!open(139,file='BCMatrix.dat',status='unknown')
	!open(140,file='StiffMatrix.dat',status='unknown')
	open(141,file='2DNumerical.dat',status='unknown')
	
	do i=1,Nx*Ny
	!	write(140,*) BCMatrix2D(i,:)
	!	write(139,*) BCMatrix2D(i,:)
	! 	write(140,*) SysMatrix(i,:)
		write(141,*) x_2Dreal(i,:)/lambdax,u_x(i,Ntime)
	end do
		
		
	
			
End Program BrainTumor