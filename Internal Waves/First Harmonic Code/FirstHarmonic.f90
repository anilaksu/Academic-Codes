
Program FirstHarmonic
	!this program is written to solve 2D first harmonic generation equation 
	! The full details of the equation will be published on Physics of Fluids Journal by Anil Aksu
	
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
	! the reference points for the incident and the reflecting internal wave beam
	real*8, allocatable:: xi_ref(:),xr_ref(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!               Dependent Variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the velocity field for the incident and the reflecting internal wave beam
	real*8, allocatable::u_inc(:,:),v_inc(:,:),u_ref(:,:),v_ref(:,:)
	! the amplitude field for the incindent and the reflecting internal wave beam
	real*8, allocatable::Ampx_inc(:),Ampz_inc(:),Ampx_ref(:),Ampz_ref(:)
	! the velocity field for the first harmonic wave field
	real*8, allocatable::u_har(:,:),v_har(:,:)
	! the famplitude for the first harmonic wave field
	real*8, allocatable::Ampx_har(:),Ampz_har(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The dimensions of dependent and 			  !
	!             independent variables                   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the number of points in grid and the order of legendre polynomial
	integer N,Nl,Ngrid
	! the number of subgrids and the number of time steps
	integer Nsub,Ntime
	! the length of domain in each direction
	real*8 :: Lx,Ly
	! the number of points in each subgrid in x and y direction
	integer Nx,Ny
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Material Properties 			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the kinematic viscosity and dynamic viscosity
	real*8 xnu, xmu
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!             Internal Wave Parameters			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the buoyancy frequency
	real*8 BV
	! the group velocities and the viscous dissipation coefficients and amplitude half-width and the wave length in x direction
	real*8 Cgx,Cgy,alpha,sigma,lambdax
	! the viscous dissipation for first harmonic and the wave numbers and the inital amplitude of the incident beam
	real*8 alpha_2,kx,kz,A_0
	! the parameter for the first higher harmonic
	real*8 kx_har,kz_har,omega_har
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Numerical Operators			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 1-D Discontinuous SEM Laplacian Operator and Boundary Condition Matrix
	real*8, allocatable:: LapDis1D(:,:),SEMBCMatrix1D(:,:)
	! the identity matrix 
	real*8, allocatable:: ID(:,:),ID2D(:,:)
	! 2-D Laplacian Operator and Boundary Condition Matrix
	real*8, allocatable:: VisAdvMatrix(:,:),BCMatrix2D(:,:)
	
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
	! Test matrix 
	real*8, allocatable:: A(:,:),Ainv(:,:)
	! Boundary condition, Governing Equation and System Matrix and inverse of the system matrix
	real*8, allocatable:: BCMatrix(:,:),GEMatrix(:,:),SysMatrix(:,:),InvSysMatrix(:,:)
	! Test Function
	real*8, allocatable:: F(:),dF(:),FB(:),LegenData(:)
	! the right hand side of the equation
	real*8, allocatable:: RHS(:)
	! the differentiation matrix
	real*8, allocatable:: D1(:,:)

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
	! The Material Properties and Internal Wave Parameters!
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the kinematic viscosity
	xnu=10E-1
	! the group velocities
	Cgx=0.05
	Cgy=-0.05
	! the viscous dissipation coefficient
	alpha=0.0002
	
	!the wave numbers and 
	kx=0.0875
	kz=-0.0875
	!the inital amplitude of the incident beam
	A_0=0.3
	! the wave length of the incident internal wave beam in x direction 
	lambdax=2.*pi/kx
	! amplitude half-width
	sigma=0.5*lambdax
	
	! the first harmonic frequency and wave numbers
	kx_har=2.*kx
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!    2D Steady Advection Convection Equation 		  !
	!			Solution on GLL-GLL Grid     			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! domian properties: number of points in domain and length of domain
	Nx=20
	Ny=20
	! the domain length in x direction and y direction
	Lx=10.*lambdax
	Ly=5.*lambdax
	
	! let's allocate all the related matrices
	allocate(BCMatrix2D(Nx*Ny,Nx*Ny))
	allocate(VisAdvMatrix(Nx*Ny,Nx*Ny))
	allocate(SysMatrix(Nx*Ny,Nx*Ny))
	allocate(InvSysMatrix(Nx*Ny,Nx*Ny))
	
	allocate(F(Nx*Ny))
	allocate(RHS(Nx*Ny))
	allocate(x_2Dmother(Nx*Ny,2))
	allocate(x_2Dreal(Nx*Ny,2))
	allocate(wx(Nx))
	allocate(wy(Ny))
	
	! the indicent beam parameters
	allocate(Ampx_inc(Nx*Ny))
	allocate(Ampz_inc(Ny*Ny))
	! the reflecting beam parameters
	allocate(Ampx_ref(Nx*Ny))
	allocate(Ampz_ref(Ny*Ny))
	! the reference points for the indicent and the reflecting internal wave beams
	allocate(xi_ref(2))
	allocate(xr_ref(2))
	
	do i=1,Ny
		do j=1,Nx
			if(j==1 .and. i>10 .and. i<12) then
				F((i-1)*Nx+j)=20.
			else
				F(i)=0.
			end if		
		end do
	end do
	
	! 2D boundary cond
	call TwoDBCFirst(BCMatrix2D,0,1,1,0,Nx,Ny)
	call ViscousAdvection(VisAdvMatrix,Nx,Ny,Lx,Ly,Cgx,Cgy,alpha)
	call TwoDBCApplyFirst(BCMatrix2D,VisAdvMatrix,SysMatrix,Nx,Ny)
	call invertSVD(Nx*Ny,SysMatrix,InvSysMatrix)
	call TwoDRHS(RHS,F,Nx,Ny)
	call matvect(InvSysMatrix,RHS,F,Nx*Ny)
	! let's generate the mother grid
	call TwoDMapXY(x_2Dmother,wx,wy,Nx,Ny)
	! let's generate the real grid
	do i=1,Nx*Ny
		x_2Dreal(i,1)=(x_2Dmother(i,1)+1)*Lx/2.
		x_2Dreal(i,2)=(x_2Dmother(i,2)+1)*Ly/2.
	end do
	
	! the reference points 
	! the incident reference
	xi_ref(1)=0.
	xi_ref(2)=0.
	! the reflecting reference
	xr_ref(1)=5.*lambdax
	xr_ref(2)=5.*lambdax
	
	! let's get the incident internal wave beam field
	do i=1,Nx*Ny
		call getIncidentAmplitude(Ampx_inc(i),x_2Dreal(i,:),xi_ref,sigma, alpha, A_0, kx, kz)
		call getIncidentAmplitude(Ampx_ref(i),x_2Dreal(i,:),xr_ref,sigma, alpha, A_0, kx, -1.*kz)
	end do
	
	open(140,file='2DLaplace.dat',status='unknown')
	open(141,file='2DNumerical.dat',status='unknown')
	open(142,file='IncidentAmplitude.dat',status='unknown')
	open(143,file='ReflectingAmplitude.dat',status='unknown')
	
	do i=1,Nx*Ny
	!	write(140,*) BCMatrix2D(i,:)
	!	write(140,*) SysMatrix(i,:)
		write(141,*) x_2Dreal(i,:)/lambdax,F(i)
		write(142,*) x_2Dreal(i,:)/lambdax,Ampx_inc(i)
		write(143,*) x_2Dreal(i,:)/lambdax,Ampx_ref(i)
	end do
		
		
	
			
End Program FirstHarmonic	