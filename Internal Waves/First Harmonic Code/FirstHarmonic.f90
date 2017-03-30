
Program FirstHarmonic
	!this program is written to solve 2D first harmonic generation equation 
	! The full details of the equation will be published on Physics of Fluids Journal by Anil Aksu
	
	integer i,j,k,jstart,jend
	! the start and the end point of the grid 
	real*8 x_ini, x_end
	! the time to output data
	real*8 t_out
	! the grid points array and corresponding weight
	real*8, allocatable:: x_grid(:),w(:),x_total(:)
	! the counter
	integer counter
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     Independent Variables  (Coordinate)             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the test grid in coordinate transformation
	real*8, allocatable:: x_2Dreal(:,:),x_2Dmother(:,:)
	! the reference points for the incident and the reflecting internal wave beam
	real*8, allocatable:: xi_ref(:),xr_ref(:)
	! the coordinates to be used in interpolation
	real*8, allocatable::x_intData(:,:)
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
	real*8, allocatable::u_har(:),v_har(:)
	! the famplitude for the first harmonic wave field
	real*8, allocatable::Ampx_har(:),Ampz_har(:),E_har(:)
	! the velocity field for the first harmonic wave field in interpolation grid
	real*8, allocatable::u_int(:),v_int(:)
	! the famplitude for the first harmonic wave field in interpolation grid
	real*8, allocatable::Ampx_int(:),Ampz_int(:),E_int(:)
	! the original data to be used in interpolation
	real*8, allocatable::Ampx_intData(:),Ampz_intData(:),E_intData(:)
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
	! the number of points in interpolation grid 
	integer Nx_int,Ny_int
	! the indice array to be used in interpolation
	integer, allocatable::Indices(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                Material Properties 			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the kinematic viscosity and dynamic viscosity and density
	real*8 xnu, xmu, rho
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!             Internal Wave Parameters			      !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the buoyancy frequency
	real*8 BV, omega
	! the group velocities and the viscous dissipation coefficients and amplitude half-width and the wave length in x direction
	real*8 Cgx,Cgy,alpha,sigma,lambdax
	! the viscous dissipation for first harmonic and the wave numbers and the inital amplitude of the incident beam
	real*8 alpha_2,kx,kz,A_0
	! the parameter for the first higher harmonic
	real*8 kx_har,kz_har,omega_har,Cgx_har,Cgy_har,alpha_har,phi_har, theta_har
	
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
	! the interpolation grid and relavent interpolation matrix
	real*8, allocatable:: x_int2D(:,:),x_rot2D(:,:),IntMatrix(:,:)
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
	xnu=10E-7
	! the density
	rho=1000.
	!the wave number in x direction
	kx=2.*pi/0.875
	!the inital amplitude of the incident beam
	A_0=0.3
	! the wave length of the incident internal wave beam in x direction 
	lambdax=2.*pi/kx
	! amplitude half-width
	sigma=0.4014
	! the bouyancy frequency 
	BV=2.34
	! the incident internal wave frequency 
	omega=0.42*BV
	! the first harmonci frequency
	omega_har=2.*omega
	! the first harmonic frequency and wave numbers
	kx_har=2.*kx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!    2D Steady Advection Convection Equation 		  !
	!			Solution on GLL-GLL Grid     			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! domian properties: number of points in domain and length of domain
	Nx=10
	Ny=10
	! the domain length in x direction and y direction
	Lx=25.*lambdax
	Ly=4.3*lambdax
	
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
	! the indicent beam parameters
	allocate(Ampx_har(Nx*Ny))
	allocate(Ampz_har(Ny*Ny))
	! the reference points for the indicent and the reflecting internal wave beams
	allocate(xi_ref(2))
	allocate(xr_ref(2))
	
	! the reference points 
	! the incident reference
	xi_ref(1)=2.*lambdax
	xi_ref(2)=0.
	
	! let's calculate the first harmonic wave number in z direction
	call getWaveNumberInZ(BV,omega_har, kx_har, kz_har)
	! let's calculat the group velocities for the incident internal wave beam
	call getGroupVelocity(Cgx,Cgy,BV,omega,kx,kz)
	! let's calculate the reflection point
	call getReflectionPoint(xr_ref,xi_ref,Cgx,Cgy,Lx,Ly)
	! let's generate the mother grid
	!call TwoDMapXY(x_2Dmother,wx,wy,Nx,Ny)
	call TwoDMapXYUniform(x_2Dmother,Nx,Ny)
	! let's generate the real grid
	do i=1,Nx*Ny
		x_2Dreal(i,1)=(x_2Dmother(i,1)+1)*Lx/2.
		x_2Dreal(i,2)=(x_2Dmother(i,2)+1)*Ly/2.
	end do
	! let's calculate viscous dissipation coefficient for the first higher harmonic signal
	alpha=xnu*(BV**2.)*(kx**2.)/(2.*dsqrt(Cgx**2.+Cgy**2.)*(omega**2.))
	! let's get the incident internal wave beam field
	do i=1,Nx*Ny
		call getIncidentAmplitude(Ampx_inc(i),x_2Dreal(i,:),xi_ref,sigma, alpha, A_0, kx, kz)
		call getIncidentAmplitude(Ampx_ref(i),x_2Dreal(i,:),xr_ref,sigma, alpha, A_0, kx, -1.*kz)
		F(i)=-2.*(((kx_har*omega)**2.)/(kx*BV**2.))*Ampx_inc(i)*Ampx_ref(i)*dsin(kz_har*x_2Dreal(i,2)-2.*kz*Ly)
		!print*,F(i)
	end do
	
	open(142,file='IncidentAmplitude.dat',status='unknown')
	open(143,file='ReflectingAmplitude.dat',status='unknown')
	open(144,file='RHSForcing.dat',status='unknown')
	
	do i=1,Nx*Ny
		write(142,*) x_2Dreal(i,:)/lambdax,Ampx_inc(i)/A_0
		write(143,*) x_2Dreal(i,:)/lambdax,Ampx_ref(i)/A_0
		write(144,*) x_2Dreal(i,:)/lambdax,F(i)/(A_0**2.)
	end do
	
	! it is set to zero at boundaries
	do i=1,Ny
		do j=1,Nx
			if(j==1 .or. i==Ny) then
				F((i-1)*Nx+j)=0.
			end if		
		end do
	end do
	
	! 2D boundary cond
	call TwoDBCFirst(BCMatrix2D,0,1,1,0,Nx,Ny)
	! let's calculate the group velocities for the first harmonic frequency 
	call getGroupVelocity(Cgx_har,Cgy_har,BV,omega_har,kx_har,kz_har)
	! let's calculate viscous dissipation coefficient for the first higher harmonic signal
	alpha_har=xnu*(BV**2.)*(kx**2.)/(2.*dsqrt(Cgx_har**2.+Cgy_har**2.)*(omega**2.))
	!print*,"the group velocity of first harmonic",Cgx_har,Cgy_har,Cgx,Cgy
	call ViscousAdvectionFinite(VisAdvMatrix,Nx,Ny,Lx,Ly,Cgx_har,-1.*Cgy_har,alpha_har)
	call TwoDBCApplyFirst(BCMatrix2D,VisAdvMatrix,SysMatrix,Nx,Ny)
	call invertSVD(Nx*Ny,SysMatrix,InvSysMatrix)
	call matvect(InvSysMatrix,F,Ampx_har,Nx*Ny)
			
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!   The total velocity field, the density field and   !
	!   the mean energy field							  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the first harmonic velocity field
	allocate(u_har(Nx*Ny))
	allocate(v_har(Nx*Ny))
	! the energy flux 
	allocate(E_har(Nx*Ny))
	
	! the sample time 
	t_out=0.5
	
	! let's generate the quantities
	do i=1,Nx*Ny
		! the amplitude of the velocity in z direction
		Ampz_har(i)=-1.*Ampx_har(i)*kx_har/kz_har
		! the velocity phase
		phi_har=kx_har*x_2Dreal(i,1)+kz_har*x_2Dreal(i,2)-omega_har*t_out
		! the velocity field 
		u_har(i)=Ampx_har(i)*dcos(phi_har)
		v_har(i)=Ampz_har(i)*dcos(phi_har)
		! the energy flux 
		E_har(i)=0.5*(Ampx_har(i)**2.+Ampx_har(i)**2.)*dsqrt(Cgx_har**2.+Cgy_har**2.)
	end do
	
	open(139,file='SysMatrix.dat',status='unknown')
	open(140,file='Ampx_har.dat',status='unknown')
	open(141,file='Ampz_har.dat',status='unknown')
	open(142,file='u_har.dat',status='unknown')
	open(143,file='v_har.dat',status='unknown')
	open(144,file='E_har.dat',status='unknown')
	
	do i=1,Nx*Ny
		write(139,*) SysMatrix(i,:)
		write(140,*) x_2Dreal(i,:)/lambdax,Ampx_har(i)/(A_0**2.)
		write(141,*) x_2Dreal(i,:)/lambdax,Ampz_har(i)/(A_0**2.)
		write(142,*) x_2Dreal(i,:)/lambdax,u_har(i)/(A_0**2.)
		write(143,*) x_2Dreal(i,:)/lambdax,v_har(i)/(A_0**2.)
		write(144,*) x_2Dreal(i,:)/lambdax,E_har(i)/(A_0**4.)
	end do
	
	! just to check the orhogonality relation between wave number and group velocity
	print*, Cgx_har,Cgy_har,kx_har,kz_har
	!print*, (Cgx_har*kx_har-Cgy_har*kz_har)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!   Interpolation onto some cross-section along the	  !
	!	Ray path										  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the number of points in interpolation grid
	Nx_int=30
	! let's allocate the interpolation parameters
	allocate(x_int2D(Nx_int,2))
	allocate(x_rot2D(Nx_int,2))
	allocate(IntMatrix(Nx*Ny,Nx_int))
	! the velocity field
	allocate(u_int(Nx_int))
	allocate(v_int(Nx_int))
	allocate(Ampx_int(Nx_int))
	allocate(Ampz_int(Nx_int))
	allocate(E_int(Nx_int))
	allocate(Indices(Nx*Ny))
	! let's generate evenly distributed grid
	do i=1,3
		do j=1,10
			counter=(i-1)*10+j
			x_int2D(counter,2)=-2.*lambdax+(j-1)*1.5*lambdax/9
			x_int2D(counter,1)=0.5*i*lambdax
			!print*,x_int2D(counter,:)
		end do
	end do
	
	theta_har=datan(-1.*kx_har/kz_har)
    !print*,"the angle of rotation", theta_har
	!print*, "the reference point", xr_ref/lambdax
	
	
	
	do i=1,Nx_int
		call BackRotate2D(x_rot2D(i,:),x_int2D(i,:),xr_ref,theta_har)
		!print*,x_int2D(i,:),x_rot2D(i,:)
	end do
	
	!! data refinement for interpolation, it eliminates the data far from the region which probably increases
	!! the error
	call getInterpolationPoints(Indices,counter,x_rot2D,x_2Dreal,lambdax,Nx_int,Nx*Ny)
	! the original data to be used in interpolation
	allocate(Ampx_intData(counter))
	allocate(Ampz_intData(counter))
	allocate(E_intData(counter))
	allocate(x_intData(counter,2))
	call getDataForInterpolation(Indices(1:counter),E_har,E_intData,x_2Dreal,x_intData,counter,Nx*Ny)
	
	do i=1,counter
		print*,x_intData(i,:)
	end do
	! let's perform the interpolation
	!call interpolation2D(IntMatrix,x_rot2D(:,1),x_2Dreal(:,1),x_rot2D(:,2),x_2Dreal(:,2),Nx_int,Nx*Ny)
	!call matvectnon(IntMatrix,E_har,E_int,Nx_int,Nx*Ny)
		
	!open(145,file='E_int.dat',status='unknown')
	
	!do i=1,30
	!	write(145,*) x_rot2D(i,:)/lambdax,E_int(i)/(A_0**4.)
	!end do
		
	
End Program FirstHarmonic	