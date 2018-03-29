
Program PlasmaLaserBeam
	! this program is written by Anil Aksu to solve 3D transient Non-linear Scrödinger Equation in cartesian coordinates
	! Year: 2018
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!   This is general non-linear schrödinger equation   !
	!   solver for plasma laser wave beam. It solves that !
	!   problem for wide range of parameter.     		  !
	!													  !
	!   With great power comes great responsibility !!!   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
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
	! the number of subgrids 
	integer Nsub
	! the length of domain in each direction
	real*8 :: Lx,Ly,Lz
	! the number of points in each subgrid in x and y direction
	integer Nx,Ny,Nz
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!         The time integration parameters             !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the time integration parameter and the step size
	real*8 theta, dt,time_tot
	! the number of time steps
	integer Ntime
	! the explicit part of integration
	real*8, allocatable:: f_time(:)
	! the time array 
	real*8, allocatable:: time_s(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!   Electric Field and its time spatial derivatives   !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Electric field and its time derivative
	complex*8, allocatable::E_field(:,:),E_t(:,:)
	! the spatial derivatives of electric field
	complex*8, allocatable::E_x(:,:),E_z(:,:)
	! the differentiation matrix in z directiom
	real*8, allocatable:: Diffz(:,:)
	! the system matrix
	complex*8, allocatable:: SysMatrix(:,:)
	! the eigenvalues and eigenvector of Diffz
	!real*8, allocatable:: EValR(:),EValI(:),EVec(:,:),VL(:,:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!                 Plasma Parameters			          !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! linear and non-linear indices of refraction
	real*8 n0,n2
	! the  wavenumber and wavelength in z direction 
	real*8 kz,lambdaz
	! GVD parameter and the critical power of self-focusing
	real*8 beta,Pcr
	! the initial value of the amplitude
	real*8 A_0,sigma

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!        2-D Grid and related parameters			  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! the grid points
	real*8, allocatable:: x_grid2D(:,:),wx(:),wy(:)
	! Test Function
	real*8, allocatable:: F(:),dF(:),FB(:),LegenData(:)
	! the right hand side of the equation
	complex*8, allocatable:: RHS(:)
	
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
	!      3D Transient Schrödinger Equation 		      !
	!	   Solution on Fourier-Fourier-GLL Grid     	  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! the wave number in z direction 
	lambdaz=1.2
	! the wavenumber in z direction
	kz=2.*pi/lambdaz
	! the number of subgrids
	Nsub=5 
	! the  number of time steps
	Ntime=1000
	! domian properties: number of points in domain and length of domain
	Nx=20
	Ny=20
	Nz=20
	! the subdomain length in x direction and y direction
	Lx=20.*lambdaz
	Ly=20.*lambdaz
	Lz=50.*lambdaz
	! the time step
	dt=0.01
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
		!print*,x_2Dreal(i,:)
	end do
	! differentiation matrix in z direction
	allocate(Diffz(Nz,Nz))
	! the governing equation matrix
	allocate(SysMatrix(Nz,Nz))
	! the RHS value matrix vector operations
	allocate(RHS(Nz))
	allocate(E_field(Ntime,Nz))
	allocate(E_t(Ntime,Nz))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!      1D Transient Linear Schrödinger Equation       !
	!			Solution on GLL Grid     		    	  !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call FirstDiff(Diffz,Nz)
	! the differentiation matrix defined above is in mother interval -1 to 1 
	! the domain length in z direction Lz, therefore the differentiation matrix should be multiplied with 2/Lz 
	Diffz=2.*Diffz/Lz
	! Governing equation matrix for the second time derivative in Schrödinger equation
	SysMatrix=0.
	do i=1,Nz
		! the first row is set to 1.,0,0,.. to satisfy the initial condition
		if (i==1) then
			SysMatrix(i,i)=0.
		else
			SysMatrix(i,:)=2.*CMPLX(0.,1.)*kz*Diffz(i,:)
			! this will modified after Fourier-Fourier discretization in x-y directions
			SysMatrix(i,i)=SysMatrix(i,i)+10.
		end if
	end do
	! the value of the envelope
	A_0=3.
	E_field(1,:)=0.
	E_t(1,:)=0.
	E_field(1,1)=A_0
	! the time integration
	! the time integration will be replaced with Adam-Bashford time integration
	do i=1,(Ntime-1)
		! the second time derivative of the electric field
		call cmatvect(SysMatrix,E_field(1,:),RHS,Nz)
		! the computation of the first time derivative of the electric field
		E_t(i+1,:)=E_t(i,:)+dt*RHS
		! the electric field itself
		E_field(i+1,:)=E_field(i,:)+dt*E_t(i,:)
	end do
	
	!print*,SysMatrix(2,:)
	
	open(100,file='E_tot.dat',status='unknown')
	
	do i=1,Nz
		write(100,*) x_2Dreal(i,1)/lambdaz,abs(E_field(Ntime,i)/A_0)
		!print*, x_2Dreal(i,1)/lambdaz,E_t(Ntime,i)/A_0
	end do
			
End Program PlasmaLaserBeam