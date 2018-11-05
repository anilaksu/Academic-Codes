program FieldFilter

	implicit none
	! The nonlinear higher harmonics and mean flow generation is analyzed by postprocessing 
	! the pressure and the velocity data obtained form the main solver
	! copyright (c) Anil A. Aksu, 2018 
	! The particular program band-pass filters the pressure and the velocity data around the zeroth mode
	! the primary frequency wave and the second harmonic frequency and calculates the mean energy flux and energy density
	integer i,j,k,Nx,Nz
	! ti initial time step, tf final time step for time averaging and ttot is the total time steps
	integer ti,tf,ttot,t1
	! the number of the data files going to be used
	integer NumDat
	! the number of time steps to be used 
	integer Num
	!N array to keep size of the integration limits
	integer ,allocatable::N(:)
	!kx,kz wavenumbers
	real*8 kx,kz,pi
	!the center of the IWB
	real*8 zcen,xcen,zlen,xlen
	!the coordinates 
	real*8 , allocatable::x1(:,:),z1(:,:)
	!Nt dummy BV profile
	real*8 Nt,N1
	!the frequency of the wave and band with for bandpass filtering
	real*8 w0 ,bw,nor
	!The actual Pressure and velocities interpolated from data files in time
	real, allocatable::Press(:,:),u(:,:),v(:,:)
	!The actual Pressure and velocities read from data files
	real, allocatable::PressD(:,:),uD(:,:),vD(:,:)
	!The interpolated Pressure
	real, allocatable::PressP(:,:),PressH(:,:),PressH2(:,:),PressZ(:,:)
	! the fourier transform of pressure
	complex , allocatable ::PressFour(:,:)
	!The interpolation points 
	real ,allocatable::xi(:),zi(:),Pin(:)
	!The interpolated velocities
	real ,allocatable::uin(:),vin(:)
	! The array used to keep the interpolated velocities every iteration and fourier transform of it 
	real, allocatable::UInt(:,:,:),VInt(:,:,:)
	! The array used to keep primary wave and higher order harmonics data 
	real, allocatable::UP(:,:),VP(:,:),UH(:,:),VH(:,:),UH2(:,:),VH2(:,:),UZ(:,:),VZ(:,:)
	!the fourier transform of the velocities
	complex, allocatable ::UFour(:,:),VFour(:,:)
	! the character where the file names kept
	CHARACTER(20),allocatable :: Pfiles(:),Vfiles(:)
	! The array used to keep the interpolated pressure velocity product
	real, allocatable::PU(:,:),PV(:,:),PVelAbs(:,:)
	! The array used to keep the interpolated pressure velocity product of primary wave
	real, allocatable::PUP(:,:),PVP(:,:),PVelAbsP(:,:)
	! The array used to keep the interpolated pressure velocity product of primary wave
	real, allocatable::PUH(:,:),PVH(:,:),PVelAbsH(:,:)
	! The array used to keep Mean Pressure Velocity Product
	real, allocatable::PUMean(:,:),PVMean(:,:),PUPMean(:,:),PVPMean(:,:),PUZMean(:,:),PVZMean(:,:)
	! the time averaged energy and the velocity fields
	real*8, allocatable::U2Mean(:,:),V2Mean(:,:),EMean(:,:)
	! The array used to keep Mean Pressure Velocity Product
	real, allocatable::PUH1Mean(:,:),PVH1Mean(:,:),PUH2Mean(:,:),PVH2Mean(:,:)
	!time step for band passfiltering 
	real*8 dt1
	! time array read from simulation data
	real*8 ,allocatable::time(:)
	!the time steps to be used 
	integer ,allocatable::Nums(:)
	! the non-dimensionaliztion parameters
	real*8 rho,lambdax,u0,Enon,Pnon,E1non

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	rho=1E3
	lambdax=0.0875*11.63
	u0=(0.1**2.)*0.00014/lambdax
	N1 = 3.9
	Enon=2.338*rho*(lambdax**2.)*(u0**2.)*N1*10E-7
	E1non=Enon/lambdax
	Pnon=rho*lambdax*u0*N1*10E-7

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 		the number of grid points		  		!
	Nx=768
	Nz=533
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initial and final time step for time averaging! 
	ti=200
	tf=316
	! and total time step tott=tf-ti+1
	ttot=tf-ti+1
	nor=10.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! pi number
	pi=4.*datan(1.d0)
	! the wavenumber vector
	kx=-2*pi/0.0875
	kz=-1.*kx
	! the number of the pressure and velocity data files
	NumDat=200
	Num=512

	allocate(time(NumDat))
	allocate(Nums(NumDat))
	allocate(Pfiles(NumDat))
	allocate(Vfiles(NumDat))
	! the allocation of grid points
	allocate(x1(Nx,Nz))
	allocate(z1(Nx,Nz))  

	! t1 is used to specify the time to output the data 
	t1=150
	! read time data from postprocessor

	open(121,file='time.dat',status='old')

	do i=1,NumDat

		Nums(i)=i

	end do

	do i=1,NumDat

		read(121,*) j, time(i)

	end do

	! Now I observed that at time steps it is dumping out non-reasonable times I should filter them out 
	! j is the counter of the non-reasonable time steps 
	j=0
	do i=2,NumDat
		j=j+1
		if (time(i-1)>time(i))  then 
			do k=i,NumDat-1
				! we are fixing the time array
				time(k)=time(k+1)
				Nums(k)=Nums(k+1)
			end do

			j=j-1
			goto 50
		end if 

		!Nums(j)=i-1

		50 continue
	end do

	allocate(Press(j,Nx*Nz))
	allocate(u(j,Nx*Nz))
	allocate(v(j,Nx*Nz))


	!print*,"NumDat"
	!print*,j

	!  do k=1,j
	! print*,k,time(k),Nums(k)
	! end do

	allocate(PressD(Num,Nx*Nz))
	allocate(uD(Num,Nx*Nz))
	allocate(vD(Num,Nx*Nz))

	! primary , harmonic and frequency domain pressure 
	! the third  harmonic
	allocate(PressH2(Num,Nx*Nz))
	! the second harmonic
	allocate(PressH(Num,Nx*Nz))
	! the primary frequency
	allocate(PressP(Num,Nx*Nz))
	! mean pressure
	allocate(PressZ(Num,Nx*Nz))
	! Fourier transform of Pressure
	allocate(PressFour(Num,Nx*Nz))

	! primary , harmonic and frequency domain u velocity 
	! the third  harmonic
	allocate(UH2(Num,Nx*Nz))
	! the second  harmonic
	allocate(UH(Num,Nx*Nz))
	! the first harmonic
	allocate(UP(Num,Nx*Nz))
	! mean u velocity
	allocate(UZ(Num,Nx*Nz))
	! Fourier transform of u velocity
	allocate(UFour(Num,Nx*Nz))

	! primary , harmonic and frequency domain u velocity 
	allocate(VH2(Num,Nx*Nz))
	allocate(VH(Num,Nx*Nz))
	allocate(VP(Num,Nx*Nz))
	allocate(VZ(Num,Nx*Nz))
	allocate(VFour(Num,Nx*Nz))


	! pressure velocity products and their L2 norms

	allocate(PU(Num,Nx*Nz))
	allocate(PV(Num,Nx*Nz))
	allocate(PVelAbs(Num,Nx*Nz))

	! pressure velocity products and their L2 norms of primary wave
	allocate(PUP(Num,Nx*Nz))
	allocate(PVP(Num,Nx*Nz))
	allocate(PVelAbsP(Num,Nx*Nz))

	! pressure velocity products and their L2 norms of primary wave
	allocate(PUH(Num,Nx*Nz))
	allocate(PVH(Num,Nx*Nz))
	allocate(PVelAbsH(Num,Nx*Nz))


	!Mean Quantities
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!     				Mean Quantities  	              !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	allocate(PUMean(Nx,Nz))
	allocate(PVMean(Nx,Nz))  

	allocate(PUZMean(Nx,Nz))
	allocate(PVZMean(Nx,Nz))  

	allocate(PUPMean(Nx,Nz))
	allocate(PVPMean(Nx,Nz))  

	allocate(PUH1Mean(Nx,Nz))
	allocate(PVH1Mean(Nx,Nz))  

	allocate(PUH2Mean(Nx,Nz))
	allocate(PVH2Mean(Nx,Nz))  

	allocate(U2Mean(Nx,Nz))
	allocate(V2Mean(Nx,Nz))
	allocate(EMean(Nx,Nz))


	! here we generate the file names for Pressure
	call pread(Pfiles,NumDat)
	! here we generate the file names for Velocity
	call vread(Vfiles,NumDat)

	! here I read all data from all files

	do i=1,j

		print*,Nums(i)
		call DataRead(pfiles(Nums(i)),vfiles(Nums(i)),Press(i,:),u(i,:),v(i,:),Nx,Nz,x1,z1)

	end do
	print*,v(1,35*Nx+250)


	! that routine interpolates data in time with uniform time step size
	call TimeInter(PressD,uD,vD,Press,u,v,time,j,Num,Nx,Nz,dt1)
	print*,vD(1,35*Nx+250)

	! the dispersion relation
	w0=-1.*N1*kx/dsqrt(kx**2.+kz**2.)

	! dt1=2.*pi/(w0*8.)
	! dt1=0.1
	! band-width where the filter function is centered around
	bw=w0*dt1/2.

	print*,"w,bw,dt1"
	print*,w0,bw,dt1
	print*,"okey0"

	call BandPassFilter(VP(:,35*Nx+250),vD(:,35*Nx+250),VFour(:,35*Nx+250),w0*dt1/nor,bw/nor,Num)
	
	open(20,file='Unfiltered.dat',status='unknown')
	open(21,file='filtered.dat',status='unknown')
	open(22,file='fourier.dat',status='unknown')

	do i=1,Num
		write(21,*) i,vP(i,35*Nx+250)  
		write(22,*) i,abs(VFour(i,35*Nx+250))   
		write(20,*) i,vD(i,35*Nx+250)  
	end do

	do i=1,Nx*Nz

		!zeroth mode Wave
		!  call BandPassFilter(PressZ(:,i),PressD(:,i),PressFour(:,i),0.,bw/nor,Num)

		!   call BandPassFilter(UZ(:,i),uD(:,i),UFour(:,i),0.,bw/nor,Num)

		!   call BandPassFilter(VZ(:,i),vD(:,i),VFour(:,i),0.,bw/nor,Num)

		!Primary Wave
		!  call BandPassFilter(PressP(:,i),PressD(:,i),PressFour(:,i),w0*dt1/nor,bw/nor,Num)

		call BandPassFilter(UP(:,i),uD(:,i),UFour(:,i),w0*dt1/nor,bw/nor,Num)

		call BandPassFilter(VP(:,i),vD(:,i),VFour(:,i),w0*dt1/nor,bw/nor,Num)

		!  First Harmonic
		!  call BandPassFilter(PressH(:,i),PressD(:,i),PressFour(:,i),2.*w0*dt1/nor,bw/nor,Num)

		!  call BandPassFilter(UH(:,i),uD(:,i),UFour(:,i),2.*w0*dt1/nor,bw/nor,Num)

		!  call BandPassFilter(VH(:,i),vD(:,i),VFour(:,i),2.*w0*dt1/nor,bw/nor,Num)

		! Second Harmonic
		!   call BandPassFilter(PressH2(:,i),PressD(:,i),PressFour(:,i),3.*w0*dt1/5.,bw/5.,Num)

		!   call BandPassFilter(UH2(:,i),uD(:,i),UFour(:,i),3.*w0*dt1/5.,bw/5.,Num)

		!   call BandPassFilter(VH2(:,i),vD(:,i),VFour(:,i),3.*w0*dt1/5.,bw/5.,Num)


	end do

	print*,"okey"

	!Pressure output for entire field total primary and harmonic

	open(19,file='FourPressuret100.dat',status='unknown')
	open(20,file='TPressuret100.dat',status='unknown')
	open(21,file='PPressuret100.dat',status='unknown') 
	open(22,file='HPressuret100.dat',status='unknown') 

	! let's output the frequency data of the pressure
	do i=1,Num

	write(19,*) abs(PressFour(i,119579))

	end do

	do j=1,Nz
		do i=1,Nx

			!  write(19,*) x1(i,j)/lambdax,z1(i,j)/lambdax,abs(PressFour(t1,(j-1)*Nx+i))
			write(20,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressD(t1,(j-1)*Nx+i)/Pnon
			write(21,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressP(t1,(j-1)*Nx+i)/Pnon
			write(22,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressH(t1,(j-1)*Nx+i)/Pnon

		end do 
	end do


	!U velocity output for entire field total primary and harmonic
	open(23,file='TUt100.dat',status='unknown')
	open(24,file='PUt100.dat',status='unknown') 
	open(25,file='HUt100.dat',status='unknown') 

	do j=1,Nz
		do i=1,Nx

			write(23,*) x1(i,j)/lambdax,z1(i,j)/lambdax,uD(t1,(j-1)*Nx+i)/u0
			write(24,*) x1(i,j)/lambdax,z1(i,j)/lambdax,UP(t1,(j-1)*Nx+i)/u0
			write(25,*) x1(i,j)/lambdax,z1(i,j)/lambdax,UH(t1,(j-1)*Nx+i)/u0

		end do 
	end do

	!V velocity output for entire field total primary and harmonic
	open(26,file='TVt100.dat',status='unknown')
	open(27,file='PVt100.dat',status='unknown') 
	open(28,file='HVt100.dat',status='unknown') 

	do j=1,Nz
		do i=1,Nx

			write(26,*) x1(i,j)/lambdax,z1(i,j)/lambdax,vD(t1,(j-1)*Nx+i)/u0
			write(27,*) x1(i,j)/lambdax,z1(i,j)/lambdax,VP(t1,(j-1)*Nx+i)/u0
			write(28,*) x1(i,j)/lambdax,z1(i,j)/lambdax,VH(t1,(j-1)*Nx+i)/u0

		end do 
	end do
	
	!PU product for total primary and harmonic
	open(29,file='TPUt100.dat',status='unknown')
	open(30,file='PPUt100.dat',status='unknown') 
	open(31,file='HPUt100.dat',status='unknown') 

	do j=1,Nz
		do i=1,Nx

		!  write(29,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressD(t1,(j-1)*Nx+i)*uD(t1,(j-1)*Nx+i)/Enon
		!  write(30,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressP(t1,(j-1)*Nx+i)*UP(t1,(j-1)*Nx+i)/Enon
		!  write(31,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressH(t1,(j-1)*Nx+i)*UH(t1,(j-1)*Nx+i)/Enon

		end do 
	end do

	!PV product for total primary and harmonic
	open(32,file='TPVt100.dat',status='unknown')
	open(33,file='PPVt100.dat',status='unknown') 
	open(34,file='HPVt100.dat',status='unknown') 

	do j=1,Nz
		do i=1,Nx

		!  write(32,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressD(t1,(j-1)*Nx+i)*vD(t1,(j-1)*Nx+i)/Enon
		!  write(33,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressP(t1,(j-1)*Nx+i)*VP(t1,(j-1)*Nx+i)/Enon
		!  write(34,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressH(t1,(j-1)*Nx+i)*VH(t1,(j-1)*Nx+i)/Enon

		end do 
	end do

	print*,"okey1"

	do j=1,Nz
		do i=1,Nx

			! we are sending 107 time steps since it is equal to one wave period

			k=(j-1)*Nx+i

			! PU flux
			! call  MeanPU(PUZMean(i,j),PressZ(ti:tf,k),UZ(ti:tf,k),UZ(ti:tf,k),1.,0.,ttot)
			! call  MeanPU(PUMean(i,j),PressD(ti:tf,k),uD(ti:tf,k),vD(ti:tf,k),1.,0.,ttot)
			!  call  MeanPU(PUPMean(i,j),PressP(ti:tf,k),UP(ti:tf,k),VP(ti:tf,k),1.,0.,ttot)
			! call  MeanPU(PUH1Mean(i,j),PressH(ti:tf,k),UH(ti:tf,k),VH(ti:tf,k),1.,0.,ttot)
			! call  MeanPU(PUH2Mean(i,j),PressH2(ti:tf,k),UH2(ti:tf,k),VH2(ti:tf,k),1.,0.,ttot)

			!PV flux
			! call  MeanPU(PVZMean(i,j),PressZ(ti:tf,k),UZ(ti:tf,k),UZ(ti:tf,k),0.,1.,ttot)
			! call  MeanPU(PVMean(i,j),PressD(ti:tf,k),uD(ti:tf,k),vD(ti:tf,k),0.,1.,ttot)
			!  call  MeanPU(PVPMean(i,j),PressP(ti:tf,k),UP(ti:tf,k),VP(ti:tf,k),0.,1.,ttot)
			! call  MeanPU(PVH1Mean(i,j),PressH(ti:tf,k),UH(ti:tf,k),VH(ti:tf,k),0.,1.,ttot)
			! call  MeanPU(PVH2Mean(i,j),PressH2(ti:tf,k),UH2(ti:tf,k),VH2(ti:tf,k),0.,1.,ttot)

			!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!        Mean Kinetic Enegy is Calculated            !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			call  MeanPU(U2Mean(i,j),UP(ti:tf,k),UP(ti:tf,k),VP(ti:tf,k),1.,0.,ttot)
			call  MeanPU(V2Mean(i,j),VP(ti:tf,k),UP(ti:tf,k),VP(ti:tf,k),0.,1.,ttot)

			EMean(i,j)=0.5*(U2Mean(i,j)+V2Mean(i,j))

		end do 
	end do
	
	! it is a flag to see the program proceed to this step 
	print*,"okey2"

	open(118,file='EMean.dat',status='unknown') 
	open(99,file='ZPUMean.dat',status='unknown')
	open(100,file='TPUMean.dat',status='unknown')
	open(101,file='PPUMean.dat',status='unknown') 
	open(102,file='1HPUMean.dat',status='unknown') 
	open(103,file='2HPUMean.dat',status='unknown') 
	open(108,file='ZPVMean.dat',status='unknown')        
	open(104,file='TPVMean.dat',status='unknown')
	open(105,file='PPVMean.dat',status='unknown') 
	open(106,file='1HPVMean.dat',status='unknown') 
	open(107,file='2HPVMean.dat',status='unknown') 


	do j=1,Nz
		do i=1,Nx

			! PU flux
			write(99,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PUZMean(i,j)/E1non
			write(100,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PUMean(i,j)/E1non
			write(101,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PUPMean(i,j)/E1non
			write(102,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PUH1Mean(i,j)/E1non
			write(103,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PUH2Mean(i,j)/E1non

			!PV flux
			write(108,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PVZMean(i,j)/E1non
			write(104,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PVMean(i,j)/E1non
			write(105,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PVPMean(i,j)/E1non
			write(106,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PVH1Mean(i,j)/E1non
			write(107,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PVH2Mean(i,j)/E1non

			! Mean Kinetic Energy 
			write(118,*) x1(i,j)/lambdax,z1(i,j)/lambdax,EMean(i,j)

		end do 
	end do

end program FieldFilter

