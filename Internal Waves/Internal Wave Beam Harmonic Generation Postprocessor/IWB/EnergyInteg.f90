program EnergyIntegral

    implicit none

  !p1,p2 are the pointers for the nodal distance btwn the center of IWB and pycnocline in z axis
  integer i,j,k,p1,p2,nxr,nzr,Npyc,Nref,Ntri,Ntrih,Ntrih2,ti,tf,ttot
  ! the number of the data files going to be used
  integer NumDat
  ! the number of time steps to be used 
  integer Num
  ! the number of points to be used in each integration line
  integer il
  !N array to keep size of the integration limits
  integer ,allocatable::N(:)
  !kx,kz wavenumbers
  real*8 kx,kz,pi
  !r=Nmax/N1 kinamatic viscosity
  real*8 BV,N1,r,z0,arg,xnu
  !measure of the pycnocline and the pycnocline thickness
  real*8 dt,h,dummy,dummy1
  !the center of the IWB
  real*8 zcen,xcen,zlen,xlen
  !the center of the reflection 
  real*8 zcenr,xcenr 
  !the coordinates 
  real*8 , allocatable::x(:,:),z(:,:),brunt(:)
  ! the coordinates in pycnocline 
  real*8,allocatable ::xpyc(:,:),zpyc(:,:) 
  ! the coordinates in triangular region 
  real*8,allocatable ::xtri(:,:),ztri(:,:),xtrih(:,:),ztrih(:,:),xtrih2(:,:),ztrih2(:,:)
  ! the coordinates in reflection region
  real*8,allocatable ::xref(:,:),zref(:,:) 
  ! epsilon the error measure 
  real*8 eps ,E1
  ! ds the time integration difference 
  real*8 ds ,dx,dz,dl
  !the group velocities
  real*8 ,allocatable::Cgx(:,:),Cgz(:,:)
   !the group velocities in pycnocline
  real*8 ,allocatable::CgxPyc(:,:),CgzPyc(:,:)
    !the group velocities in triangular region 
  real*8 ,allocatable::Cgxtri(:,:),Cgztri(:,:),Cgxtrih(:,:),Cgztrih(:,:),Cgxtrih2(:,:),Cgztrih2(:,:)
    
   !the group velocities in reflection region
  real*8 ,allocatable::CgxRef(:,:),CgzRef(:,:)
  !the starting points 
  real*8 , allocatable :: Qx(:),Qz(:),x1(:,:),z1(:,:)
   !Nt dummy BV profile
  real*8 Nt
   !the frequency of the wave and band with for bandpass filtering and the normalization 
  real*8 w0 ,bw,nor
   ! the thickness of the integration line 
  real*8 zdim 
   !the indexes of the interpolation points 
  integer ,allocatable::Point(:,:,:),PointPyc(:,:,:),PointRef(:,:,:),PointTri(:,:,:),PointTrih(:,:,:)
   !The actual Pressure and velocities interpolated from data files in time
  real*8, allocatable::Press(:,:),u(:,:),v(:,:)
   !The actual Pressure and velocities read from data files
  real*8, allocatable::PressD(:,:),uD(:,:),vD(:,:)
  !The interpolated Pressure
  real*8, allocatable::PressI(:,:,:),PressIP(:,:,:),PressIH1(:,:,:),PressIH2(:,:,:)
  !The interpolated pressure in pycnocline
  real*8,allocatable:: PressIntPyc(:,:,:),PressIntPycP(:,:,:),PressIntPycH(:,:,:),PressIntPycH2(:,:,:)
  !The interpolated pressure in pycnocline
  real*8,allocatable:: PressIntRef(:,:,:),PressIntRefP(:,:,:),PressIntRefH(:,:,:)
   !The interpolated pressure in pycnocline
  real*8,allocatable:: PressIntTri(:,:,:),PressIntTriP(:,:,:),PressIntTriH(:,:,:),PressIntTriH2(:,:,:)
    !The interpolated pressure in pycnocline
  real*8,allocatable:: PressTri(:,:,:),PressTriP(:,:,:),PressTriH1(:,:,:),PressTriH2(:,:,:)
  ! the fourier transform of pressure
  complex , allocatable ::PressIntRefFour(:,:,:),PressIFour(:,:,:),PressIntPycFour(:,:,:),PressTriFour(:,:,:)
   !The interpolation points 
  real*8,allocatable::xi(:),zi(:),Pin(:)
   !The interpolated velocities
  real*8,allocatable::uin(:),vin(:)
   ! The array used to keep the interpolated velocities every iteration and fourier transform of it 
   real*8, allocatable::UI(:,:,:),VI(:,:,:)
   ! The array used to keep primary wave and higher order harmonics data 
   real*8, allocatable::UIP(:,:,:),VIP(:,:,:),UIH1(:,:,:),VIH1(:,:,:),UIH2(:,:,:),VIH2(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntPyc(:,:,:),VIntPyc(:,:,:)
  ! The array used to keep the interpolated pressure velocity produts in triangular zone 
   real*8, allocatable::PUTriT(:,:),PUTriZ(:,:),PUTriP(:,:),PUTriH1(:,:),PUTriH2(:,:)
   ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UTri(:,:,:),VTri(:,:,:),UTriP(:,:,:),VTriP(:,:,:)
    ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UTriH1(:,:,:),VTriH1(:,:,:),UTriH2(:,:,:),VTriH2(:,:,:)
   ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntPycP(:,:,:),VIntPycP(:,:,:),UIntPycH(:,:,:),VIntPycH(:,:,:),UIntPycH2(:,:,:),VIntPycH2(:,:,:)
! The array used to keep the interpolated velocities in triangular region
   real*8, allocatable::UIntTriP(:,:,:),VIntTriP(:,:,:),UIntTriH(:,:,:),VIntTriH(:,:,:),UIntTriH2(:,:,:),VIntTriH2(:,:,:)
!The array use to keep the fourier transform of the velociies in the triangular region
  complex ,allocatable::UTriFour(:,:,:),VTriFour(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntRef(:,:,:),VIntRef(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntRefP(:,:,:),VIntRefP(:,:,:),UIntRefH(:,:,:),VIntRefH(:,:,:)
  ! the fourier transform of the velocities
   complex , allocatable ::UIntRefFour(:,:,:),VIntRefFour(:,:,:),UIntPycFour(:,:,:),VIntPycFour(:,:,:)
  !the fourier transform of the velocities
   complex , allocatable ::UIFour(:,:,:),VIFour(:,:,:)
   !the integration path ,path1 for incident ,path2 for pycnocline ,path3 fo ref
   real*8,allocatable::path1(:),path2(:),path3(:)
   !the energy flux of the beam 
   real*8,allocatable::Energy(:),EnergyP(:),EnergyH1(:),EnergyH2(:)
   ! the dummy energy 
   real*8,allocatable::DumE(:)
   ! dummy z to be used as a dummy variable in triangular region 
   real*8 zdum
   !the energ flux of the beam in pycnocline
   real*8,allocatable:: EnergyPyc(:,:),EnergyPycP(:,:),EnergyPycH(:,:),EnergyPycH2(:,:),EnergyPycZ(:,:)
   !the energ flux of the beam in triangular region
   real*8,allocatable:: EnergyTriP(:),EnergyTri(:),EnergyTriH1(:),EnergyTriH2(:),EnergyTriZ(:)
   !the energ flux of the beam in reflection zone
   real*8,allocatable:: EnergyRef(:,:),EnergyRefP(:,:),EnergyRefH(:,:)
   
   !the array where all the interpolation matrices to be stored
   real*8,allocatable::IR(:,:,:),IRPyc(:,:,:),IRRef(:,:,:),IRTri(:,:,:),IRTrih(:,:,:)
   !the array where the number of points to be used in each interpolation stored
   integer,allocatable::Np(:)
   !the array where the number of points to be used in each interpolation stored in pycnocline
   integer,allocatable::Np1(:)
    !the array where the number of points to be used in each interpolation stored in reflection region
   integer,allocatable::Np2(:)
    !the array where the number of points to be used in each interpolation stored in triangular region 
   integer,allocatable::NpTri(:),NpTrih(:),NpTri2(:)
   ! the width of the integration path in triangular region where reflection occurs
   real*8 ,allocatable::zdimr(:),zdimrh(:),zdimrh2(:)
   
   ! the arclength path in triangular region where reflection occurs
   real*8 ,allocatable::arctri(:),arcpyc(:),arctrir(:),arcref(:),arci(:)
   ! the isoenergy line in triangular region where reflection occurs
   real*8 ,allocatable::parctri(:,:),parcpyc(:),parctrir(:),parcref(:),parci(:,:)


   ! the character where the file names kept
   CHARACTER(20),allocatable :: Pfiles(:),Vfiles(:)
   
    ! The array used to keep the interpolated pressure velocity product
   real*8, allocatable::PUInt(:,:),PVInt(:,:),PVelAbs(:,:)

     ! harmonics frequncy and wave numbers
   real*8 :: wh1,wh2,kxh1,kxh2,kzh1,kzh2

   !time step for band passfiltering 
   real*8 dt1
   ! the output time 
   integer Tout

   ! time array read from simulation data
   real*8 ,allocatable::time(:),timec(:)
   !the time steps to be used 
   integer ,allocatable::Nums(:)
   ! the non-dimensionaliztion parameters
   real*8 rho,lambdax,u0,Enon,Pnon,E1non

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      The Shifted Coordinate in Triangular Zone          !
    
   ! the shifted coordinate and Jacoian associated with it
   real*8 ,allocatable::eta(:,:),ksi(:,:),Jac(:,:),Dfx(:,:),Dfz(:,:)
    
   !  viscous dissipation associated with it
   real*8 ,allocatable::vis(:,:,:)

   !  time averaged viscous dissipation 
   real*8,allocatable::visave(:,:)
   
   ! viscous dissipation integrated
   real*8 ,allocatable::visloss(:)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rho=1E3
  
  lambdax=0.0875

  xnu=10E-8 
 
  u0=0.1*10E-3

  N1 = 2.34
  Enon=2.338*rho*(lambdax**2.)*(u0**2.)*N1*10E-7

  E1non=Enon/lambdax

  Pnon=rho*lambdax*u0*N1*10E-7

  eps=10E-3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nxr=512
  
  nzr=491


   ds=0.1

   xcen=0.2

   zcen =0.13
  

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! initial and final time step for time averaging! 
     ti=200
     tf=423
   ! and total time step tott=tf-ti+1
     ttot=tf-ti+1
   ! the filter normalization parameter
     nor=10.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !r=2.
  !dt=0.0344188927666368
   
 ! r=6
 ! dt=0.0270295880980098
  ! r=10.
  ! dt=0.025
    r=8.
    dt=0.022
 
    !0.0270295880980098
   ! r=4.

   !dt= 0.0287606935364547

   z0 = 0.2625
   
   pi=4.*datan(1.d0)
   
   kx=-2*pi/0.0875

   kz=-1.*kx

  ! zlen=0.4

  ! the number of the pressure and velocity data files
  NumDat=200

  Num=512

  ! il has to be 2**N+1 since it is used in romberg integration 
  il=33

  ! the output time 
 ! Tout=300
   
  allocate(time(NumDat))

   allocate(Nums(NumDat))

  allocate(Pfiles(NumDat))

  allocate(Vfiles(NumDat))

allocate(timec(Num))

allocate(Press(Num,nxr*nzr))

allocate(u(Num,nxr*nzr))

allocate(v(Num,nxr*nzr))

allocate(DumE(il))

allocate(x1(nxr,nzr))
allocate(z1(nxr,nzr))  

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

    ! this loop is generated to eleminate unresonable time steps
   if (time(i-1)<100.)  then 
   do k=i,NumDat-1
    ! we are fixing the time array
     time(k-1)=time(k)
     Nums(k-1)=Nums(k)
   end do
   j=j-1
   goto 50
   end if  

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



!print*,"NumDat"
!print*,j

 !  do k=1,j
    
  ! print*,k,time(k),Nums(k)

  ! end do

allocate(PressD(j,nxr*nzr))

allocate(uD(j,nxr*nzr))

allocate(vD(j,nxr*nzr))


   ! here we generate the file names for Pressure
    call pread(Pfiles,NumDat)

   ! here we generate the file names for Velocity
    call vread(Vfiles,NumDat)
   
   ! here I read all data from all files
 
    do i=1,j

    print*,Nums(i)
    call DataRead(pfiles(Nums(i)),vfiles(Nums(i)),PressD(i,:),uD(i,:),vD(i,:),nxr,nzr,x1,z1)

    end do

  
    
    ! that routine interpolates data in time with uniform time step size
 
    call TimeInter(Press,u,v,PressD,uD,vD,time,j,Num,nxr,nzr,dt1)
       
   ! the evenly disributed time stepping 

   timec=0.
    do i=1,Num-1

     timec(i+1)=timec(i)+0.01*dt1

    end do
    !print*,Press(1,35*nxr+250),PressD(1,35*nxr+250)

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Incident Beam Calculations            !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  w0=-1.*N1*kx/dsqrt(kx**2.+kz**2.)
  
  !dt1=2.*pi/(w0*8.)

  ! dt1=0.1
 
  bw=w0*dt1/2.

  !zdim=pi/dsqrt(kx**2.+kz**2.)
  
  print*,w0,dt1

  zdim=0.05

  call  bisection(arg,r)
  !Pycnocline Thickness
   h=-2.*arg*dt*0.7/2.

   print*,"pycnocline thickness",h/lambdax

allocate(Qx(il))
allocate(Qz(il))
allocate(N(il))



xlen=x1(512,1)

 call QuadPoints(Qx,Qz,kx,kz,xcen,zcen,zdim,il)

open(121,file='IntPoint.dat',status='unknown')

do i=1,il

write(121,*) Qx(i),Qz(i)

end do 


do i=1,il

 call NumPointsPyc(N(i),z0,Qz(i),491,h,kx,kz,r,dt,ds,N1,w0)

end do 

open(3,file='Npoints.dat',status='unknown')

do i=1,il

  write(3,*) N(i)
  print*,N(i)
end do
print*,"okey2"

allocate(Point(N(1),nxr*nzr,2))

allocate(path1(N(1)-1))

allocate(arci(N(il)))

allocate(parci(N(il),il))

allocate(Cgx(N(1),il))
allocate(Cgz(N(1),il))

allocate(x(N(1),il))
allocate(z(N(1),il))

allocate(Brunt(N(1)))

allocate(PressI(Num,N(1),il))

allocate(UI(Num,N(1),il))

allocate(VI(Num,N(1),il))


! fourier
allocate(PressIFour(Num,N(1),il))

allocate(UIFour(Num,N(1),il))

allocate(VIFour(Num,N(1),il))

! primary
allocate(PressIP(Num,N(1),il))

allocate(UIP(Num,N(1),il))

allocate(VIP(Num,N(1),il))

! first harmonics
allocate(PressIH1(Num,N(1),il))

allocate(UIH1(Num,N(1),il))

allocate(VIH1(Num,N(1),il))

! second harmonics
allocate(PressIH2(Num,N(1),il))

allocate(UIH2(Num,N(1),il))

allocate(VIH2(Num,N(1),il))


!allocate(PUInt(N(1),il))

!allocate(PVInt(N(1),il))

!allocate(PVelAbs(N(1),il))

allocate(Energy(N(1)))

allocate(EnergyP(N(1)))

allocate(EnergyH1(N(1)))

allocate(EnergyH2(N(1)))
! the dimension of the point array is N(1) as it is the largest 
allocate(Np(N(1)))


!print*,"okey1",N(:)

do i=1,il

 call Path(Cgx(:,i),Cgz(:,i),x(:,i),z(:,i),Qx(i),Qz(i),N1,N(i),kx,ds,r,dt,z0,w0)

end do

! the arclength starting from the center of excitation

arctrir=0.


do i=2,N(il)
 
 arci(i)=arci(i-1)+dsqrt((x(i,il)-x(i-1,il))**2.+(z(i,il)-z(i-1,il))**2.)

end do 


 
parctrir=0.

do j=1,N(il)

 do i=2,il
 
   parci(j,i)=parci(j,i-1)+dsqrt((x(j,i)-x(j,i-1))**2.+(z(j,i)-z(j,i-1))**2.)

 end do 

end do 



allocate(IR(N(1),il,nxr*nzr))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                Interpolation in IWB                                ! 
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


dz=0.5*abs(Qz(2)-Qz(1))

dx=0.5*abs(Qx(2)-Qx(1))


print*,"okey 4"
do i=1,N(il)

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,Point(i,:,:),x(i,:),z(i,:),il,nzr,nxr,Np(i),dx,dz)

 print*,Np(i)


!allocate them but it will change in every iteration 
allocate(xi(Np(i)))

allocate(zi(Np(i)))


   do j=1,Np(i)

     xi(j)=x1(Point(i,j,1),Point(i,j,2))

     zi(j)=z1(Point(i,j,1),Point(i,j,2))


   end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IR(i,:,:),x(i,:),xi,z(i,:),zi,Np(i),il)


  deallocate(xi)

  deallocate(zi)

end do

do i=1,N(il)

  do j=1,Num
 

      call interpVal(IR(i,:,:),PressI(j,i,:),Press(j,:),Point(i,:,1),Point(i,:,2),nxr,nzr,Np(i),il)
   
      call interpVal(IR(i,:,:),UI(j,i,:),u(j,:),Point(i,:,1),Point(i,:,2),nxr,nzr,Np(i),il)

      call interpVal(IR(i,:,:),VI(j,i,:),v(j,:),Point(i,:,1),Point(i,:,2),nxr,nzr,Np(i),il)

  
  end do

end do

 
! bandpass filtering in triangular region 

 !do i=1,1
  ! do j=1,il

    !Primary Wave
   !call BandPassFilter(PressIP(:,i,j),PressI(:,i,j),PressIFour(:,i,j),w0*dt1/5.,bw/5.,Num)

 
 ! end do
! end do

!open(911,file='TimeSeriesT.dat',status='unknown')

!do i=1,Num

 ! write(911,*) 2.*pi*timec(i)/w0,PressI(i,1,1)
  !print*,N(i)
!end do


!open(912,file='TimeSeriesP.dat',status='unknown')

!do i=1,Num

 ! write(912,*) 2.*pi*timec(i)/w0,PressIP(i,1,1)
  !print*,N(i)
!end do


!print*,"Frequency and time step",w0,dt1

do i=1,N(il)
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressIP(:,i,j),PressI(:,i,j),PressIFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(UIP(:,i,j),UI(:,i,j),UIFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VIP(:,i,j),VI(:,i,j),VIFour(:,i,j),w0*dt1/nor,bw/nor,Num)
 
   !  First Harmonic
  call BandPassFilter(PressIH1(:,i,j),PressI(:,i,j),PressIFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UIH1(:,i,j),UI(:,i,j),UIFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VIH1(:,i,j),VI(:,i,j),VIFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressIH2(:,i,j),PressI(:,i,j),PressIFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UIH2(:,i,j),UI(:,i,j),UIFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(VIH2(:,i,j),VI(:,i,j),VIFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  end do
 end do



! output of inlet and outlet profile at 

! k denotes output time 

k=200

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(712,file='TotPUIWB200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressI(k,j,i)*(UI(k,j,i)*Cgx(j,i)+VI(k,j,i)*Cgz(j,i))/dummy1 

  write(712,*) parci(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(713,file='PriPUIWB200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIP(k,j,i)*(UIP(k,j,i)*Cgx(j,i)+VIP(k,j,i)*Cgz(j,i))/dummy1 

  write(713,*) parci(j,i)/lambdax,dummy/E1non

  
end do


! First Harmonic Energy Flux 

open(714,file='1HPUIWB200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIH1(k,j,i)*(UIH1(k,j,i)*Cgx(j,i)+VIH1(k,j,i)*Cgz(j,i))/dummy1 

  write(714,*) parci(j,i)/lambdax,dummy/E1non

  
end do



! Second Harmonic Energy Flux 

open(715,file='2HPUIWB200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIH2(k,j,i)*(UIH2(k,j,i)*Cgx(j,i)+VIH2(k,j,i)*Cgz(j,i))/dummy1 

  write(715,*) parci(j,i)/lambdax,dummy/E1non

  
end do


! k denotes output time 

k=300

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(722,file='TotPUIWB300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressI(k,j,i)*(UI(k,j,i)*Cgx(j,i)+VI(k,j,i)*Cgz(j,i))/dummy1 

  write(722,*) parci(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(723,file='PriPUIWB300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIP(k,j,i)*(UIP(k,j,i)*Cgx(j,i)+VIP(k,j,i)*Cgz(j,i))/dummy1 

  write(723,*) parci(j,i)/lambdax,dummy/E1non

  
end do


! First Harmonic Energy Flux 

open(724,file='1HPUIWB300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIH1(k,j,i)*(UIH1(k,j,i)*Cgx(j,i)+VIH1(k,j,i)*Cgz(j,i))/dummy1 

  write(724,*) parci(j,i)/lambdax,dummy/E1non

  
end do



! Second Harmonic Energy Flux 

open(725,file='2HPUIWB300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIH2(k,j,i)*(UIH2(k,j,i)*Cgx(j,i)+VIH2(k,j,i)*Cgz(j,i))/dummy1 

  write(725,*) parci(j,i)/lambdax,dummy/E1non

  
end do


! k denotes output time 

k=400

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(732,file='TotPUIWB400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressI(k,j,i)*(UI(k,j,i)*Cgx(j,i)+VI(k,j,i)*Cgz(j,i))/dummy1 

  write(732,*) parci(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(733,file='PriPUIWB400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIP(k,j,i)*(UIP(k,j,i)*Cgx(j,i)+VIP(k,j,i)*Cgz(j,i))/dummy1 

  write(733,*) parci(j,i)/lambdax,dummy/E1non

  
end do


! First Harmonic Energy Flux 

open(734,file='1HPUIWB400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIH1(k,j,i)*(UIH1(k,j,i)*Cgx(j,i)+VIH1(k,j,i)*Cgz(j,i))/dummy1 

  write(734,*) parci(j,i)/lambdax,dummy/E1non

  
end do



! Second Harmonic Energy Flux 

open(735,file='2HPUIWB400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgx(j,i)**2.+Cgz(j,i)**2.)
dummy=PressIH2(k,j,i)*(UIH2(k,j,i)*Cgx(j,i)+VIH2(k,j,i)*Cgz(j,i))/dummy1 

  write(735,*) parci(j,i)/lambdax,dummy/E1non

  
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
!              Energy Integral for IWB                  !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! it is output time 
j=1

do i=1,N(il)


  
   ! total enegy integration along IWB 
   call PressVelocityInteg(Energy(i),PressI(j,i,:),UI(j,i,:),VI(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! Primary Energy integration along IWB
   call PressVelocityInteg(EnergyP(i),PressIP(j,i,:),UIP(j,i,:),VIP(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! First Harmonic Energy integration along IWB
   call PressVelocityInteg(EnergyH1(i),PressIH1(j,i,:),UIH1(j,i,:),VIH1(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! Second Harmonic Energy integration along IWB
   call PressVelocityInteg(EnergyH2(i),PressIH2(j,i,:),UIH2(j,i,:),VIH2(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

end do



open(300,file='EnergyIWB200.dat',status='unknown')

do i=1,N(il)

  write(300,*) arci(i)/lambdax,Energy(i)/Enon
 
end do

open(301,file='PEnergyIWB200.dat',status='unknown')

do i=1,N(il)

  write(301,*) arci(i)/lambdax,EnergyP(i)/Enon
 
end do


open(302,file='1HEnergyIWB200.dat',status='unknown')

do i=1,N(il)
  write(302,*) arci(i)/lambdax,EnergyH1(i)/Enon
 
end do


open(303,file='2HEnergyIWB200.dat',status='unknown')

do i=1,N(il)
  write(303,*) arci(i)/lambdax,EnergyH2(i)/Enon
 
end do


! it is output time 
j=300

do i=1,N(il)


  
   ! total enegy integration along IWB 
   call PressVelocityInteg(Energy(i),PressI(j,i,:),UI(j,i,:),VI(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! Primary Energy integration along IWB
   call PressVelocityInteg(EnergyP(i),PressIP(j,i,:),UIP(j,i,:),VIP(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! First Harmonic Energy integration along IWB
   call PressVelocityInteg(EnergyH1(i),PressIH1(j,i,:),UIH1(j,i,:),VIH1(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! Second Harmonic Energy integration along IWB
   call PressVelocityInteg(EnergyH2(i),PressIH2(j,i,:),UIH2(j,i,:),VIH2(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

end do


open(310,file='EnergyIWB300.dat',status='unknown')

do i=1,N(il)

  write(310,*) arci(i)/lambdax,Energy(i)/Enon
 
end do

open(311,file='PEnergyIWB300.dat',status='unknown')

do i=1,N(il)

  write(311,*) arci(i)/lambdax,EnergyP(i)/Enon
 
end do


open(312,file='1HEnergyIWB300.dat',status='unknown')

do i=1,N(il)
  write(312,*) arci(i)/lambdax,EnergyH1(i)/Enon
 
end do


open(313,file='2HEnergyIWB300.dat',status='unknown')

do i=1,N(il)
  write(313,*) arci(i)/lambdax,EnergyH2(i)/Enon
 
end do



! it is output time 
j=400

do i=1,N(il)


  
   ! total enegy integration along IWB 
   call PressVelocityInteg(Energy(i),PressI(j,i,:),UI(j,i,:),VI(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! Primary Energy integration along IWB
   call PressVelocityInteg(EnergyP(i),PressIP(j,i,:),UIP(j,i,:),VIP(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! First Harmonic Energy integration along IWB
   call PressVelocityInteg(EnergyH1(i),PressIH1(j,i,:),UIH1(j,i,:),VIH1(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

   ! Second Harmonic Energy integration along IWB
   call PressVelocityInteg(EnergyH2(i),PressIH2(j,i,:),UIH2(j,i,:),VIH2(j,i,:),il,Cgx(i,:),Cgx(i,:),zdim)

end do


open(320,file='EnergyIWB400.dat',status='unknown')

do i=1,N(il)

  write(320,*) arci(i)/lambdax,Energy(i)/Enon
 
end do

open(321,file='PEnergyIWB400.dat',status='unknown')

do i=1,N(il)

  write(321,*) arci(i)/lambdax,EnergyP(i)/Enon
 
end do


open(322,file='1HEnergyIWB400.dat',status='unknown')

do i=1,N(il)
  write(322,*) arci(i)/lambdax,EnergyH1(i)/Enon
 
end do


open(323,file='2HEnergyIWB400.dat',status='unknown')

do i=1,N(il)
  write(323,*) arci(i)/lambdax,EnergyH2(i)/Enon
 
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!        Triangular Reflection Region            !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! to track isoline in that region I use group velocities

! for primary frequency and the total energy Number of points 

! it is a critical region which requires bit more sophisticated calculation 

! first length of the integration should be determined 

! those isolines are perpendicular to group velocities therefore they are tangential to phase vectors in internal waves 

! call GroupVel(Cgx,Cgz,N1,z,z0,dt,r,kx,w0)
 
! the number of isoline along the triangular path 


!Ntri is the number of points of primary frequency in triangular region 
Ntri=N(1)-N(il)

allocate(xtri(Ntri,il))

allocate(ztri(Ntri,il))


allocate(zdimr(Ntri))

allocate(Cgxtri(Ntri,il))

allocate(Cgztri(Ntri,il))

allocate(arctri(Ntri))

allocate(parctri(Ntri,il))

  


print*," okey 3",Ntri
ds=0.0001

do i=1,Ntri

 zdimr(i)=0.

 !x(N(il)+i,1)

 zdum=z(N(il)+i,1)


 ! in this loop , we find the length of the integration path 
  do j=1,2000 
 
   zdum=zdum+ds

   call BVprofile(BV,N1,r,zdum,z0,dt)

    if(zdum>z(N(il),il)+h .or. BV<N1) then 

    goto 600

    end if 

    
     
     
      kz=dsqrt((BV**2.)*(kx**2.)/(w0**2.)-kx**2.)

     zdimr(i)=zdimr(i)+ds*dsqrt(1.+(kx/kz)**2.)
   !if(i==1) then
   !print*,zdimr(i),j
   !end if 


   end do

   600 continue 
    
    print*,zdimr(i),z(N(il),il)+h,zdum
   
    xtri(i,1)=x(N(il)+i,1)

    ztri(i,1)=z(N(il)+i,1)
  
    call GroupVel(Cgxtri(i,1),Cgztri(i,1),N1,ztri(i,1),z0,dt,r,kx,w0)

   do j=1,il-1

      call BVprofile(BV,N1,r,ztri(i,j),z0,dt)

      if( BV < w0) then 

       !print*,j,BV 

        BV=w0
       
      end if 

      kz=dsqrt((BV**2.)*(kx**2.)/(w0**2.)-kx**2.)

      xtri(i,j+1)=xtri(i,j)+(zdimr(i)/(il-1))*kx/dsqrt(kx**2.+kz**2.)

      ztri(i,j+1)=ztri(i,j)+(zdimr(i)/(il-1))*kz/dsqrt(kx**2.+kz**2.)
       

      if(zdum < ztri(i,il)) then 

         Cgxtri(i,j+1)=0.
     
         Cgztri(i,j+1)=1.

      end if 

     
      
      if(zdum > ztri(i,il)) then 
   
         call GroupVel(Cgxtri(i,j+1),Cgztri(i,j+1),N1,ztri(i,j+1),z0,dt,r,kx,w0)

      end if


   end do 


   !if (zdum < ztri(i,il)) then 

    ! zdimr(i)=zdimr(i)-eps 

     !goto 600     

   !end if 

 end do

print*,"starting point",x(N(il),1),z(N(il),1)

open(110,file='TriPath.dat',status='unknown')

do i=1,il

  write(110,*) xtri(1,i),ztri(1,i)
  !print*,N(i)
end do

open(111,file='FPointLineP.dat',status='unknown')

do i=1,Ntri

  write(111,*) xtri(i,1),ztri(i,1)
  !print*,N(i)
end do





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
! Interpolation in Triangular region for primary frequency !
!                                                          !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(Nptri(Ntri))

allocate(IRTri(Ntri,il,nxr*nzr))

allocate(PointTri(Ntri,nxr*nzr,2))

! total presurre and velocity
allocate(PressTri(Num,Ntri,il))

allocate(UTri(Num,Ntri,il))

allocate(VTri(Num,Ntri,il))

! primary presurre and velocity
allocate(PressTriP(Num,Ntri,il))

allocate(UTriP(Num,Ntri,il))

allocate(VTriP(Num,Ntri,il))

! first harmonic presurre and velocity
allocate(PressTriH1(Num,Ntri,il))

allocate(UTriH1(Num,Ntri,il))

allocate(VTriH1(Num,Ntri,il))

! second harmonic presurre and velocity
allocate(PressTriH2(Num,Ntri,il))

allocate(UTriH2(Num,Ntri,il))

allocate(VTriH2(Num,Ntri,il))

! Fourier transform of pressure and velocities
allocate(PressTriFour(Num,Ntri,il))

allocate(UTriFour(Num,Ntri,il))

allocate(VTriFour(Num,Ntri,il))


! energy arrays 

allocate(EnergyTri(Ntri))

allocate(EnergyTriP(Ntri))

allocate(EnergyTriZ(Ntri))

allocate(EnergyTriH1(Ntri))

allocate(EnergyTriH2(Ntri))

!allocate(EnergyTotTri(Ntri))

dx=0.5*abs(xtri(1,1)-xtri(1,2))

dz=0.5*abs(ztri(1,1)-ztri(1,2))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      The Shifted Coordinate in Triangular Zone          !
    

  allocate(ksi(Ntri,il))

  allocate(eta(Ntri,il))

  allocate(Jac(Ntri,il))

  allocate(Dfx(Ntri*il,Ntri*il))

  allocate(Dfz(Ntri*il,Ntri*il))

  allocate(vis(Num,Ntri,il))
  
  allocate(visave(Ntri,il))
  
  allocate(visloss(Ntri))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 Call CoorShift(xtri,ztri,eta,ksi,Ntri,il)
  
 Call Diffx(Dfx,xtri,ztri,Ntri,il)

 Call Diffx(Dfz,ztri,xtri,Ntri,il)

 Call Jacobian(Jac,Dfx,Dfz,eta,ksi,Ntri,il)

print*,"dx and dz",dx,dz
do i=1,Ntri


!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointTri(i,:,:),xtri(i,:),ztri(i,:),il,nzr,nxr,Nptri(i),dx,dz)

print*,Nptri(i)


!allocate them but it will change in every iteration 
allocate(xi(Nptri(i)))

allocate(zi(Nptri(i)))


   do j=1,Nptri(i)

     xi(j)=x1(PointTri(i,j,1),PointTri(i,j,2))

     zi(j)=z1(PointTri(i,j,1),PointTri(i,j,2))
 ! print*,xi(j),zi(j)

   end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IRTri(i,:,:),xtri(i,:),xi,ztri(i,:),zi,Nptri(i),il)


  deallocate(xi)

  deallocate(zi)

end do

print*,"okey tri"
!print*,IRTri(1,:,1)

do i=1,Ntri

  do j=1,Num
 

      call interpVal(IRTri(i,:,:),PressTri(j,i,:),Press(j,:),PointTri(i,:,1),PointTri(i,:,2),nxr,nzr,Nptri(i),il)
   
      call interpVal(IRTri(i,:,:),UTri(j,i,:),u(j,:),PointTri(i,:,1),PointTri(i,:,2),nxr,nzr,Nptri(i),il)

      call interpVal(IRTri(i,:,:),VTri(j,i,:),v(j,:),PointTri(i,:,1),PointTri(i,:,2),nxr,nzr,Nptri(i),il)

  
  end do

end do

 ! let's calculate viscous dissipation

 do j=1,Num

  call VisDis(vis(j,:,:),Dfx,Dfz,UTri(j,:,:),VTri(j,:,:),Ntri,il,xnu)

 end do
  
 ! let's take average of viscous dissipation over one wave period
 Call VisAver(visave,vis(ti:tf,:,:),ttot,Ntri,il)

 ! let's integrate it throughout the triangular region

 Call DisInteg(visloss,visave,Jac,Ntri,il,1.,1.)



!print*,"okey tri 2"
! the arclength starting from the entry of triangular region 

arctrir=0.


do i=2,Ntri
 
 arctri(i)=arctri(i-1)+dsqrt((xtri(i,1)-xtri(i-1,1))**2.+(ztri(i,1)-ztri(i-1,1))**2.)

end do 


 
parctrir=0.

do j=1,Ntri

 do i=2,il
 
   parctri(j,i)=parctri(j,i-1)+dsqrt((xtri(j,i)-xtri(j,i-1))**2.+(ztri(j,i)-ztri(j,i-1))**2.)

 end do 

end do 


  open(391,file='Visloss.dat',status='unknown')
 do i=1,Ntri
   write(391,*)  arctri(i)/lambdax,visloss(i)/Enon

 end do


print*,"okey tri 1"

! bandpass filtering in triangular region 

 do i=1,Ntri
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressTriP(:,i,j),PressTri(:,i,j),PressTriFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(UTriP(:,i,j),UTri(:,i,j),UTriFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VTriP(:,i,j),VTri(:,i,j),VTriFour(:,i,j),w0*dt1/nor,bw/nor,Num)
 
   !  First Harmonic
  call BandPassFilter(PressTriH1(:,i,j),PressTri(:,i,j),PressTriFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UTriH1(:,i,j),UTri(:,i,j),UTriFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(VTriH1(:,i,j),VTri(:,i,j),VTriFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressTriH2(:,i,j),PressTri(:,i,j),PressTriFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UTriH2(:,i,j),UTri(:,i,j),UTriFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VTriH2(:,i,j),VTri(:,i,j),VTriFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)
 
  end do
 end do

print*,"okey tri 2"



! the arclength starting from the entry of triangular region 

arctri=0.

do i=2,Ntri
 
 arctri(i)=arctri(i-1)+dsqrt((xtri(i,1)-xtri(i-1,1))**2.+(ztri(i,1)-ztri(i-1,1))**2.)

end do 


! output of inlet and outlet profile at t=200 

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(812,file='TotPU1Tri200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTri(200,1,i)*(UTri(200,1,i)*Cgxtri(1,i)+VTri(200,1,i)*Cgztri(1,i))/dummy1 

  write(812,*) parctri(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(813,file='PriPU1Tri200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriP(200,1,i)*(UTriP(200,1,i)*Cgxtri(1,i)+VTriP(200,1,i)*Cgztri(1,i))/dummy1 

  write(813,*) parctri(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(814,file='1HPU1Trir200.dat',status='unknown')

do i=1,il
dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriH1(200,1,i)*(UTriH1(200,1,i)*Cgxtri(1,i)+VTriH1(200,1,i)*Cgztri(1,i))/dummy1 

  write(814,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(815,file='2HPU1Tri200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriH2(200,1,i)*(UTriH2(200,1,i)*Cgxtri(1,i)+VTriH2(200,1,i)*Cgztri(1,i))/dummy1 

  write(815,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


j=Ntri

! Total Energy Flux at the outlet of the reflection path

open(816,file='TotPU1Tri200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTri(200,N(il),i)*(UTri(200,N(il),i)*Cgxtri(N(il),i)+VTri(200,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(816,*)  parctri(j,i)/lambdax,dummy/E1non
  
end do

! Primary Energy Flux

open(817,file='PriPU1Tri200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriP(200,N(il),i)*(UTriP(200,N(il),i)*Cgxtri(N(il),i)+VTriP(200,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(817,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(818,file='1HPU1Trir200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriH1(200,N(il),i)*(UTriH1(200,N(il),i)*Cgxtri(N(il),i)+VTriH1(200,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(818,*)  parctri(j,i)/lambdax,dummy/E1non
  
end do


! Second Harmonic Energy Flux 

open(819,file='2HPU1Tri200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriH2(200,N(il),i)*(UTriH2(200,N(il),i)*Cgxtri(N(il),i)+VTriH2(200,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(819,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do







! output of inlet and outlet profile at t=300 

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(822,file='TotPU1Tri300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTri(300,1,i)*(UTri(300,1,i)*Cgxtri(1,i)+VTri(300,1,i)*Cgztri(1,i))/dummy1 

  write(822,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(823,file='PriPU1Tri300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriP(300,1,i)*(UTriP(300,1,i)*Cgxtri(1,i)+VTriP(300,1,i)*Cgztri(1,i))/dummy1 

  write(823,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(824,file='1HPU1Tri300.dat',status='unknown')

do i=1,il
dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriH1(300,1,i)*(UTriH1(300,1,i)*Cgxtri(1,i)+VTriH1(300,1,i)*Cgztri(1,i))/dummy1 

  write(824,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(825,file='2HPU1Tri300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriH2(300,1,i)*(UTriH2(300,1,i)*Cgxtri(1,i)+VTriH2(300,1,i)*Cgztri(1,i))/dummy1 

  write(825,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


j=Ntri

! Total Energy Flux at the outlet of the reflection path

open(826,file='TotPU1Tri300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTri(300,N(il),i)*(UTri(300,N(il),i)*Cgxtri(N(il),i)+VTri(300,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(826,*)  parctri(j,i)/lambdax,dummy/E1non
  
end do

! Primary Energy Flux

open(827,file='PriPU1Trir300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriP(300,N(il),i)*(UTriP(300,N(il),i)*Cgxtri(N(il),i)+VTriP(300,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(827,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(828,file='1HPU1Tri300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriH1(300,N(il),i)*(UTriH1(300,N(il),i)*Cgxtri(N(il),i)+VTriH1(300,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(828,*)  parctri(j,i)/lambdax,dummy/E1non
  
end do


! Second Harmonic Energy Flux 

open(829,file='2HPU1Tri300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriH2(300,N(il),i)*(UTriH2(300,N(il),i)*Cgxtri(N(il),i)+VTriH2(300,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(829,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do





! output of inlet and outlet profile at t=400 

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(832,file='TotPU1Tri400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTri(400,1,i)*(UTri(400,1,i)*Cgxtri(1,i)+VTri(400,1,i)*Cgztri(1,i))/dummy1 

  write(832,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(833,file='PriPU1Tri400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriP(400,1,i)*(UTriP(400,1,i)*Cgxtri(1,i)+VTriP(400,1,i)*Cgztri(1,i))/dummy1 

  write(833,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(834,file='1HPU1Tri400.dat',status='unknown')

do i=1,il
dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriH1(400,1,i)*(UTriH1(400,1,i)*Cgxtri(1,i)+VTriH1(400,1,i)*Cgztri(1,i))/dummy1 

  write(834,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(835,file='2HPU1Tri400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(1,i)**2.+Cgztri(1,i)**2.)
dummy=PressTriH2(400,1,i)*(UTriH2(400,1,i)*Cgxtri(1,i)+VTriH2(400,1,i)*Cgztri(1,i))/dummy1 

  write(835,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


j=Ntri

! Total Energy Flux at the outlet of the reflection path

open(836,file='TotPU1Tri400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTri(400,N(il),i)*(UTri(400,N(il),i)*Cgxtri(N(il),i)+VTri(400,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(836,*)  parctri(j,i)/lambdax,dummy/E1non
  
end do

! Primary Energy Flux

open(837,file='PriPU1Tri400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriP(400,N(il),i)*(UTriP(400,N(il),i)*Cgxtri(N(il),i)+VTriP(400,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(837,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(838,file='1HPU1Tri400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriH1(400,N(il),i)*(UTriH1(400,N(il),i)*Cgxtri(N(il),i)+VTriH1(400,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(838,*)  parctri(j,i)/lambdax,dummy/E1non
  
end do


! Second Harmonic Energy Flux 

open(839,file='2HPU1Tri400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtri(N(il),i)**2.+Cgztri(N(il),i)**2.)
dummy=PressTriH2(400,N(il),i)*(UTriH2(400,N(il),i)*Cgxtri(N(il),i)+VTriH2(400,N(il),i)*Cgztri(N(il),i))/dummy1 

  write(839,*)  parctri(j,i)/lambdax,dummy/E1non

  
end do


print*,"okey tri 3"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                        !
!         Energy Integration for Triangular Region                       !
!                                                                        !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Integration of the energy  of t=200

print*,"okey tri 3"

j=200

do i=1,Ntri

 call PressVelocityInteg(EnergyTri(i),PressTri(j,i,:),UTri(j,i,:),VTri(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriP(i),PressTriP(j,i,:),UTriP(j,i,:),VTriP(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriH1(i),PressTriH1(j,i,:),UTriH1(j,i,:),VTriH1(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriH2(i),PressTriH2(j,i,:),UTriH2(j,i,:),VTriH2(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))


end do 



open(906,file='TriEnergyT200.dat',status='unknown')

do i=1,Ntri

  write(906,*) arctri(i)/lambdax,EnergyTri(i)/Enon
  
end do

open(907,file='TriEnergyP200.dat',status='unknown')

do i=1,Ntri

  write(907,*) arctri(i)/lambdax,EnergyTriP(i)/Enon
  
end do

open(908,file='TriEnergy1H200.dat',status='unknown')

do i=1,Ntri

  write(908,*) arctri(i)/lambdax,EnergyTriH1(i)/Enon
  
end do

open(909,file='TriEnergy2H200.dat',status='unknown')

do i=1,Ntri

  write(909,*) arctri(i)/lambdax,EnergyTriH2(i)/Enon
  
end do

!Integration of the energy  of t=300

j=300

do i=1,Ntri

 call PressVelocityInteg(EnergyTri(i),PressTri(j,i,:),UTri(j,i,:),VTri(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriP(i),PressTriP(j,i,:),UTriP(j,i,:),VTriP(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriH1(i),PressTriH1(j,i,:),UTriH1(j,i,:),VTriH1(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriH2(i),PressTriH2(j,i,:),UTriH2(j,i,:),VTriH2(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))


end do 



print*,"okey tri 4"

open(916,file='TriEnergyT300.dat',status='unknown')

do i=1,Ntri

  write(916,*) arctri(i)/lambdax,EnergyTri(i)/Enon
  
end do

open(917,file='TriEnergyP300.dat',status='unknown')

do i=1,Ntri

  write(917,*) arctri(i)/lambdax,EnergyTriP(i)/Enon

end do

open(918,file='TriEnergy1H300.dat',status='unknown')

do i=1,Ntri

  write(918,*) arctri(i)/lambdax,EnergyTriH1(i)/Enon

end do

open(919,file='TriEnergy2H300.dat',status='unknown')

do i=1,Ntri

  write(919,*) arctri(i)/lambdax,EnergyTriH2(i)/Enon
  
end do

!Integration of the energy  of t=400

j=400

do i=1,Ntri

 call PressVelocityInteg(EnergyTri(i),PressTri(j,i,:),UTri(j,i,:),VTri(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriP(i),PressTriP(j,i,:),UTriP(j,i,:),VTriP(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriH1(i),PressTriH1(j,i,:),UTriH1(j,i,:),VTriH1(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))

 call PressVelocityInteg(EnergyTriH2(i),PressTriH2(j,i,:),UTriH2(j,i,:),VTriH2(j,i,:),il,Cgxtri(i,:),Cgztri(i,:),0.5*zdimr(i))


end do 


print*,"okey tri 5"

open(926,file='TriEnergyT400.dat',status='unknown')

do i=1,Ntri

  write(926,*) arctri(i)/lambdax,EnergyTri(i)/Enon
  
end do

open(927,file='TriEnergyP400.dat',status='unknown')

do i=1,Ntri

  write(927,*) arctri(i)/lambdax,EnergyTriP(i)/Enon
  
end do

open(928,file='TriEnergy1H400.dat',status='unknown')

do i=1,Ntri

  write(928,*) arctri(i)/lambdax,EnergyTriH1(i)/Enon
  
end do

open(929,file='TriEnergy2H400.dat',status='unknown')

do i=1,Ntri

  write(929,*) arctri(i)/lambdax,EnergyTriH2(i)/Enon

end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          ! 
!  Averaged Energy Flux over One Wave Period  in triangular region         ! 
!                                                                          !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(PUTriT(Ntri,il))

allocate(PUTriP(Ntri,il))

allocate(PUTriZ(Ntri,il))

allocate(PUTriH1(Ntri,il))

allocate(PUTriH2(Ntri,il))



! let's generate the averaged pressure and the velocity product

do i=1,Ntri

 do j=1,il
 
  ! we are sending 107 time steps since it is equla to one wave period
  call  MeanPU(PUTriT(i,j),PressTri(ti:tf,i,j),UTri(ti:tf,i,j),VTri(ti:tf,i,j),Cgxtri(i,j),Cgztri(i,j),ttot)

  call  MeanPU(PUTriP(i,j),PressTriP(ti:tf,i,j),UTriP(ti:tf,i,j),VTriP(ti:tf,i,j),Cgxtri(i,j),Cgztri(i,j),ttot)

  call  MeanPU(PUTriH1(i,j),PressTriH1(ti:tf,i,j),UTriH1(ti:tf,i,j),VTriH1(ti:tf,i,j),Cgxtri(i,j),Cgztri(i,j),ttot)
  
  call  MeanPU(PUTriH2(i,j),PressTriH2(ti:tf,i,j),UTriH2(ti:tf,i,j),VTriH2(ti:tf,i,j),Cgxtri(i,j),Cgztri(i,j),ttot)

 end do 

end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!     Total Energy flux a entry region to Triangular Zone       !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


open(932,file='Tot1PUTriAve.dat',status='unknown')

do i=1,il

  write(932,*) parctri(1,i)/lambdax,PUTriT(1,i)/E1non
 
end do



! Total Energy flux at the outlet of the triangular zone 

open(933,file='Tot2PUTriAve.dat',status='unknown')

do i=1,il

  write(933,*) parctri(Ntri,i)/lambdax,PUTriT(Ntri,i)/E1non
 
end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for primary frequency           !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(934,file='P1PUTriAve.dat',status='unknown')

do i=1,il

  write(934,*) parctri(1,i)/lambdax,PUTriP(1,i)/E1non
 
end do



! Primary Energy flux at the outlet of the triangular zone 

open(935,file='P2PUTriAve.dat',status='unknown')

do i=1,il

  write(935,*) parctri(Ntri,i)/lambdax,PUTriP(Ntri,i)/E1non
 
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for first Harmonic              !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(936,file='1Har1PUTriAve.dat',status='unknown')

do i=1,il

  write(936,*) parctri(1,i)/lambdax,PUTriH1(1,i)/E1non
 
end do



! First harmonic Energy flux after one lambdax 

open(938,file='1Har2PUTriAve.dat',status='unknown')

do i=1,il

  write(938,*) parctri(Ntri,i)/lambdax,PUTriH1(Ntri,i)/E1non
 
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for second Harmonic             !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(939,file='2Har1PUTriAve.dat',status='unknown')

do i=1,il

  write(939,*) parctri(1,i)/lambdax,PUTriH2(1,i)/E1non
 
end do



!Second harmonic Energy flux after one lambdax 

open(940,file='2Har2PUTriAve.dat',status='unknown')

do i=1,il

  write(940,*) parctri(Ntri,i)/lambdax,PUTriH2(Ntri,i)/E1non
 
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for zeroth mode                 !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(941,file='Z1PUTriAve.dat',status='unknown')

do i=1,il

  write(941,*) parctri(1,i)/lambdax,PUTriZ(1,i)/E1non
 
end do



! First harmonic Energy flux after one lambdax 

open(942,file='Z2PUTriAve.dat',status='unknown')

do i=1,il

  write(942,*) parctri(Ntri,i)/lambdax,PUTriZ(Ntri,i)/E1non
 
end do





! Integral of the averaged energy 

do i=1,Ntri

 !Total Energy 
  call RomBerg(PUTriT(i,:),EnergyTri(i),0.5*zdimr(i),il)

 ! Primary Frequency Energy 
 call RomBerg(PUTriP(i,:),EnergyTriP(i),0.5*zdimr(i),il)

 ! First Harmonic Energy
 call RomBerg(PUTriH1(i,:),EnergyTriH1(i),0.5*zdimr(i),il)
 
 ! Second Harmonic Energy
 call RomBerg(PUTriH2(i,:),EnergyTriH2(i),0.5*zdimr(i),il)
 
 ! Zeroth Mode Energy
! call RomBerg(PUTriZ(i,:),EnergyTriZ(i),0.5*zdimr(i),il)
 

end do 



open(950,file='EnergyTriAve.dat',status='unknown')

do i=1,Ntri

  write(950,*) arctri(i)/lambdax,EnergyTri(i)/Enon
 
end do


open(951,file='PEnergyTriAve.dat',status='unknown')

do i=1,Ntri

  write(951,*) arctri(i)/lambdax,EnergyTriP(i)/Enon
 
end do


open(952,file='1HEnergyTriAve.dat',status='unknown')

do i=1,Ntri

  write(952,*) arctri(i)/lambdax,EnergyTriH1(i)/Enon
 
end do


open(953,file='2HEnergyTriAve.dat',status='unknown')

do i=1,Ntri

  write(953,*) arctri(i)/lambdax,EnergyTriH2(i)/Enon
 
end do

!open(822,file='ZEnergyPycAve.dat',status='unknown')

!do i=1,Npyc

 ! write(822,*) arcpyc(i)/lambdax,EnergyPycZ(i)/Enon
 
!end do

end program EnergyIntegral
