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
  !r=Nmax/N1
  real*8 BV,N1,r,z0,arg
  !measure of the pycnocline and the pycnocline thickness
  real*8 dt,h,dummy,dummy1,dummy2
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
  real*8, allocatable::PressInt(:,:,:),PressIntP(:,:,:),PressIntH(:,:,:)
  !The interpolated pressure in pycnocline
  real*8,allocatable:: PressPyc(:,:,:),PressPycP(:,:,:),PressPycH1(:,:,:),PressPycH2(:,:,:),PressPycM(:,:),PressPycZ(:,:,:)
  !The interpolated pressure in pycnocline
  real*8,allocatable:: PressIntRef(:,:,:),PressIntRefP(:,:,:),PressIntRefH(:,:,:)
    !The interpolated pressure in pycnocline
  real*8,allocatable:: PressTri(:,:,:),PressTriP(:,:,:),PressTriH1(:,:,:),PressTriH2(:,:,:)
  ! the fourier transform of pressure
  complex , allocatable ::PressIntRefFour(:,:,:),PressIntFour(:,:,:),PressPycFour(:,:,:),PressTriFour(:,:,:)
   !The interpolation points 
  real*8,allocatable::xi(:),zi(:),Pin(:)
   !The interpolated velocities
  real*8,allocatable::uin(:),vin(:)
   ! The array used to keep the interpolated velocities every iteration and fourier transform of it 
   real*8, allocatable::UInt(:,:,:),VInt(:,:,:)
   ! The array used to keep primary wave and higher order harmonics data 
   real*8, allocatable::UIntP(:,:,:),VIntP(:,:,:),UIntH(:,:,:),VIntH(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline and the mean properties
   real*8, allocatable::UPyc(:,:,:),VPyc(:,:,:),UPycM(:,:),VPycM(:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntTri(:,:,:),VIntTri(:,:,:)
   ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UTri(:,:,:),VTri(:,:,:),UTriP(:,:,:),VTriP(:,:,:)
    ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UTriH1(:,:,:),VTriH1(:,:,:),UTriH2(:,:,:),VTriH2(:,:,:)
   ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UPycP(:,:,:),VPycP(:,:,:),UPycH1(:,:,:),VPycH1(:,:,:),UPycH2(:,:,:),VPycH2(:,:,:),UPycZ(:,:,:),VPycZ(:,:,:)
!The array use to keep the fourier transform of the velociies in the triangular region
  real*8 ,allocatable::UTriFour(:,:,:),VTriFour(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntRef(:,:,:),VIntRef(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntRefP(:,:,:),VIntRefP(:,:,:),UIntRefH(:,:,:),VIntRefH(:,:,:)
  ! the fourier transform of the velocities
   complex , allocatable ::UIntRefFour(:,:,:),VIntRefFour(:,:,:),UPycFour(:,:,:),VPycFour(:,:,:)
  !the fourier transform of the velocities
   complex , allocatable ::UIntFour(:,:,:),VIntFour(:,:,:)
   !the integration path ,path1 for incident ,path2 for pycnocline ,path3 fo ref
   real*8,allocatable::path1(:),path2(:),path3(:)
   !the energy flux of the beam 
   real*8,allocatable::Energy(:),EnergyP(:),EnergyH(:)
   ! the dummy energy 
   real*8,allocatable::DumE(:)
   ! dummy z to be used as a dummy variable in triangular region 
   real*8 zdum
   !the energ flux of the beam in pycnocline
   real*8,allocatable:: EnergyPyc(:),EnergyPycP(:),EnergyPycH1(:),EnergyPycH2(:),EnergyPycM(:),EnergyPycZ(:)
   ! the averaged pressure velocity product inpycnocline 
   real*8,allocatable::PUT(:,:),PUZ(:,:),PUP(:,:),PUH1(:,:),PUH2(:,:)
    ! the averaged pressure velocity product radiated from pycnocline 
   real*8,allocatable::PURadT(:),PURadZ(:),PURadP(:),PURadH1(:),PURadH2(:)
    !the energ flux of the beam in pycnocline
   real*8,allocatable:: EnergyRad(:),EnergyRadP(:),EnergyRadH1(:),EnergyRadH2(:),EnergyRadZ(:)
     !the energ flux of the beam in pycnocline
   real*8,allocatable:: DEnergyRad(:),DEnergyRadP(:),DEnergyRadH1(:),DEnergyRadH2(:)
   !the energ flux of the beam in triangular region
   real*8,allocatable:: EnergyTriP(:),EnergyTri(:),EnergyTriH1(:),EnergyTriH2(:)
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
   real*8 ,allocatable::arctri(:),arcpyc(:),arctrir(:),arcref(:)
   ! the isoenergy line in triangular region where reflection occurs
   real*8 ,allocatable::parctri(:),parcpyc(:),parctrir(:),parcref(:)

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rho=1E3
  
  lambdax=0.0875

  u0=0.1*10E-3

  N1 = 3.9

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

  !r=2.
  !dt=0.0344188927665368

  !  r=4.  
  !  dt=   0.0287606935364547
 
  ! r=6.
  ! dt=0.0270295880980098
   r=8.
   dt=0.022
   !r=10.
   !dt=0.025

   z0 = 0.2625
   
   pi=4.*datan(1.d0)
   
   kx=-2*pi/0.0875

   kz=-1.*kx

   zlen=0.4

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! initial and final time step for time averaging! 
     ti=200
     tf=423
   ! and total time step tott=tf-ti+1
     ttot=tf-ti+1

   ! filter normalization parameter
    nor=10.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       
    ! the evenly distributed time 

    timec=0.

    do i=1,Num-1

      timec(i+1)=timec(i) + 0.01*dt1
    
    end do 

 
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

   print*,"pycnocline thickness",h

allocate(Qx(il))
allocate(Qz(il))
allocate(N(il))



xlen=x1(384,1)

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


allocate(Cgx(N(1),il))
allocate(Cgz(N(1),il))

allocate(x(N(1),il))
allocate(z(N(1),il))

allocate(Brunt(N(1)))

allocate(PressInt(Num,N(1),il))

allocate(UInt(Num,N(1),il))

allocate(VInt(Num,N(1),il))

allocate(PressIntFour(Num,N(1),il))

allocate(UIntFour(Num,N(1),il))

allocate(VIntFour(Num,N(1),il))

allocate(PressIntP(Num,N(1),il))

allocate(UIntP(Num,N(1),il))

allocate(VIntP(Num,N(1),il))

allocate(PressIntH(Num,N(1),il))

allocate(UIntH(Num,N(1),il))

allocate(VIntH(Num,N(1),il))

allocate(PUInt(N(1),il))

allocate(PVInt(N(1),il))

allocate(PVelAbs(N(1),il))

allocate(Energy(N(1)))

allocate(EnergyP(N(1)))

allocate(EnergyH(N(1)))
! the dimension of the point array is N(1) as it is the largest 
allocate(Np(N(1)))


!print*,"okey1",N(:)

do i=1,il

 call Path(Cgx(:,i),Cgz(:,i),x(:,i),z(:,i),Qx(i),Qz(i),N1,N(i),kx,ds,r,dt,z0,w0)

end do

allocate(IR(N(1),il,nxr*nzr))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Pycnocline Calculations               !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!this do loop is used to calculate the number of points in pycnocline

ds=0.0025


do i=1,1000

if (ds*i > xlen-x(N(1),1)) then

Npyc=i-1

goto 500
end if





end do

500 continue

print*,"Number of points in pycnocline" ,Npyc,xlen-x(N(il),il)
!Npyc=3

! to calculate the quadrature points in pycnocline 

! we start from x(N(5),5) to reduce the computational effort 

! And also our beam enters to pycnocline at that point 
 call QuadPoints(Qx,Qz,0.,1.,x(N(1),1),z(N(1),1)+h/2,h/2,il) 

!print*,Qx
!print*,Qz
! now lets find the interpolation points in pycnocline


allocate(xpyc(Npyc,il))
allocate(zpyc(Npyc,il))

allocate(parcpyc(il))
allocate(arcpyc(Npyc))


! the interpolated data in pycnocline
allocate(PressPyc(Num,Npyc,il))

allocate(UPyc(Num,Npyc,il))

allocate(VPyc(Num,Npyc,il))

! fourier transformed data
allocate(PressPycFour(Num,Npyc,il))

allocate(UPycFour(Num,Npyc,il))

allocate(VPycFour(Num,Npyc,il))

! primary frequency data 
allocate(PressPycP(Num,Npyc,il))

allocate(UPycP(Num,Npyc,il))

allocate(VPycP(Num,Npyc,il))

! first harmonic data
allocate(PressPycH1(Num,Npyc,il))

allocate(UPycH1(Num,Npyc,il))

allocate(VPycH1(Num,Npyc,il))

! second harmonic data
allocate(PressPycH2(Num,Npyc,il))

allocate(UPycH2(Num,Npyc,il))

allocate(VPycH2(Num,Npyc,il))

! zeroth mode
allocate(PressPycZ(Num,Npyc,il))

allocate(UPycZ(Num,Npyc,il))

allocate(VPycZ(Num,Npyc,il))


! group velocities to denote the unit normal in integration
allocate(CgxPyc(Npyc,il))

allocate(CgzPyc(Npyc,il))

! energy in pycnocline
allocate(EnergyPyc(Npyc))

allocate(EnergyPycP(Npyc))

allocate(EnergyPycH1(Npyc))

allocate(EnergyPycH2(Npyc))

allocate(EnergyPycZ(Npyc))


 CgxPyc=1.
 CgzPyc=0.

!ds=2.*lambdax
do i=1,Npyc 

do j=1,il
! ds increment in x 
xpyc(i,j)=(i)*ds+Qx(j)

!z is same along the path 

zpyc(i,j)=Qz(j)

end do


end do




allocate(Np1(Npyc))

allocate(IRPyc(Npyc,il,nxr*nzr))

allocate(PointPyc(Npyc,nxr*nzr,2))


print*,"okey 5"
! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once


do i=1,Npyc

dz=5.*abs(Qz(2)-Qz(1))

dx=dz

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointPyc(i,:,:),xpyc(i,:),zpyc(i,:),il,nzr,nxr,Np1(i),dx,dz)
!print*,"Okey",Np1(i),dx,dz,h



!allocate them but it will change in every iteration 
allocate(xi(Np1(i)))

allocate(zi(Np1(i)))


do j=1,Np1(i)
 !if (i==1) then
 !print*,PointPyc(i,j,1),PointPyc(i,j,2)
 !print*,x1(PointPyc(i,j,1),PointPyc(i,j,2)),z1(PointPyc(i,j,1),PointPyc(i,j,2))
 !end if 
 xi(j)=x1(PointPyc(i,j,1),PointPyc(i,j,2))

 zi(j)=z1(PointPyc(i,j,1),PointPyc(i,j,2))

 
end do 

 call interpolation2D(IRPyc(i,:,:),xpyc(i,:),xi,zpyc(i,:),zi,Np1(i),il)

 
deallocate(xi)

deallocate(zi)



end do

print*,"Int Matrix"
print*,IRPyc(1,:,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!   Pycnocline Pressure ,Velocity,Energy         !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
do i=1,Npyc
 
  do j=1,Num


   call interpVal(IRPyc(i,:,:),PressPyc(j,i,:),Press(j,:),PointPyc(i,:,1),PointPyc(i,:,2),nxr,nzr,Np1(i),il)
 
   call interpVal(IRPyc(i,:,:),UPyc(j,i,:),u(j,:),PointPyc(i,:,1),PointPyc(i,:,2),nxr,nzr,Np1(i),il)

   call interpVal(IRPyc(i,:,:),VPyc(j,i,:),v(j,:),PointPyc(i,:,1),PointPyc(i,:,2),nxr,nzr,Np1(i),il)

  end do

 end do


print*,"okey pyc"

! here we calculate the arclength along the energy flux line 

parcpyc=0.

do i=2,il
 
 parcpyc(i)=parcpyc(i-1)+dsqrt((xpyc(1,i)-xpyc(1,i-1))**2.+(zpyc(1,i)-zpyc(1,i-1))**2.)

end do 

arcpyc=0.

do i=2,Npyc
 
 arcpyc(i)=arcpyc(i-1)+dsqrt((xpyc(i,1)-xpyc(i-1,1))**2.+(zpyc(i,1)-zpyc(i-1,1))**2.)

end do 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                             !
!      the band-pass filtering of the interpolated data       !
!                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,Npyc
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressPycP(:,i,j),PressPyc(:,i,j),PressPycFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(UPycP(:,i,j),UPyc(:,i,j),UPycFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VPycP(:,i,j),VPyc(:,i,j),VPycFour(:,i,j),w0*dt1/nor,bw/nor,Num)
 
   !  First Harmonic
  call BandPassFilter(PressPycH1(:,i,j),PressPyc(:,i,j),PressPycFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UPycH1(:,i,j),UPyc(:,i,j),UPycFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(VPycH1(:,i,j),VPyc(:,i,j),VPycFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressPycH2(:,i,j),PressPyc(:,i,j),PressPycFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UPycH2(:,i,j),UPyc(:,i,j),UPycFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VPycH2(:,i,j),VPyc(:,i,j),VPycFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

 !  Zeroth Mode
  call BandPassFilter(PressPycZ(:,i,j),PressPyc(:,i,j),PressPycFour(:,i,j),0.,bw/nor,Num)

  call BandPassFilter(UPycZ(:,i,j),UPyc(:,i,j),UPycFour(:,i,j),0.,bw/nor,Num)

   call BandPassFilter(VPycZ(:,i,j),VPyc(:,i,j),VPycFour(:,i,j),0.,bw/nor,Num)
 
 
  end do
 end do
 

print*,"okey pyc 2"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!             Total Energy flux a entry region                  !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(412,file='Tot1PUPyc200.dat',status='unknown')

do i=1,il

  write(412,*) parcpyc(i)/lambdax,PressPyc(200,1,i)*UPyc(200,1,i)/E1non
 
end do

open(413,file='Tot1PUPyc300.dat',status='unknown')

do i=1,il

  write(413,*) parcpyc(i)/lambdax,PressPyc(300,1,i)*UPyc(300,1,i)/E1non
  
end do

open(414,file='Tot1PUPyc400.dat',status='unknown')

do i=1,il

  write(414,*) parcpyc(i)/lambdax,PressPyc(400,1,i)*UPyc(400,1,i)/E1non
  
end do


! Total Energy flux after one lambdax 

open(415,file='Tot2PUPyc200.dat',status='unknown')

do i=1,il

  write(415,*) parcpyc(i)/lambdax,PressPyc(200,19,i)*UPyc(200,19,i)/E1non
 
end do

open(416,file='Tot2PUPyc300.dat',status='unknown')

do i=1,il

  write(416,*) parcpyc(i)/lambdax,PressPyc(300,19,i)*UPyc(300,19,i)/E1non
  
end do

open(417,file='Tot2PUPyc400.dat',status='unknown')

do i=1,il

  write(417,*) parcpyc(i)/lambdax,PressPyc(400,19,i)*UPyc(400,19,i)/E1non
  
end do


! Total Energy flux after two lambdax 

open(418,file='Tot3PUPyc200.dat',status='unknown')

do i=1,il

  write(418,*) parcpyc(i)/lambdax,PressPyc(200,37,i)*UPyc(200,37,i)/E1non
 
end do

open(419,file='Tot3PUPyc300.dat',status='unknown')

do i=1,il

  write(419,*) parcpyc(i)/lambdax,PressPyc(300,37,i)*UPyc(300,37,i)/E1non
  
end do

open(420,file='Tot3PUPyc400.dat',status='unknown')

do i=1,il

  write(420,*) parcpyc(i)/lambdax,PressPyc(400,37,i)*UPyc(400,37,i)/E1non
  
end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for primary frequency           !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(402,file='Pri1PUPyc200.dat',status='unknown')

do i=1,il

  write(402,*) parcpyc(i)/lambdax,PressPycP(200,1,i)*UPycP(200,1,i)/E1non
 
end do

open(403,file='Pri1PUPyc300.dat',status='unknown')

do i=1,il

  write(403,*) parcpyc(i)/lambdax,PressPycP(300,1,i)*UPycP(300,1,i)/E1non
  
end do

open(404,file='Pri1PUPyc400.dat',status='unknown')

do i=1,il

  write(404,*) parcpyc(i)/lambdax,PressPycP(400,1,i)*UPycP(400,1,i)/E1non
  
end do


! Energy flux after one lambdax 

open(405,file='Pri2PUPyc200.dat',status='unknown')

do i=1,il

  write(405,*) parcpyc(i)/lambdax,PressPycP(200,19,i)*UPycP(200,19,i)/E1non
 
end do

open(406,file='Pri2PUPyc300.dat',status='unknown')

do i=1,il

  write(406,*) parcpyc(i)/lambdax,PressPycP(300,19,i)*UPycP(300,19,i)/E1non
  
end do

open(407,file='Pri2PUPyc400.dat',status='unknown')

do i=1,il

  write(407,*) parcpyc(i)/lambdax,PressPycP(400,19,i)*UPycP(400,19,i)/E1non
  
end do


! Energy flux after two lambdax 

open(408,file='Pri3PUPyc200.dat',status='unknown')

do i=1,il

  write(408,*) parcpyc(i)/lambdax,PressPycP(200,37,i)*UPycP(200,37,i)/E1non
 
end do

open(409,file='Pri3PUPyc300.dat',status='unknown')

do i=1,il

  write(409,*) parcpyc(i)/lambdax,PressPycP(300,37,i)*UPycP(300,37,i)/E1non
  
end do

open(410,file='Pri3PUPyc400.dat',status='unknown')

do i=1,il

  write(410,*) parcpyc(i)/lambdax,PressPycP(400,37,i)*UPycP(400,37,i)/E1non
  
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for first Harmonic              !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(422,file='1Har1PUPyc200.dat',status='unknown')

do i=1,il

  write(422,*) parcpyc(i)/lambdax,PressPycH1(200,1,i)*UPycH1(200,1,i)/E1non
 
end do

open(423,file='1Har1PUPyc300.dat',status='unknown')

do i=1,il

  write(423,*) parcpyc(i)/lambdax,PressPycH1(300,1,i)*UPycH1(300,1,i)/E1non
  
end do

open(424,file='1Har1PUPyc400.dat',status='unknown')

do i=1,il

  write(424,*) parcpyc(i)/lambdax,PressPycH1(400,1,i)*UPycH1(400,1,i)/E1non
  
end do


! Energy flux after one lambdax 

open(425,file='1Har2PUPyc200.dat',status='unknown')

do i=1,il

  write(425,*) parcpyc(i)/lambdax,PressPycH1(200,19,i)*UPycH1(200,19,i)/E1non
 
end do

open(426,file='1Har2PUPyc300.dat',status='unknown')

do i=1,il

  write(426,*) parcpyc(i)/lambdax,PressPycH1(300,19,i)*UPycH1(300,19,i)/E1non
  
end do

open(427,file='1Har2PUPyc400.dat',status='unknown')

do i=1,il

  write(427,*) parcpyc(i)/lambdax,PressPycH1(400,19,i)*UPycH1(400,19,i)/E1non
  
end do


! Energy flux after two lambdax 

open(428,file='1Har3PUPyc200.dat',status='unknown')

do i=1,il

  write(428,*) parcpyc(i)/lambdax,PressPycH1(200,37,i)*UPycH1(200,37,i)/E1non
 
end do

open(429,file='1Har3PUPyc300.dat',status='unknown')

do i=1,il

  write(429,*) parcpyc(i)/lambdax,PressPycH1(300,37,i)*UPycH1(300,37,i)/E1non
  
end do

open(430,file='1Har3PUPyc400.dat',status='unknown')

do i=1,il

  write(430,*) parcpyc(i)/lambdax,PressPycH1(400,37,i)*UPycH1(400,37,i)/E1non
  
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for second Harmonic             !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(432,file='2Har1PUPyc200.dat',status='unknown')

do i=1,il

  write(432,*) parcpyc(i)/lambdax,PressPycH2(200,1,i)*UPycH2(200,1,i)/E1non
 
end do

open(433,file='2Har1PUPyc300.dat',status='unknown')

do i=1,il

  write(433,*) parcpyc(i)/lambdax,PressPycH2(300,1,i)*UPycH2(300,1,i)/E1non
  
end do

open(434,file='2Har1PUPyc400.dat',status='unknown')

do i=1,il

  write(434,*) parcpyc(i)/lambdax,PressPycH2(400,1,i)*UPycH2(400,1,i)/E1non
  
end do


! Energy flux after one lambdax 

open(435,file='2Har2PUPyc200.dat',status='unknown')

do i=1,il

  write(435,*) parcpyc(i)/lambdax,PressPycH2(200,19,i)*UPycH2(200,19,i)/E1non
 
end do

open(436,file='2Har2PUPyc300.dat',status='unknown')

do i=1,il

  write(436,*) parcpyc(i)/lambdax,PressPycH2(300,19,i)*UPycH2(300,19,i)/E1non
  
end do

open(437,file='2Har2PUPyc400.dat',status='unknown')

do i=1,il

  write(437,*) parcpyc(i)/lambdax,PressPycH2(400,19,i)*UPycH2(400,19,i)/E1non
  
end do


! Energy flux after two lambdax 

open(438,file='2Har3PUPyc200.dat',status='unknown')

do i=1,il

  write(438,*) parcpyc(i)/lambdax,PressPycH2(200,37,i)*UPycH2(200,37,i)/E1non
 
end do

open(439,file='2Har3PUPyc300.dat',status='unknown')

do i=1,il

  write(439,*) parcpyc(i)/lambdax,PressPycH2(300,37,i)*UPycH2(300,37,i)/E1non
  
end do

open(440,file='2Har3PUPyc400.dat',status='unknown')

do i=1,il

  write(440,*) parcpyc(i)/lambdax,PressPycH2(400,37,i)*UPycH2(400,37,i)/E1non
  
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for zeroth mode                 !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(452,file='Z1PUPyc200.dat',status='unknown')

do i=1,il

  write(452,*) parcpyc(i)/lambdax,PressPycZ(200,1,i)*UPycZ(200,1,i)/E1non
 
end do

open(453,file='Z1PUPyc300.dat',status='unknown')

do i=1,il

  write(453,*) parcpyc(i)/lambdax,PressPycZ(300,1,i)*UPycZ(300,1,i)/E1non
  
end do

open(454,file='Z1PUPyc400.dat',status='unknown')

do i=1,il

  write(454,*) parcpyc(i)/lambdax,PressPycZ(400,1,i)*UPycZ(400,1,i)/E1non
  
end do


! Energy flux after one lambdax 

open(455,file='Z2PUPyc200.dat',status='unknown')

do i=1,il

  write(455,*) parcpyc(i)/lambdax,PressPycZ(200,19,i)*UPycZ(200,19,i)/E1non
 
end do

open(456,file='Z2PUPyc300.dat',status='unknown')

do i=1,il

  write(456,*) parcpyc(i)/lambdax,PressPycZ(300,19,i)*UPycZ(300,19,i)/E1non
  
end do

open(457,file='Z2PUPyc400.dat',status='unknown')

do i=1,il

  write(457,*) parcpyc(i)/lambdax,PressPycZ(400,19,i)*UPycZ(400,19,i)/E1non
  
end do


! Energy flux after two lambdax 

open(458,file='Z3PUPyc200.dat',status='unknown')

do i=1,il

  write(458,*) parcpyc(i)/lambdax,PressPycZ(200,37,i)*UPycZ(200,37,i)/E1non
 
end do

open(459,file='Z3PUPyc300.dat',status='unknown')

do i=1,il

  write(459,*) parcpyc(i)/lambdax,PressPycZ(300,37,i)*UPycZ(300,37,i)/E1non
  
end do

open(460,file='Z3PUPyc400.dat',status='unknown')

do i=1,il

  write(460,*) parcpyc(i)/lambdax,PressPycZ(400,37,i)*UPycZ(400,37,i)/E1non
  
end do


print*,"okey pyc3"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
!            Energy Integral along Pycnocline                    !
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Integral for t=200

do i=1,Npyc

 !Total Energy 
 call PressVelocityInteg(EnergyPyc(i),PressPyc(200,i,:),UPyc(200,i,:),VPyc(200,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

 ! Primary Frequency Energy 
 call PressVelocityInteg(EnergyPycP(i),PressPycP(200,i,:),UPycP(200,i,:),VPycP(200,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

 ! First Harmonic Energy
 call PressVelocityInteg(EnergyPycH1(i),PressPycH1(200,i,:),UPycH1(200,i,:),VPycH1(200,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
 ! Second Harmonic Energy
 call PressVelocityInteg(EnergyPycH2(i),PressPycH2(200,i,:),UPycH2(200,i,:),VPycH2(200,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
 ! Zeroth Mode Energy
 call PressVelocityInteg(EnergyPycZ(i),PressPycZ(200,i,:),UPycZ(200,i,:),VPycZ(200,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
end do 



open(500,file='EnergyPyc200.dat',status='unknown')

do i=1,Npyc

  write(500,*) arcpyc(i)/lambdax,EnergyPyc(i)/Enon
 
end do


open(501,file='PEnergyPyc200.dat',status='unknown')

do i=1,Npyc

  write(501,*) arcpyc(i)/lambdax,EnergyPycP(i)/Enon
 
end do


open(502,file='1HEnergyPyc200.dat',status='unknown')

do i=1,Npyc

  write(502,*) arcpyc(i)/lambdax,EnergyPycH1(i)/Enon
 
end do


open(503,file='2HEnergyPyc200.dat',status='unknown')

do i=1,Npyc

  write(503,*) arcpyc(i)/lambdax,EnergyPycH2(i)/Enon
 
end do


open(520,file='ZEnergyPyc200.dat',status='unknown')

do i=1,Npyc

  write(520,*) arcpyc(i)/lambdax,EnergyPycZ(i)/Enon
 
end do

! Integral for t=300

do i=1,Npyc

 !Total Energy 
 call PressVelocityInteg(EnergyPyc(i),PressPyc(300,i,:),UPyc(300,i,:),VPyc(300,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

 ! Primary Frequency Energy 
 call PressVelocityInteg(EnergyPycP(i),PressPycP(300,i,:),UPycP(300,i,:),VPycP(300,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

 ! First Harmonic Energy
 call PressVelocityInteg(EnergyPycH1(i),PressPycH1(300,i,:),UPycH1(300,i,:),VPycH1(300,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
 ! Second Harmonic Energy
 call PressVelocityInteg(EnergyPycH2(i),PressPycH2(300,i,:),UPycH2(300,i,:),VPycH2(300,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
 ! Zeroth Mode Energy
 call PressVelocityInteg(EnergyPycZ(i),PressPycZ(300,i,:),UPycZ(300,i,:),VPycZ(300,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

end do 



open(504,file='EnergyPyc300.dat',status='unknown')

do i=1,Npyc

  write(504,*) arcpyc(i)/lambdax,EnergyPyc(i)/Enon
 
end do


open(505,file='PEnergyPyc300.dat',status='unknown')

do i=1,Npyc

  write(505,*) arcpyc(i)/lambdax,EnergyPycP(i)/Enon
 
end do


open(507,file='1HEnergyPyc300.dat',status='unknown')

do i=1,Npyc

  write(507,*) arcpyc(i)/lambdax,EnergyPycH1(i)/Enon
 
end do


open(508,file='2HEnergyPyc300.dat',status='unknown')

do i=1,Npyc

  write(508,*) arcpyc(i)/lambdax,EnergyPycH2(i)/Enon
 
end do

open(521,file='ZEnergyPyc300.dat',status='unknown')

do i=1,Npyc

  write(521,*) arcpyc(i)/lambdax,EnergyPycZ(i)/Enon
 
end do


! Integral for t=400

do i=1,Npyc

 !Total Energy 
 call PressVelocityInteg(EnergyPyc(i),PressPyc(400,i,:),UPyc(400,i,:),VPyc(400,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

 ! Primary Frequency Energy 
 call PressVelocityInteg(EnergyPycP(i),PressPycP(400,i,:),UPycP(400,i,:),VPycP(400,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)

 ! First Harmonic Energy
 call PressVelocityInteg(EnergyPycH1(i),PressPycH1(400,i,:),UPycH1(400,i,:),VPycH1(400,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
 ! Second Harmonic Energy
 call PressVelocityInteg(EnergyPycH2(i),PressPycH2(400,i,:),UPycH2(400,i,:),VPycH2(400,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 
 ! Zeroth Mode Energy
 call PressVelocityInteg(EnergyPycZ(i),PressPycZ(400,i,:),UPycZ(400,i,:),VPycZ(400,i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)
 

end do 



open(509,file='EnergyPyc400.dat',status='unknown')

do i=1,Npyc

  write(509,*) arcpyc(i)/lambdax,EnergyPyc(i)/Enon
 
end do


open(510,file='PEnergyPyc400.dat',status='unknown')

do i=1,Npyc

  write(510,*) arcpyc(i)/lambdax,EnergyPycP(i)/Enon
 
end do


open(511,file='1HEnergyPyc400.dat',status='unknown')

do i=1,Npyc

  write(511,*) arcpyc(i)/lambdax,EnergyPycH1(i)/Enon
 
end do


open(512,file='2HEnergyPyc400.dat',status='unknown')

do i=1,Npyc

  write(512,*) arcpyc(i)/lambdax,EnergyPycH2(i)/Enon
 
end do

open(522,file='ZEnergyPyc400.dat',status='unknown')

do i=1,Npyc

  write(522,*) arcpyc(i)/lambdax,EnergyPycZ(i)/Enon
 
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
!            The Calculation of Reradiating Energy               !  
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! In order to calculate the reradiatin energy which is basically 
! pressure and z-velocity at the lower bound of the pycnocline 

allocate(EnergyRad(Npyc))
allocate(EnergyRadP(Npyc))
allocate(EnergyRadH1(Npyc))
allocate(EnergyRadH2(Npyc))

allocate(DEnergyRad(Npyc))
allocate(DEnergyRadP(Npyc))
allocate(DEnergyRadH1(Npyc))
allocate(DEnergyRadH2(Npyc))



! Right now smple trapezoid rule 

! t=200

EnergyRad=0.

DEnergyRad=0.

EnergyRadP=0.

DEnergyRadP=0.

EnergyRadH1=0.

DEnergyRadH1=0.

EnergyRadH2=0.

DEnergyRadH2=0.

! total energy 

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRad(i)= EnergyRad(i-1)-0.5*ds*(PressPyc(200,i,1)*VPyc(200,i,1)+PressPyc(200,i-1,1)*VPyc(200,i-1,1))
 

end do 

! energy of primary frequency

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadP(i)= EnergyRadP(i-1)-0.5*ds*(PressPycP(200,i,1)*VPycP(200,i,1)+PressPycP(200,i-1,1)*VPycP(200,i-1,1))
 

end do 

! energy of first harmonic

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH1(i)= EnergyRadH1(i-1)-0.5*ds*(PressPycH1(200,i,1)*VPycH1(200,i,1)+PressPycH1(200,i-1,1)*VPycH1(200,i-1,1))
 

end do 

! energy of second harmonic

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH2(i)= EnergyRadH2(i-1)-0.5*ds*(PressPycH2(200,i,1)*VPycH2(200,i,1)+PressPycH2(200,i-1,1)*VPycH2(200,i-1,1))
 

end do 

! the derivative of the reradiating energy 

! total energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRad(i)= -1.*PressPyc(200,i,1)*VPyc(200,i,1)
 

end do 

! primary energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadP(i)= -1.*PressPycP(200,i,1)*VPycP(200,i,1)
 

end do 

! first harmonics energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadH1(i)= -1.*PressPycH1(200,i,1)*VPycH1(200,i,1)
 

end do 

! second harmonics energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadH2(i)= -1.*PressPycH2(200,i,1)*VPycH2(200,i,1)
 

end do 


open(600,file='EnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(600,*) arcpyc(i)/lambdax,EnergyRad(i)/Enon
  
end do

open(601,file='DEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(601,*) arcpyc(i)/lambdax,DEnergyRad(i)/E1non
 
end do


open(602,file='PEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(602,*) arcpyc(i)/lambdax,EnergyRadP(i)/Enon
  
end do

open(603,file='PDEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(603,*) arcpyc(i)/lambdax,DEnergyRadP(i)/E1non
 
end do


open(604,file='1HEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(604,*) arcpyc(i)/lambdax,EnergyRadH1(i)/Enon
  
end do

open(605,file='1HDEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(605,*) arcpyc(i)/lambdax,DEnergyRadH1(i)/E1non
 
end do


open(606,file='2HEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(606,*) arcpyc(i)/lambdax,EnergyRadH2(i)/Enon
  
end do

open(607,file='2HDEnergyRad200.dat',status='unknown')

do i=1,Npyc

  write(607,*) arcpyc(i)/lambdax,DEnergyRadH2(i)/E1non
 
end do


! t=300

EnergyRad=0.

DEnergyRad=0.

EnergyRadP=0.

DEnergyRadP=0.

EnergyRadH1=0.

DEnergyRadH1=0.

EnergyRadH2=0.

DEnergyRadH2=0.

! total energy 

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRad(i)= EnergyRad(i-1)-0.5*ds*(PressPyc(300,i,1)*VPyc(300,i,1)+PressPyc(300,i-1,1)*VPyc(300,i-1,1))
 

end do 

! energy of primary frequency

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadP(i)= EnergyRadP(i-1)-0.5*ds*(PressPycP(300,i,1)*VPycP(300,i,1)+PressPycP(300,i-1,1)*VPycP(300,i-1,1))
 

end do 

! energy of first harmonic

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH1(i)= EnergyRadH1(i-1)-0.5*ds*(PressPycH1(300,i,1)*VPycH1(300,i,1)+PressPycH1(300,i-1,1)*VPycH1(300,i-1,1))
 

end do 

! energy of second harmonic

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH2(i)= EnergyRadH2(i-1)-0.5*ds*(PressPycH2(300,i,1)*VPycH2(300,i,1)+PressPycH2(300,i-1,1)*VPycH2(300,i-1,1))
 

end do 

! the derivative of the reradiating energy 

! total energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRad(i)= -1.*PressPyc(300,i,1)*VPyc(300,i,1)
 

end do 

! primary energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadP(i)= -1.*PressPycP(300,i,1)*VPycP(300,i,1)
 

end do 

! first harmonics energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadH1(i)= -1.*PressPycH1(300,i,1)*VPycH1(300,i,1)
 

end do 

! second harmonics energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadH2(i)= -1.*PressPycH2(300,i,1)*VPycH2(300,i,1)
 

end do 


open(610,file='EnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(610,*) arcpyc(i)/lambdax,EnergyRad(i)/Enon
  
end do

open(611,file='DEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(611,*) arcpyc(i)/lambdax,DEnergyRad(i)/E1non
 
end do


open(612,file='PEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(612,*) arcpyc(i)/lambdax,EnergyRadP(i)/Enon
  
end do

open(613,file='PDEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(613,*) arcpyc(i)/lambdax,DEnergyRadP(i)/E1non
 
end do


open(614,file='1HEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(614,*) arcpyc(i)/lambdax,EnergyRadH1(i)/Enon
  
end do

open(615,file='1HDEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(615,*) arcpyc(i)/lambdax,DEnergyRadH1(i)/E1non
 
end do


open(616,file='2HEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(616,*) arcpyc(i)/lambdax,EnergyRadH2(i)/Enon
  
end do

open(617,file='2HDEnergyRad300.dat',status='unknown')

do i=1,Npyc

  write(617,*) arcpyc(i)/lambdax,DEnergyRadH2(i)/E1non
 
end do





! t=400

EnergyRad=0.

DEnergyRad=0.

EnergyRadP=0.

DEnergyRadP=0.

EnergyRadH1=0.

DEnergyRadH1=0.

EnergyRadH2=0.

DEnergyRadH2=0.

! total energy 

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRad(i)= EnergyRad(i-1)-0.5*ds*(PressPyc(400,i,1)*VPyc(400,i,1)+PressPyc(400,i-1,1)*VPyc(400,i-1,1))
 

end do 

! energy of primary frequency

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadP(i)= EnergyRadP(i-1)-0.5*ds*(PressPycP(400,i,1)*VPycP(400,i,1)+PressPycP(400,i-1,1)*VPycP(400,i-1,1))
 

end do 

! energy of first harmonic

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH1(i)= EnergyRadH1(i-1)-0.5*ds*(PressPycH1(400,i,1)*VPycH1(400,i,1)+PressPycH1(400,i-1,1)*VPycH1(400,i-1,1))
 

end do 

! energy of second harmonic

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH2(i)= EnergyRadH2(i-1)-0.5*ds*(PressPycH2(400,i,1)*VPycH2(400,i,1)+PressPycH2(400,i-1,1)*VPycH2(400,i-1,1))
 

end do 

! the derivative of the reradiating energy 

! total energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRad(i)= -1.*PressPyc(400,i,1)*VPyc(400,i,1)
 

end do 

! primary energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadP(i)= -1.*PressPycP(400,i,1)*VPycP(400,i,1)
 

end do 

! first harmonics energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadH1(i)= -1.*PressPycH1(400,i,1)*VPycH1(400,i,1)
 

end do 

! second harmonics energy 

do i=1,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  DEnergyRadH2(i)= -1.*PressPycH2(400,i,1)*VPycH2(400,i,1)
 

end do 


open(620,file='EnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(620,*) arcpyc(i)/lambdax,EnergyRad(i)/Enon
  
end do

open(621,file='DEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(621,*) arcpyc(i)/lambdax,DEnergyRad(i)/E1non
 
end do


open(622,file='PEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(622,*) arcpyc(i)/lambdax,EnergyRadP(i)/Enon
  
end do

open(623,file='PDEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(623,*) arcpyc(i)/lambdax,DEnergyRadP(i)/E1non
 
end do


open(624,file='1HEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(624,*) arcpyc(i)/lambdax,EnergyRadH1(i)/Enon
  
end do

open(625,file='1HDEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(625,*) arcpyc(i)/lambdax,DEnergyRadH1(i)/E1non
 
end do


open(626,file='2HEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(626,*) arcpyc(i)/lambdax,EnergyRadH2(i)/Enon
  
end do

open(627,file='2HDEnergyRad400.dat',status='unknown')

do i=1,Npyc

  write(627,*) arcpyc(i)/lambdax,DEnergyRadH2(i)/E1non
 
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           !
!      The Calculations of the Mean Properties              !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! the interpolated data in pycnocline
allocate(PressPycM(Npyc,il))

allocate(UPycM(Npyc,il))

allocate(VPycM(Npyc,il))

! the mean energy flux

allocate(EnergyPycM(Npyc))

do i=1,il

  do j=1,Npyc
   
   dummy=0.

   dummy1=0.

   dummy2=0.

     do k=1,Num-1

        dummy=dummy+0.5*(PressPyc(k,j,i)+PressPyc(k,j,i+1))*(timec(i+1)-timec(i))

        dummy1=dummy1+0.5*(UPyc(k,j,i)+UPyc(k,j,i+1))*(timec(i+1)-timec(i))
       
        dummy2=dummy2+0.5*(VPyc(k,j,i)+VPyc(k,j,i+1))*(timec(i+1)-timec(i))

     end do

  PressPycM(j,i)=dummy/timec(Num) 

  UPycM(j,i)=dummy1/timec(Num) 

  VPycM(j,i)=dummy2/timec(Num)   

  end do 

end do

! the integration of the mean energy

do i=1,Npyc

 ! Mean Energy 
 call PressVelocityInteg(EnergyPycM(i),PressPycM(i,:),UPycM(i,:),VPycM(i,:),il,CgxPyc(i,:),CgzPyc(i,:),h)


end do 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!             Energy flux of The mean flow                      !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(902,file='1PUPycM.dat',status='unknown')

do i=1,il

  write(902,*) parcpyc(i)/lambdax,PressPycM(1,i)*UPycM(1,i)/E1non
 
end do

open(903,file='2PUPycM.dat',status='unknown')

do i=1,il

  write(903,*) parcpyc(i)/lambdax,PressPycM(19,i)*UPycM(19,i)/E1non
  
end do

open(904,file='3PUPycM.dat',status='unknown')

do i=1,il

  write(904,*) parcpyc(i)/lambdax,PressPycM(37,i)*UPycM(37,i)/E1non
  
end do

! the energy integral of the mean flow

open(905,file='EnergyMean.dat',status='unknown')

do i=1,Npyc

  write(905,*) arcpyc(i)/lambdax,EnergyPycM(i)/Enon
  
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           ! 
!           Averaged Energy Flux over One Wave Period       ! 
!                                                           !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(PUT(Npyc,il))

allocate(PUZ(Npyc,il))

allocate(PUP(Npyc,il))

allocate(PUH1(Npyc,il))

allocate(PUH2(Npyc,il))



! let's generate the averaged pressure and the velocity product

do i=1,Npyc

 do j=1,il
 
  ! we are sending 107 time steps since it is equla to one wave period
  call  MeanPU(PUT(i,j),PressPyc(ti:tf,i,j),UPyc(ti:tf,i,j),VPyc(ti:tf,i,j),CgxPyc(i,j),CgzPyc(i,j),ttot)
  
  call  MeanPU(PUZ(i,j),PressPycZ(ti:tf,i,j),UPycZ(ti:tf,i,j),VPycZ(ti:tf,i,j),CgxPyc(i,j),CgzPyc(i,j),ttot)

  call  MeanPU(PUP(i,j),PressPycP(ti:tf,i,j),UPycP(ti:tf,i,j),VPycP(ti:tf,i,j),CgxPyc(i,j),CgzPyc(i,j),ttot)

  call  MeanPU(PUH1(i,j),PressPycH1(ti:tf,i,j),UPycH1(ti:tf,i,j),VPycH1(ti:tf,i,j),CgxPyc(i,j),CgzPyc(i,j),ttot)
  
  call  MeanPU(PUH2(i,j),PressPycH2(ti:tf,i,j),UPycH2(ti:tf,i,j),VPycH2(ti:tf,i,j),CgxPyc(i,j),CgzPyc(i,j),ttot)

 end do 

end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!             Total Energy flux a entry region                  !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(912,file='Tot1PUPycAve.dat',status='unknown')

do i=1,il

  write(912,*) parcpyc(i)/lambdax,PUT(1,i)/E1non
 
end do



! Total Energy flux after one lambdax 

open(913,file='Tot2PUPycAve.dat',status='unknown')

do i=1,il

  write(913,*) parcpyc(i)/lambdax,PUT(19,i)/E1non
 
end do



! Total Energy flux after two lambdax 

open(914,file='Tot3PUPycAve.dat',status='unknown')

do i=1,il

  write(914,*) parcpyc(i)/lambdax,PUT(37,i)/E1non
 
end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for primary frequency           !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(915,file='P1PUPycAve.dat',status='unknown')

do i=1,il

  write(915,*) parcpyc(i)/lambdax,PUP(1,i)/E1non
 
end do



! Primary Energy flux after one lambdax 

open(916,file='P2PUPycAve.dat',status='unknown')

do i=1,il

  write(916,*) parcpyc(i)/lambdax,PUP(19,i)/E1non
 
end do



! Primary Energy flux after two lambdax 

open(917,file='P3PUPycAve.dat',status='unknown')

do i=1,il

  write(917,*) parcpyc(i)/lambdax,PUP(37,i)/E1non
 
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for first Harmonic              !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(918,file='1Har1PUPycAve.dat',status='unknown')

do i=1,il

  write(918,*) parcpyc(i)/lambdax,PUH1(1,i)/E1non
 
end do



! First harmonic Energy flux after one lambdax 

open(919,file='1Har2PUPycAve.dat',status='unknown')

do i=1,il

  write(919,*) parcpyc(i)/lambdax,PUH1(19,i)/E1non
 
end do



! First harmonic Energy flux after two lambdax 

open(920,file='1Har3PUPycAve.dat',status='unknown')

do i=1,il

  write(920,*) parcpyc(i)/lambdax,PUH1(37,i)/E1non
 
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for second Harmonic             !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(921,file='2Har1PUPycAve.dat',status='unknown')

do i=1,il

  write(921,*) parcpyc(i)/lambdax,PUH2(1,i)/E1non
 
end do



!Second harmonic Energy flux after one lambdax 

open(922,file='2Har2PUPycAve.dat',status='unknown')

do i=1,il

  write(922,*) parcpyc(i)/lambdax,PUH2(19,i)/E1non
 
end do



! Second harmonic Energy flux after two lambdax 

open(922,file='2Har3PUPycAve.dat',status='unknown')

do i=1,il

  write(922,*) parcpyc(i)/lambdax,PUH2(37,i)/E1non
 
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for zeroth mode                 !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(923,file='Z1PUPycAve.dat',status='unknown')

do i=1,il

  write(923,*) parcpyc(i)/lambdax,PUZ(1,i)/E1non
 
end do



! First harmonic Energy flux after one lambdax 

open(924,file='Z2PUPycAve.dat',status='unknown')

do i=1,il

  write(924,*) parcpyc(i)/lambdax,PUZ(19,i)/E1non
 
end do



! First harmonic Energy flux after two lambdax 

open(925,file='Z3PUPycAve.dat',status='unknown')

do i=1,il

  write(925,*) parcpyc(i)/lambdax,PUZ(37,i)/E1non
 
end do


! Integral of the averaged energy 

do i=1,Npyc

 !Total Energy 
  call RomBerg(PUT(i,:),EnergyPyc(i),0.5*h,il)

 ! Primary Frequency Energy 
 call RomBerg(PUP(i,:),EnergyPycP(i),0.5*h,il)

 ! First Harmonic Energy
 call RomBerg(PUH1(i,:),EnergyPycH1(i),0.5*h,il)
 
 ! Second Harmonic Energy
 call RomBerg(PUH2(i,:),EnergyPycH2(i),0.5*h,il)
 
 ! Zeroth Mode Energy
 call RomBerg(PUZ(i,:),EnergyPycZ(i),0.5*h,il)
 

end do 



open(809,file='EnergyPycAve.dat',status='unknown')

do i=1,Npyc

  write(809,*) arcpyc(i)/lambdax,EnergyPyc(i)/Enon
 
end do


open(810,file='PEnergyPycAve.dat',status='unknown')

do i=1,Npyc

  write(810,*) arcpyc(i)/lambdax,EnergyPycP(i)/Enon
 
end do


open(811,file='1HEnergyPycAve.dat',status='unknown')

do i=1,Npyc

  write(811,*) arcpyc(i)/lambdax,EnergyPycH1(i)/Enon
 
end do


open(812,file='2HEnergyPycAve.dat',status='unknown')

do i=1,Npyc

  write(812,*) arcpyc(i)/lambdax,EnergyPycH2(i)/Enon
 
end do

open(822,file='ZEnergyPycAve.dat',status='unknown')

do i=1,Npyc

  write(822,*) arcpyc(i)/lambdax,EnergyPycZ(i)/Enon
 
end do

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
!    Averaged Radiated Energy from Pycnocline              !
!                                                          ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



allocate(PURadT(Npyc))

allocate(PURadZ(Npyc))

allocate(PURadP(Npyc))

allocate(PURadH1(Npyc))

allocate(PURadH2(Npyc))


! let's generate the averaged pressure and the velocity product

do i=1,Npyc

  ! we are sending 107 time steps since it is equla to one wave period
  call  MeanPU(PURadT(i),PressPyc(ti:tf,i,1),UPyc(ti:tf,i,1),VPyc(ti:tf,i,1),0.,1.,ttot)
  
  call  MeanPU(PURadZ(i),PressPycZ(ti:tf,i,1),UPycZ(ti:tf,i,1),VPycZ(ti:tf,i,1),0.,1.,ttot)

  call  MeanPU(PURadP(i),PressPycP(ti:tf,i,1),UPycP(ti:tf,i,1),VPycP(ti:tf,i,1),0.,1.,ttot)

  call  MeanPU(PURadH1(i),PressPycH1(ti:tf,i,1),UPycH1(ti:tf,i,1),VPycH1(ti:tf,i,1),0.,1.,ttot)
  
  call  MeanPU(PURadH2(i),PressPycH2(ti:tf,i,1),UPycH2(ti:tf,i,1),VPycH2(ti:tf,i,1),0.,1.,ttot)


end do



EnergyRad=0.

EnergyRadP=0.

EnergyRadH1=0.

EnergyRadH2=0.


do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRad(i)= EnergyRad(i-1)-0.5*ds*(PURadT(i-1)+PURadT(i))
 

end do 

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadP(i)= EnergyRadP(i-1)-0.5*ds*(PURadP(i-1)+PURadP(i))
 

end do 

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH1(i)= EnergyRadH1(i-1)-0.5*ds*(PURadH1(i-1)+PURadH1(i))
 

end do 

do i=2,Npyc
! number 1 is used for last indice as it is the lower bound of the Pycnocline 
  EnergyRadH2(i)= EnergyRadH2(i-1)-0.5*ds*(PURadH2(i-1)+PURadH2(i))
 

end do 


open(830,file='EnergyRadAve.dat',status='unknown')

do i=1,Npyc

  write(830,*) arcpyc(i)/lambdax,EnergyRad(i)/Enon
  
end do


open(831,file='PEnergyRadAve.dat',status='unknown')

do i=1,Npyc

  write(831,*) arcpyc(i)/lambdax,EnergyRadP(i)/Enon
  
end do


open(832,file='1HEnergyRadAve.dat',status='unknown')

do i=1,Npyc

  write(832,*) arcpyc(i)/lambdax,EnergyRadH1(i)/Enon
  
end do

open(833,file='2HEnergyRadAve.dat',status='unknown')

do i=1,Npyc

  write(833,*) arcpyc(i)/lambdax,EnergyRadH2(i)/Enon
  
end do



 


end program EnergyIntegral
