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
   ! the coordinates in triangular region reflection
  real*8,allocatable ::xtrir(:,:),ztrir(:,:)
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
    
 !the group velocities in triangular region 
  real*8 ,allocatable::Cgxtrir(:,:),Cgztrir(:,:)

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
  integer ,allocatable::Point(:,:,:),PointPyc(:,:,:),PointRef(:,:,:),PointTri(:,:,:),PointTriR(:,:,:)
   !The actual Pressure and velocities interpolated from data files in time
  real*8, allocatable::Press(:,:),u(:,:),v(:,:)
   !The actual Pressure and velocities read from data files
  real*8, allocatable::PressD(:,:),uD(:,:),vD(:,:)
  !The interpolated Pressure
  real*8, allocatable::PressInt(:,:,:),PressIntP(:,:,:),PressIntH(:,:,:)
  !The interpolated pressure in pycnocline
  real*8,allocatable:: PressIntPyc(:,:,:),PressIntPycP(:,:,:),PressIntPycH(:,:,:),PressIntPycH2(:,:,:)
  !The interpolated pressure in pycnocline
  real*8,allocatable:: PressRef(:,:,:),PressRefP(:,:,:),PressRefH1(:,:,:),PressRefH2(:,:,:)
   !The interpolated pressure in pycnocline
  real*8,allocatable:: PressIntTri(:,:,:),PressIntTriP(:,:,:),PressIntTriH(:,:,:),PressIntTriH2(:,:,:)
    !The interpolated pressure in triangular region
  real*8,allocatable:: PressTri(:,:,:),PressTriP(:,:,:),PressTriH1(:,:,:),PressTriH2(:,:,:)
     !The interpolated pressure in triangular region
  real*8,allocatable:: PressTriR(:,:,:),PressTriRP(:,:,:),PressTriRH1(:,:,:),PressTriRH2(:,:,:)  
  ! the fourier transform of pressure
  complex , allocatable ::PressRefFour(:,:,:),PressIntFour(:,:,:),PressIntPycFour(:,:,:)
   ! the fourier transform of pressure
  complex , allocatable :: PressTriFour(:,:,:),PressTriRFour(:,:,:)

   !The interpolation points 
  real*8,allocatable::xi(:),zi(:),Pin(:)
   !The interpolated velocities
  real*8,allocatable::uin(:),vin(:)
   ! The array used to keep the interpolated velocities every iteration and fourier transform of it 
   real*8, allocatable::UInt(:,:,:),VInt(:,:,:)
   ! The array used to keep primary wave and higher order harmonics data 
   real*8, allocatable::UIntP(:,:,:),VIntP(:,:,:),UIntH(:,:,:),VIntH(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntPyc(:,:,:),VIntPyc(:,:,:)
  ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntTri(:,:,:),VIntTri(:,:,:)
   ! The array used to keep the interpolated velocities in triangular region 
   real*8, allocatable::UTri(:,:,:),VTri(:,:,:),UTriP(:,:,:),VTriP(:,:,:)
    ! The array used to keep the interpolated velocities in triangular region 
   real*8, allocatable::UTriH1(:,:,:),VTriH1(:,:,:),UTriH2(:,:,:),VTriH2(:,:,:)

  ! The array used to keep the interpolated velocities in triangular reflection region 
   real*8, allocatable::UTriR(:,:,:),VTriR(:,:,:),UTriRP(:,:,:),VTriRP(:,:,:)
    ! The array used to keep the interpolated velocities in triangular reflection region 
   real*8, allocatable::UTriRH1(:,:,:),VTriRH1(:,:,:),UTriRH2(:,:,:),VTriRH2(:,:,:)

   ! The array used to keep the interpolated pressure velocity produts in triangular zone 
   real*8, allocatable::PUTriRT(:,:),PUTriRZ(:,:),PUTriRP(:,:),PUTriRH1(:,:),PUTriRH2(:,:)
! The array used to keep the interpolated velocities in triangular region
   real*8, allocatable::UIntTriP(:,:,:),VIntTriP(:,:,:),UIntTriH(:,:,:),VIntTriH(:,:,:),UIntTriH2(:,:,:),VIntTriH2(:,:,:)
!The array use to keep the fourier transform of the velociies in the triangular region
  real*8 ,allocatable::UTriFour(:,:,:),VTriFour(:,:,:),UTriRFour(:,:,:),VTriRFour(:,:,:)
  ! The array used to keep the interpolated velocities in reflection region
   real*8, allocatable::URef(:,:,:),VRef(:,:,:)
  ! The array used to keep the interpolated velocities in reflection region
   real*8, allocatable::URefP(:,:,:),VRefP(:,:,:),URefH1(:,:,:),VRefH1(:,:,:),URefH2(:,:,:),VRefH2(:,:,:)
  ! the fourier transform of the velocities
   complex , allocatable ::URefFour(:,:,:),VRefFour(:,:,:),UIntPycFour(:,:,:),VIntPycFour(:,:,:)
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
   real*8,allocatable:: EnergyPyc(:,:),EnergyPycP(:,:),EnergyPycH(:,:),EnergyPycH2(:,:)
   !the energ flux of the beam in triangular region
   real*8,allocatable:: EnergyTriP(:),EnergyTri(:),EnergyTriH1(:),EnergyTriH2(:)
   !the energ flux of the beam in triangular region
   real*8,allocatable:: EnergyTriRP(:),EnergyTriR(:),EnergyTriRH1(:),EnergyTriRH2(:)
   !the energ flux of the beam in reflection zone
   real*8,allocatable:: EnergyRef(:),EnergyRefP(:),EnergyRefH1(:),EnergyRefH2(:)
   
   !the array where all the interpolation matrices to be stored
   real*8,allocatable::IR(:,:,:),IRPyc(:,:,:),IRRef(:,:,:),IRTri(:,:,:),IRTriR(:,:,:)
   !the array where the number of points to be used in each interpolation stored
   integer,allocatable::Np(:)
   !the array where the number of points to be used in each interpolation stored in pycnocline
   integer,allocatable::Np1(:)
    !the array where the number of points to be used in each interpolation stored in reflection region
   integer,allocatable::Np2(:)
    !the array where the number of points to be used in each interpolation stored in triangular region 
   integer,allocatable::NpTri(:),Npref(:),NpTri2(:),NpTriR(:)
   ! the width of the integration path in triangular region where reflection occurs
   real*8 ,allocatable::zdimr(:),zdimrh(:),dr(:)
   ! the arclength path in triangular region where reflection occurs
   real*8 ,allocatable::arctri(:),arctrih(:),arctrir(:),arcref(:)
   ! the isoenergy line in triangular region where reflection occurs
   real*8 ,allocatable::parctri(:,:),parctrih(:),parctrir(:,:),parcref(:)
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
   real*8 ,allocatable::time(:)
   !the time steps to be used 
   integer ,allocatable::Nums(:)
   ! the non-dimensionaliztion parameters
   real*8 rho,lambdax,u0,Enon,Pnon,E1non

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rho=1E3
  
  lambdax=0.0875

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

   !r=2.
   !dt=0.0344188927665368

   !r=4.
   !dt=   0.0287606935364547

    !r=6.
    !dt=0.0270295880980098 
    
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

  ! the filter normalization parameter
     nor=10.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! the number of the pressure and velocity data files
  NumDat=200

  Num=512

  ! il has to be 2**N+1 since it is used in romberg integration 
  il=17

  ! the output time 
 ! Tout=300
   
  allocate(time(NumDat))

   allocate(Nums(NumDat))

  allocate(Pfiles(NumDat))

  allocate(Vfiles(NumDat))


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

   print*,"pycnocline thickness",h

allocate(Qx(il))
allocate(Qz(il))
allocate(N(il))



xlen=x1(nxr,1)

 ! here we generate the initial quadrature points for
 ! for the wave

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

print*," okey 3"
ds=0.0001

do i=1,Ntri

 zdimr(i)=0.

 !x(N(il)+i,1)

 zdum=z(N(il)+i,1)


 ! in this loop , we find the length of the integration path 
  do j=1,2000 
 
   zdum=zdum+ds

   call BVprofile(BV,N1,r,zdum,z0,dt)

    if(zdum>z(N(il),il)) then 

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

open(300,file='TriPath40.dat',status='unknown')

do i=1,il

  write(300,*) xtri(40,i),ztri(40,i)
  !print*,N(i)
end do

open(301,file='FPointLineP.dat',status='unknown')

do i=1,Ntri

  write(301,*) xtri(i,1),ztri(i,1)
  !print*,N(i)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       ! 
!     The reflection lines in triangular region         !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
allocate(xtrir(Ntri,il))

allocate(ztrir(Ntri,il))


allocate(Cgxtrir(Ntri,il))

allocate(Cgztrir(Ntri,il))

allocate(arctrir(Ntri))

allocate(parctrir(Ntri,il))

allocate(dr(Ntri))

!refletion path after triangular region 
allocate(xref(N(il),il))

allocate(zref(N(il),il))


allocate(CgxRef(N(il),il))

allocate(CgzRef(N(il),il))

allocate(arcref(N(il)))

allocate(parcref(il))

! the center for reflection to be used to locate reflection lines

! in that reflection region points will be symmetric wtih respec to xcenr

xcenr=0.5*(x(N(1),1)+x(N(il),il))

do i=1,Ntri 
 do j=1,il

  ! that reflection only affects the x-coordinate

  xtrir(i,j)=2.*xcenr-xtri(Ntri-i+1,j)

  ztrir(i,j)=ztri(Ntri-i+1,j)

  ! at reflection region Cgx vector points out the same direction  

  Cgxtrir(i,j)=Cgxtri(Ntri-i+1,j)

  Cgztrir(i,j)=-1.*Cgztri(Ntri-i+1,j)

 end do

 dr(i)=zdimr(Ntri-i+1)/2.

end do 


do i=1,N(il)
 do j=1,il

  ! that reflection only affects the x-coordinate

  xref(i,j)=2.*xcenr-x(N(il)-i+1,j)

  zref(i,j)=z(N(il)-i+1,j)

  ! at reflection region Cgx vector points out the same direction  

  CgxRef(i,j)=Cgx(N(il)-i+1,j)

  CgzRef(i,j)=-1.*Cgz(N(il)-i+1,j)

 end do
end do 


open(400,file='TriPathR40.dat',status='unknown')

do i=1,Ntri

  write(400,*) xtrir(i,1),ztrir(i,1)
  !print*,N(i)
end do

open(401,file='RefPathR40.dat',status='unknown')

do i=1,N(il)

  write(401,*) x(i,1),z(i,1)
  !print*,N(i)
end do



print*,"Primary and harmonic parameters"

print*,"wave numbers",kx

print*,"frequencies",w0,Ntri


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
!    Interpolation in Triangular Reflection Region         !
!                                                          !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(Nptrir(N(il)))

allocate(IRTriR(Ntri,il,nxr*nzr))

allocate(PointTriR(Ntri,nxr*nzr,2))

! total presurre and velocity
allocate(PressTriR(Num,Ntri,il))

allocate(UTriR(Num,Ntri,il))

allocate(VTriR(Num,Ntri,il))

! primary presurre and velocity
allocate(PressTriRP(Num,Ntri,il))

allocate(UTriRP(Num,Ntri,il))

allocate(VTriRP(Num,Ntri,il))

! first harmonic presurre and velocity
allocate(PressTriRH1(Num,Ntri,il))

allocate(UTriRH1(Num,Ntri,il))

allocate(VTriRH1(Num,Ntri,il))

! second harmonic presurre and velocity
allocate(PressTriRH2(Num,Ntri,il))

allocate(UTriRH2(Num,Ntri,il))

allocate(VTriRH2(Num,Ntri,il))

! Fourier transform of pressure and velocities
allocate(PressTriRFour(Num,Ntri,il))

allocate(UTriRFour(Num,Ntri,il))

allocate(VTriRFour(Num,Ntri,il))


! energy arrays 

allocate(EnergyTriR(Ntri))

allocate(EnergyTriRP(Ntri))

allocate(EnergyTriRH1(Ntri))

allocate(EnergyTriRH2(Ntri))

!allocate(EnergyTotTri(Ntri))

dx=abs(xtrir(1,1)-xtrir(2,1))

dz=abs(ztrir(1,1)-ztrir(2,1))

print*,"dx and dz",dx,dz

print*,zdimr

Ntri=Ntri-1

do i=1,Ntri


!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointTriR(i,:,:),xtrir(i,:),ztrir(i,:),il,nzr,nxr,Nptrir(i),dx,dz)

print*,"Nptrir",Nptrir(i)


!allocate them but it will change in every iteration 
allocate(xi(Nptrir(i)))

allocate(zi(Nptrir(i)))


   do j=1,Nptrir(i)

     xi(j)=x1(PointTriR(i,j,1),PointTriR(i,j,2))

     zi(j)=z1(PointTriR(i,j,1),PointTriR(i,j,2))

  !print*,xi(j),zi(j)

   end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IRTriR(i,:,:),xtrir(i,:),xi,ztrir(i,:),zi,Nptrir(i),il)


  deallocate(xi)

  deallocate(zi)

end do

print*,"okey tri"
!print*,IRTri(1,:,1)

do i=1,Ntri

  do j=1,Num
 

      call interpVal(IRTriR(i,:,:),PressTriR(j,i,:),Press(j,:),PointTriR(i,:,1),PointTriR(i,:,2),nxr,nzr,Nptrir(i),il)
   
      call interpVal(IRTriR(i,:,:),UTriR(j,i,:),u(j,:),PointTriR(i,:,1),PointTriR(i,:,2),nxr,nzr,Nptrir(i),il)

      call interpVal(IRTriR(i,:,:),VTriR(j,i,:),v(j,:),PointTriR(i,:,1),PointTriR(i,:,2),nxr,nzr,Nptrir(i),il)

  
  end do

end do



! the arclength starting from the entry of triangular region 

arctrir=0.


do i=2,Ntri
 
 arctrir(i)=arctrir(i-1)+dsqrt((xtrir(i,il)-xtrir(i-1,il))**2.+(ztrir(i,il)-ztrir(i-1,il))**2.)

end do 


 
parctrir=0.

do j=1,Ntri

 do i=2,il
 
   parctrir(j,i)=parctrir(j,i-1)+dsqrt((xtrir(j,i)-xtrir(j,i-1))**2.+(ztrir(j,i)-ztrir(j,i-1))**2.)

 end do 

end do 

print*,"okey tri 1"



! bandpass filtering in triangular region 

 do i=1,Ntri
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressTriRP(:,i,j),PressTriR(:,i,j),PressTriRFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(UTriRP(:,i,j),UTriR(:,i,j),UTriRFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VTriRP(:,i,j),VTriR(:,i,j),VTriRFour(:,i,j),w0*dt1/nor,bw/nor,Num)
 
   !  First Harmonic
  call BandPassFilter(PressTriRH1(:,i,j),PressTriR(:,i,j),PressTriRFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UTriRH1(:,i,j),UTriR(:,i,j),UTriRFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(VTriRH1(:,i,j),VTriR(:,i,j),VTriRFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressTriRH2(:,i,j),PressTriR(:,i,j),PressTriRFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(UTriRH2(:,i,j),UTriR(:,i,j),UTriRFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(VTriRH2(:,i,j),VTriR(:,i,j),VTriRFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)
 
  end do
 end do

print*, PressTriR(200,1,:)

print*,"okey tri 2"


! output of inlet and outlet profile at t=200 

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(312,file='TotPU1Trir200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriR(200,1,i)*(UTriR(200,1,i)*Cgxtrir(1,i)+VTriR(200,1,i)*Cgztrir(1,i))/dummy1 

  write(312,*) parctrir(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(313,file='PriPU1Trir200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRP(200,1,i)*(UTriRP(200,1,i)*Cgxtrir(1,i)+VTriRP(200,1,i)*Cgztrir(1,i))/dummy1 

  write(313,*) parctrir(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(314,file='1HPU1Trir200.dat',status='unknown')

do i=1,il
dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRH1(200,1,i)*(UTriRH1(200,1,i)*Cgxtrir(1,i)+VTriRH1(200,1,i)*Cgztrir(1,i))/dummy1 

  write(314,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(315,file='2HPU1Trir200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRH2(200,1,i)*(UTriRH2(200,1,i)*Cgxtrir(1,i)+VTriRH2(200,1,i)*Cgztrir(1,i))/dummy1 

  write(315,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


j=Ntri

! Total Energy Flux at the outlet of the reflection path

open(316,file='TotPU1Trir200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriR(200,N(il),i)*(UTriR(200,N(il),i)*Cgxtrir(N(il),i)+VTriR(200,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(316,*)  parctrir(j,i)/lambdax,dummy/E1non
  
end do

! Primary Energy Flux

open(317,file='PriPU1Trir200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRP(200,N(il),i)*(UTriRP(200,N(il),i)*Cgxtrir(N(il),i)+VTriRP(200,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(317,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(318,file='1HPU1Trir200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRH1(200,N(il),i)*(UTriRH1(200,N(il),i)*Cgxtrir(N(il),i)+VTriRH1(200,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(318,*)  parctrir(j,i)/lambdax,dummy/E1non
  
end do


! Second Harmonic Energy Flux 

open(319,file='2HPU1Trir200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRH2(200,N(il),i)*(UTriRH2(200,N(il),i)*Cgxtrir(N(il),i)+VTriRH2(200,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(319,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do







! output of inlet and outlet profile at t=300 

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(322,file='TotPU1Trir300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriR(300,1,i)*(UTriR(300,1,i)*Cgxtrir(1,i)+VTriR(300,1,i)*Cgztrir(1,i))/dummy1 

  write(322,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(323,file='PriPU1Trir300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRP(300,1,i)*(UTriRP(300,1,i)*Cgxtrir(1,i)+VTriRP(300,1,i)*Cgztrir(1,i))/dummy1 

  write(323,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(324,file='1HPU1Trir300.dat',status='unknown')

do i=1,il
dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRH1(300,1,i)*(UTriRH1(300,1,i)*Cgxtrir(1,i)+VTriRH1(300,1,i)*Cgztrir(1,i))/dummy1 

  write(324,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(325,file='2HPU1Trir300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRH2(300,1,i)*(UTriRH2(300,1,i)*Cgxtrir(1,i)+VTriRH2(300,1,i)*Cgztrir(1,i))/dummy1 

  write(325,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


j=Ntri

! Total Energy Flux at the outlet of the reflection path

open(326,file='TotPU1Trir300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriR(300,N(il),i)*(UTriR(300,N(il),i)*Cgxtrir(N(il),i)+VTriR(300,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(326,*)  parctrir(j,i)/lambdax,dummy/E1non
  
end do

! Primary Energy Flux

open(327,file='PriPU1Trir300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRP(300,N(il),i)*(UTriRP(300,N(il),i)*Cgxtrir(N(il),i)+VTriRP(300,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(327,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(328,file='1HPU1Trir300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRH1(300,N(il),i)*(UTriRH1(300,N(il),i)*Cgxtrir(N(il),i)+VTriRH1(300,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(328,*)  parctrir(j,i)/lambdax,dummy/E1non
  
end do


! Second Harmonic Energy Flux 

open(329,file='2HPU1Trir300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRH2(300,N(il),i)*(UTriRH2(300,N(il),i)*Cgxtrir(N(il),i)+VTriRH2(300,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(329,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do





! output of inlet and outlet profile at t=400 

j=1

! Total Energy Flux at the inlet of the reflection path triangular 

open(332,file='TotPU1Trir400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriR(400,1,i)*(UTriR(400,1,i)*Cgxtrir(1,i)+VTriR(400,1,i)*Cgztrir(1,i))/dummy1 

  write(332,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(333,file='PriPU1Trir400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRP(400,1,i)*(UTriRP(400,1,i)*Cgxtrir(1,i)+VTriRP(400,1,i)*Cgztrir(1,i))/dummy1 

  write(333,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(334,file='1HPU1Trir400.dat',status='unknown')

do i=1,il
dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRH1(400,1,i)*(UTriRH1(400,1,i)*Cgxtrir(1,i)+VTriRH1(400,1,i)*Cgztrir(1,i))/dummy1 

  write(334,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(335,file='2HPU1Trir400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(1,i)**2.+Cgztrir(1,i)**2.)
dummy=PressTriRH2(400,1,i)*(UTriRH2(400,1,i)*Cgxtrir(1,i)+VTriRH2(400,1,i)*Cgztrir(1,i))/dummy1 

  write(335,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


j=Ntri

! Total Energy Flux at the outlet of the reflection path

open(336,file='TotPU1Trir400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriR(400,N(il),i)*(UTriR(400,N(il),i)*Cgxtrir(N(il),i)+VTriR(400,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(336,*)  parctrir(j,i)/lambdax,dummy/E1non
  
end do

! Primary Energy Flux

open(337,file='PriPU1Trir400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRP(400,N(il),i)*(UTriRP(400,N(il),i)*Cgxtrir(N(il),i)+VTriRP(400,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(337,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(338,file='1HPU1Trir400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRH1(400,N(il),i)*(UTriRH1(400,N(il),i)*Cgxtrir(N(il),i)+VTriRH1(400,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(338,*)  parctrir(j,i)/lambdax,dummy/E1non
  
end do


! Second Harmonic Energy Flux 

open(339,file='2HPU1Trir400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(Cgxtrir(N(il),i)**2.+Cgztrir(N(il),i)**2.)
dummy=PressTriRH2(400,N(il),i)*(UTriRH2(400,N(il),i)*Cgxtrir(N(il),i)+VTriRH2(400,N(il),i)*Cgztrir(N(il),i))/dummy1 

  write(339,*)  parctrir(j,i)/lambdax,dummy/E1non

  
end do


print*,"okey tri 3"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
!       Energy Integral along Triangular Reflection              !
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Integral for t=200

do i=1,Ntri

 !Total Energy 
 call PressVelocityInteg(EnergyTriR(i),PressTriR(200,i,:),UTriR(200,i,:),VTriR(200,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !Primary Energy 
 call PressVelocityInteg(EnergyTriRP(i),PressTriRP(200,i,:),UTriRP(200,i,:),VTriRP(200,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !First Harmonic Energy 
 call PressVelocityInteg(EnergyTriRH1(i),PressTriRH1(200,i,:),UTriRH1(200,i,:),VTriRH1(200,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !Second Energy 
 call PressVelocityInteg(EnergyTriRH2(i),PressTriRH2(200,i,:),UTriRH2(200,i,:),VTriRH2(200,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 

end do 

open(600,file='EnergyTriR200.dat',status='unknown')

do i=1,Ntri

  write(600,*) arctrir(i)/lambdax,EnergyTriR(i)/Enon
 
end do

open(601,file='PEnergyTriR200.dat',status='unknown')

do i=1,Ntri

  write(601,*) arctrir(i)/lambdax,EnergyTriRP(i)/Enon
 
end do


open(602,file='1HEnergyTriR200.dat',status='unknown')

do i=1,Ntri
  write(602,*) arctrir(i)/lambdax,EnergyTriRH1(i)/Enon
 
end do


open(603,file='2HEnergyTriR200.dat',status='unknown')

do i=1,Ntri

  write(603,*) arctrir(i)/lambdax,EnergyTriRH2(i)/Enon
 
end do


print*,"okey tri 4"

! Integral for t=300

do i=1,Ntri

 !Total Energy 
 call PressVelocityInteg(EnergyTriR(i),PressTriR(300,i,:),UTriR(300,i,:),VTriR(300,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !Primary Energy 
 call PressVelocityInteg(EnergyTriRP(i),PressTriRP(300,i,:),UTriRP(300,i,:),VTriRP(300,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !First Harmonic Energy 
 call PressVelocityInteg(EnergyTriRH1(i),PressTriRH1(300,i,:),UTriRH1(300,i,:),VTriRH1(300,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !Second Energy 
 call PressVelocityInteg(EnergyTriRH2(i),PressTriRH2(300,i,:),UTriRH2(300,i,:),VTriRH2(300,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 

end do 

open(610,file='EnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(610,*) arctrir(i)/lambdax,EnergyTriR(i)/Enon
 
end do

open(611,file='PEnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(611,*) arctrir(i)/lambdax,EnergyTriRP(i)/Enon
 
end do


open(612,file='1HEnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(612,*) arctrir(i)/lambdax,EnergyTriRH1(i)/Enon
 
end do


open(613,file='2HEnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(613,*) arctrir(i)/lambdax,EnergyTriRH2(i)/Enon
 
end do





! Integral for t=400

do i=1,Ntri

 !Total Energy 
 call PressVelocityInteg(EnergyTriR(i),PressTriR(400,i,:),UTriR(400,i,:),VTriR(400,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !Primary Energy 
 call PressVelocityInteg(EnergyTriRP(i),PressTriRP(400,i,:),UTriRP(400,i,:),VTriRP(400,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !First Harmonic Energy 
 call PressVelocityInteg(EnergyTriRH1(i),PressTriRH1(400,i,:),UTriRH1(400,i,:),VTriRH1(400,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 !Second Energy 
 call PressVelocityInteg(EnergyTriRH2(i),PressTriRH2(400,i,:),UTriRH2(400,i,:),VTriRH2(400,i,:),il,Cgxtrir(i,:),Cgztrir(i,:),dr(i))

 

end do 

open(620,file='EnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(620,*) arctrir(i)/lambdax,EnergyTriR(i)/Enon
 
end do

open(621,file='PEnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(621,*) arctrir(i)/lambdax,EnergyTriRP(i)/Enon
 
end do


open(622,file='1HEnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(622,*) arctrir(i)/lambdax,EnergyTriRH1(i)/Enon
 
end do


open(623,file='2HEnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(623,*) arctrir(i)/lambdax,EnergyTriRH2(i)/Enon
 
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                  !
!   Interpolation in Reflection Region after triangular region     !
!                                                                  !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(Npref(N(il)))

allocate(IRRef(N(il),il,nxr*nzr))

allocate(PointRef(N(il),nxr*nzr,2))

! total presurre and velocity
allocate(PressRef(Num,N(il),il))

allocate(URef(Num,N(il),il))

allocate(VRef(Num,N(il),il))

! primary presurre and velocity
allocate(PressRefP(Num,N(il),il))

allocate(URefP(Num,N(il),il))

allocate(VRefP(Num,N(il),il))

! first harmonic presurre and velocity
allocate(PressRefH1(Num,N(il),il))

allocate(URefH1(Num,N(il),il))

allocate(VRefH1(Num,N(il),il))

! second harmonic presurre and velocity
allocate(PressRefH2(Num,N(il),il))

allocate(URefH2(Num,N(il),il))

allocate(VRefH2(Num,N(il),il))

! Fourier transform of pressure and velocities
allocate(PressRefFour(Num,N(il),il))

allocate(URefFour(Num,N(il),il))

allocate(VRefFour(Num,N(il),il))


! energy arrays 

allocate(EnergyRef(N(il)))

allocate(EnergyRefP(N(il)))

allocate(EnergyRefH1(N(il)))

allocate(EnergyRefH2(N(il)))

!allocate(EnergyTotTri(Ntri))

dx=0.25*abs(xref(1,1)-xref(1,2))

dz=0.25*abs(zref(1,1)-zref(1,2))

print*,"dx and dz",dx,dz

do i=1,N(il)


!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointRef(i,:,:),xref(i,:),zref(i,:),il,nzr,nxr,Npref(i),dx,dz)

print*,Npref(i)


!allocate them but it will change in every iteration 
allocate(xi(Npref(i)))

allocate(zi(Npref(i)))


   do j=1,Npref(i)

     xi(j)=x1(PointRef(i,j,1),PointRef(i,j,2))

     zi(j)=z1(PointRef(i,j,1),PointRef(i,j,2))

  !print*,xi(j),zi(j)

   end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IRRef(i,:,:),xref(i,:),xi,zref(i,:),zi,Npref(i),il)


  deallocate(xi)

  deallocate(zi)

end do

print*,"okey tri"
!print*,IRTri(1,:,1)

do i=1,N(il)

  do j=1,Num
 

      call interpVal(IRRef(i,:,:),PressRef(j,i,:),Press(j,:),PointRef(i,:,1),PointRef(i,:,2),nxr,nzr,Npref(i),il)
   
      call interpVal(IRRef(i,:,:),URef(j,i,:),u(j,:),PointRef(i,:,1),PointRef(i,:,2),nxr,nzr,Npref(i),il)

      call interpVal(IRRef(i,:,:),URef(j,i,:),u(j,:),PointRef(i,:,1),PointRef(i,:,2),nxr,nzr,Npref(i),il)

  
  end do

end do

! here we calculate the arclength along the energy flux line 

parcref=0.

do i=2,il
 
 parcref(i)=parcref(i-1)+dsqrt((xref(1,i)-xref(1,i-1))**2.+(zref(1,i)-zref(1,i-1))**2.)

end do 

arcref=0.

do i=2,N(il)
 
 arcref(i)=arcref(i-1)+dsqrt((xref(i,1)-xref(i-1,1))**2.+(zref(i,1)-zref(i-1,1))**2.)

end do 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                             !
!      the band-pass filtering of the interpolated data       !
!                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,N(il)
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressRefP(:,i,j),PressRef(:,i,j),PressRefFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(URefP(:,i,j),URef(:,i,j),URefFour(:,i,j),w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VRefP(:,i,j),VRef(:,i,j),VRefFour(:,i,j),w0*dt1/nor,bw/nor,Num)
 
   !  First Harmonic
  call BandPassFilter(PressRefH1(:,i,j),PressRef(:,i,j),PressRefFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(URefH1(:,i,j),URef(:,i,j),URefFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(VRefH1(:,i,j),VRef(:,i,j),VRefFour(:,i,j),2.*w0*dt1/nor,bw/nor,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressRefH2(:,i,j),PressRef(:,i,j),PressRefFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

  call BandPassFilter(URefH2(:,i,j),URef(:,i,j),URefFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)

   call BandPassFilter(VRefH2(:,i,j),VRef(:,i,j),VRefFour(:,i,j),3.*w0*dt1/nor,bw/nor,Num)
 
  end do
 end do
 

print*,"okey tri 2"

! output of inlet and outlet profile at t=200 

! Total Energy Flux at the inlet of the reflection path

open(412,file='TotPU1Ref200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRef(200,1,i)*(URef(200,1,i)*CgxRef(1,i)+VRef(200,1,i)*CgzRef(1,i))/dummy1 

  write(412,*) parcref(i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(413,file='PriPU1Ref200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefP(200,1,i)*(URefP(200,1,i)*CgxRef(1,i)+VRefP(200,1,i)*CgzRef(1,i))/dummy1 

  write(413,*) parcref(i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(414,file='1HPU1Ref200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefH1(200,1,i)*(URefH1(200,1,i)*CgxRef(1,i)+VRefH1(200,1,i)*CgzRef(1,i))/dummy1 

  write(414,*) parcref(i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(415,file='2HPU1Ref200.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefH2(200,1,i)*(URefH2(200,1,i)*CgxRef(1,i)+VRefH2(200,1,i)*CgzRef(1,i))/dummy1 

  write(415,*) parcref(i)/lambdax,dummy/E1non

  
end do



! Total Energy Flux at the outlet of the reflection path

open(416,file='TotPU1Ref200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRef(200,N(il),i)*(URef(200,N(il),i)*CgxRef(N(il),i)+VRef(200,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(416,*) parcref(i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(417,file='PriPU1Ref200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefP(200,N(il),i)*(URefP(200,N(il),i)*CgxRef(N(il),i)+VRefP(200,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(417,*) parcref(i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(418,file='1HPU1Ref200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefH1(200,N(il),i)*(URefH1(200,N(il),i)*CgxRef(N(il),i)+VRefH1(200,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(418,*) parcref(i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(419,file='2HPU1Ref200out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefH2(200,N(il),i)*(URefH2(200,N(il),i)*CgxRef(N(il),i)+VRefH2(200,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(419,*) parcref(i)/lambdax,dummy/E1non

  
end do









! output of inlet and outlet profile at t=300 

! Total Energy Flux at the inlet of the reflection path

open(402,file='TotPU1Ref300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRef(1,1,i)*(URef(1,1,i)*CgxRef(1,i)+VRef(1,1,i)*CgzRef(1,i))/dummy1 

  write(402,*) parcref(i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(403,file='PriPU1Ref300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefP(1,1,i)*(URefP(1,1,i)*CgxRef(1,i)+VRefP(1,1,i)*CgzRef(1,i))/dummy1 

  write(403,*) parcref(i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(404,file='1HPU1Ref300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefH1(1,1,i)*(URefH1(1,1,i)*CgxRef(1,i)+VRefH1(1,1,i)*CgzRef(1,i))/dummy1 

  write(404,*) parcref(i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(405,file='2HPU1Ref300.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefH2(1,1,i)*(URefH2(1,1,i)*CgxRef(1,i)+VRefH2(1,1,i)*CgzRef(1,i))/dummy1 

  write(405,*) parcref(i)/lambdax,dummy/E1non

  
end do



! Total Energy Flux at the outlet of the reflection path

open(406,file='TotPU1Ref300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRef(1,N(il),i)*(URef(1,N(il),i)*CgxRef(N(il),i)+VRef(1,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(406,*) parcref(i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(407,file='PriPU1Ref300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefP(1,N(il),i)*(URefP(1,N(il),i)*CgxRef(N(il),i)+VRefP(1,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(407,*) parcref(i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(408,file='1HPU1Ref300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefH1(1,N(il),i)*(URefH1(1,N(il),i)*CgxRef(N(il),i)+VRefH1(1,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(408,*) parcref(i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(409,file='2HPU1Ref300out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefH2(1,N(il),i)*(URefH2(1,N(il),i)*CgxRef(N(il),i)+VRefH2(1,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(409,*) parcref(i)/lambdax,dummy/E1non

  
end do





! output of inlet and outlet profile at t=400 

! Total Energy Flux at the inlet of the reflection path

open(422,file='TotPU1Ref400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRef(400,1,i)*(URef(400,1,i)*CgxRef(1,i)+VRef(400,1,i)*CgzRef(1,i))/dummy1 

  write(422,*) parcref(i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(423,file='PriPU1Ref400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefP(400,1,i)*(URefP(400,1,i)*CgxRef(1,i)+VRefP(400,1,i)*CgzRef(1,i))/dummy1 

  write(423,*) parcref(i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(424,file='1HPU1Ref400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefH1(400,1,i)*(URefH1(400,1,i)*CgxRef(1,i)+VRefH1(400,1,i)*CgzRef(1,i))/dummy1 

  write(424,*) parcref(i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(425,file='2HPU1Ref400.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(1,i)**2.+CgzRef(1,i)**2.)
dummy=PressRefH2(400,1,i)*(URefH2(400,1,i)*CgxRef(1,i)+VRefH2(400,1,i)*CgzRef(1,i))/dummy1 

  write(425,*) parcref(i)/lambdax,dummy/E1non

  
end do



! Total Energy Flux at the outlet of the reflection path

open(426,file='TotPU1Ref400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRef(400,N(il),i)*(URef(400,N(il),i)*CgxRef(N(il),i)+VRef(400,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(426,*) parcref(i)/lambdax,dummy/E1non

  
end do

! Primary Energy Flux

open(427,file='PriPU1Ref400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefP(400,N(il),i)*(URefP(400,N(il),i)*CgxRef(N(il),i)+VRefP(400,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(427,*) parcref(i)/lambdax,dummy/E1non

  
end do

! First Harmonic Energy Flux 

open(428,file='1HPU1Ref400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefH1(400,N(il),i)*(URefH1(400,N(il),i)*CgxRef(N(il),i)+VRefH1(400,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(428,*) parcref(i)/lambdax,dummy/E1non

  
end do


! Second Harmonic Energy Flux 

open(429,file='2HPU1Ref400out.dat',status='unknown')

do i=1,il

dummy1=dsqrt(CgxRef(N(il),i)**2.+CgzRef(N(il),i)**2.)
dummy=PressRefH2(400,N(il),i)*(URefH2(400,N(il),i)*CgxRef(N(il),i)+VRefH2(400,N(il),i)*CgzRef(N(il),i))/dummy1 

  write(429,*) parcref(i)/lambdax,dummy/E1non

  
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                !
!            Energy Integral along Reflection Path               !
!                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Integral for t=200

do i=1,N(il)

 !Total Energy 
 call PressVelocityInteg(EnergyRef(i),PressRef(200,i,:),URef(200,i,:),VRef(200,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !Primary Energy 
 call PressVelocityInteg(EnergyRefP(i),PressRefP(200,i,:),URefP(200,i,:),VRefP(200,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !First Harmonic Energy 
 call PressVelocityInteg(EnergyRefH1(i),PressRefH1(200,i,:),URefH1(200,i,:),VRefH1(200,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !Second Energy 
 call PressVelocityInteg(EnergyRefH2(i),PressRefH2(200,i,:),URefH2(200,i,:),VRefH2(200,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 

end do 

open(500,file='EnergyRef200.dat',status='unknown')

do i=1,N(il)

  write(500,*) arcref(i)/lambdax,EnergyRef(i)/Enon
 
end do

open(501,file='PEnergyRef200.dat',status='unknown')

do i=1,N(il)

  write(501,*) arcref(i)/lambdax,EnergyRefP(i)/Enon
 
end do


open(502,file='1HEnergyRef200.dat',status='unknown')

do i=1,N(il)

  write(502,*) arcref(i)/lambdax,EnergyRefH1(i)/Enon
 
end do


open(503,file='2HEnergyRef200.dat',status='unknown')

do i=1,N(il)

  write(503,*) arcref(i)/lambdax,EnergyRefH2(i)/Enon
 
end do


! Integral for t=300

do i=1,N(il)

 !Total Energy 
 call PressVelocityInteg(EnergyRef(i),PressRef(300,i,:),URef(300,i,:),VRef(300,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !Primary Energy 
 call PressVelocityInteg(EnergyRefP(i),PressRefP(300,i,:),URefP(300,i,:),VRefP(300,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !First Harmonic Energy 
 call PressVelocityInteg(EnergyRefH1(i),PressRefH1(300,i,:),URefH1(300,i,:),VRefH1(300,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !Second Energy 
 call PressVelocityInteg(EnergyRefH2(i),PressRefH2(300,i,:),URefH2(300,i,:),VRefH2(300,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 

end do 

open(504,file='EnergyRef300.dat',status='unknown')

do i=1,N(il)

  write(504,*) arcref(i)/lambdax,EnergyRef(i)/Enon
 
end do

open(505,file='PEnergyRef300.dat',status='unknown')

do i=1,N(il)

  write(505,*) arcref(i)/lambdax,EnergyRefP(i)/Enon
 
end do


open(506,file='1HEnergyRef300.dat',status='unknown')

do i=1,N(il)

  write(506,*) arcref(i)/lambdax,EnergyRefH1(i)/Enon
 
end do


open(507,file='2HEnergyRef300.dat',status='unknown')

do i=1,N(il)

  write(507,*) arcref(i)/lambdax,EnergyRefH2(i)/Enon
 
end do



! Integral for t=400

do i=1,N(il)

 !Total Energy 
 call PressVelocityInteg(EnergyRef(i),PressRef(400,i,:),URef(400,i,:),VRef(400,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !Primary Energy 
 call PressVelocityInteg(EnergyRefP(i),PressRefP(400,i,:),URefP(400,i,:),VRefP(400,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !First Harmonic Energy 
 call PressVelocityInteg(EnergyRefH1(i),PressRefH1(400,i,:),URefH1(400,i,:),VRefH1(400,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 !Second Energy 
 call PressVelocityInteg(EnergyRefH2(i),PressRefH2(400,i,:),URefH2(400,i,:),VRefH2(400,i,:),il,CgxRef(i,:),CgzRef(i,:),zdim)

 

end do 

open(508,file='EnergyRef400.dat',status='unknown')

do i=1,N(il)

  write(508,*) arcref(i)/lambdax,EnergyRef(i)/Enon
 
end do

open(509,file='PEnergyRef400.dat',status='unknown')

do i=1,N(il)

  write(509,*) arcref(i)/lambdax,EnergyRefP(i)/Enon
 
end do


open(510,file='1HEnergyRef400.dat',status='unknown')

do i=1,N(il)

  write(510,*) arcref(i)/lambdax,EnergyRefH1(i)/Enon
 
end do


open(511,file='2HEnergyRef400.dat',status='unknown')

do i=1,N(il)

  write(511,*) arcref(i)/lambdax,EnergyRefH2(i)/Enon
 
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                     ! 
!    Averaged Energy Flux over One Wave Period  in triangular reflection region       ! 
!                                                                                     !   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(PUTriRT(Ntri,il))

allocate(PUTriRP(Ntri,il))

allocate(PUTriRZ(Ntri,il))

allocate(PUTriRH1(Ntri,il))

allocate(PUTriRH2(Ntri,il))



! let's generate the averaged pressure and the velocity product

do i=1,Ntri

 do j=1,il
 
  ! we are sending 107 time steps since it is equla to one wave period
  call  MeanPU(PUTriRT(i,j),PressTriR(ti:tf,i,j),UTriR(ti:tf,i,j),VTriR(ti:tf,i,j),Cgxtrir(i,j),Cgztrir(i,j),ttot)

  call  MeanPU(PUTriRP(i,j),PressTriRP(ti:tf,i,j),UTriRP(ti:tf,i,j),VTriRP(ti:tf,i,j),Cgxtrir(i,j),Cgztrir(i,j),ttot)

  call  MeanPU(PUTriRH1(i,j),PressTriRH1(ti:tf,i,j),UTriRH1(ti:tf,i,j),VTriRH1(ti:tf,i,j),Cgxtrir(i,j),Cgztrir(i,j),ttot)
  
  call  MeanPU(PUTriRH2(i,j),PressTriRH2(ti:tf,i,j),UTriRH2(ti:tf,i,j),VTriRH2(ti:tf,i,j),Cgxtrir(i,j),Cgztrir(i,j),ttot)

 end do 

end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!     Total Energy flux a entry region to Triangular Zone       !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


open(932,file='Tot1PUTrirAve.dat',status='unknown')

do i=1,il

  write(932,*) parctrir(1,i)/lambdax,PUTriRT(1,i)/E1non
 
end do



! Total Energy flux at the outlet of the triangular zone 

open(933,file='Tot2PUTrirAve.dat',status='unknown')

do i=1,il

  write(933,*) parctrir(Ntri,i)/lambdax,PUTriRT(Ntri,i)/E1non
 
end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for primary frequency           !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(934,file='P1PUTrirAve.dat',status='unknown')

do i=1,il

  write(934,*) parctrir(1,i)/lambdax,PUTriRP(1,i)/E1non
 
end do



! Primary Energy flux at the outlet of the triangular zone 

open(935,file='P2PUTrirAve.dat',status='unknown')

do i=1,il

  write(935,*) parctrir(Ntri,i)/lambdax,PUTriRP(Ntri,i)/E1non
 
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for first Harmonic              !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(936,file='1Har1PUTrirAve.dat',status='unknown')

do i=1,il

  write(936,*) parctrir(1,i)/lambdax,PUTriRH1(1,i)/E1non
 
end do



! First harmonic Energy flux after one lambdax 

open(938,file='1Har2PUTrirAve.dat',status='unknown')

do i=1,il

  write(938,*) parctrir(Ntri,i)/lambdax,PUTriRH1(Ntri,i)/E1non
 
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                               !
!    Energy flux a entry region for second Harmonic             !
!                                                               !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(939,file='2Har1PUTrirAve.dat',status='unknown')

do i=1,il

  write(939,*) parctrir(1,i)/lambdax,PUTriRH2(1,i)/E1non
 
end do



!Second harmonic Energy flux after one lambdax 

open(940,file='2Har2PUTrirAve.dat',status='unknown')

do i=1,il

  write(940,*) parctrir(Ntri,i)/lambdax,PUTriRH2(Ntri,i)/E1non
 
end do




! Integral of the averaged energy 

do i=1,Ntri

 !Total Energy 
  call RomBerg(PUTriRT(i,:),EnergyTriR(i),dr(i),il)

 ! Primary Frequency Energy 
 call RomBerg(PUTriRP(i,:),EnergyTriRP(i),dr(i),il)

 ! First Harmonic Energy
 call RomBerg(PUTriRH1(i,:),EnergyTriRH1(i),dr(i),il)
 
 ! Second Harmonic Energy
 call RomBerg(PUTriRH2(i,:),EnergyTriRH2(i),dr(i),il)
 
 ! Zeroth Mode Energy
! call RomBerg(PUTriZ(i,:),EnergyTriZ(i),0.5*zdimr(i),il)
 

end do 



open(950,file='EnergyTrirAve.dat',status='unknown')

do i=1,Ntri

  write(950,*) arctrir(i)/lambdax,EnergyTriR(i)/Enon
 
end do


open(951,file='PEnergyTrirAve.dat',status='unknown')

do i=1,Ntri

  write(951,*) arctrir(i)/lambdax,EnergyTriRP(i)/Enon
 
end do


open(952,file='1HEnergyTrirAve.dat',status='unknown')

do i=1,Ntri

  write(952,*) arctrir(i)/lambdax,EnergyTriRH1(i)/Enon
 
end do


open(953,file='2HEnergyTrirAve.dat',status='unknown')

do i=1,Ntri

  write(953,*) arctrir(i)/lambdax,EnergyTriRH2(i)/Enon
 
end do

!open(822,file='ZEnergyPycAve.dat',status='unknown')

!do i=1,Npyc

 ! write(822,*) arcpyc(i)/lambdax,EnergyPycZ(i)/Enon
 
!end do



end program EnergyIntegral
