program EnergyIntegral

    implicit none

  !p1,p2 are the pointers for the nodal distance btwn the center of IWB and pycnocline in z axis
  integer i,j,k,p1,p2,nxr,nzr,Npyc,Nref,Ntri,Ntrih,Ntrih2,Npri,Nh1,Nh2,RefP,RefH1,RefH2
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
  real*8 ,allocatable::Cgx(:,:),Cgz(:,:),Cgxpri(:,:),Cgzpri(:,:)
  !the group velocities
  real*8 ,allocatable::Cgxh1(:,:),Cgzh1(:,:),Cgxh2(:,:),Cgzh2(:,:)
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
   !the frequency of the wave and band with for bandpass filtering
  real*8 w0 ,bw
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

   ! The array used to keep the interpolated velocities in pycnocline
   real*8, allocatable::UIntPycP(:,:,:),VIntPycP(:,:,:),UIntPycH(:,:,:),VIntPycH(:,:,:),UIntPycH2(:,:,:),VIntPycH2(:,:,:)
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
   ! the ray paths in pycnocline
   real*8,allocatable::xpri(:,:),zpri(:,:),xh1(:,:),zh1(:,:),xh2(:,:),zh2(:,:)
   
  
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

  N1 = 3.9

  Enon=rho*(lambdax**2.)*(u0**2.)*N1*10E-7

  E1non=Enon/lambdax

  Pnon=rho*lambdax*u0*N1*10E-7

  eps=10E-4
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  nxr=512
  
  nzr=491

   ds=0.1

   xcen=0.2

   zcen =0.13

   !r=2.
    !r=10.
   !dt=0.025
    r=6. 
    dt= 0.0270295880980098
   ! r=2.
   ! 0.0344188927666368 
   ! r=4.
   ! 0.0287606935364547
  
   ! r=8.
   ! dt=0.022

   z0 = 0.2625

   pi=4.*datan(1.d0)
   
   kx=-2*pi/0.0875

   kz=-1.*kx

   zlen=0.4

 

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
   h=-2.*arg*dt*0.7/2

   print*,"pycnocline thickness",h/lambdax

allocate(Qx(il))
allocate(Qz(il))
allocate(N(il))



xlen=0.7

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
ds=0.00001

do i=1,Ntri

 zdimr(i)=0.

 !x(N(il)+i,1)

 zdum=z(N(il)+i,1)


 ! in this loop , we find the length of the integration path 
  do j=1,2000 
 
   zdum=zdum+ds

   call BVprofile(BV,N1,r,zdum,z0,dt)

    if( BV<N1/1.412) then 
  ! zdum>z(N(il),il)+ 2.*h .or.
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


 end do

print*,"starting point",x(N(il),1),z(N(il),1)

! IWB lines

open(201,file='IWBPath1.dat',status='unknown')

do i=1,N(1)

  write(201,*) x(i,1),z(i,1)
  
end do


open(202,file='IWBPath5.dat',status='unknown')

do i=1,N(5)

  write(202,*) x(i,5),z(i,5)
  
end do


open(203,file='IWBPath9.dat',status='unknown')

do i=1,N(9)

  write(203,*) x(i,9),z(i,9)
  
end do


open(204,file='IWBPath13.dat',status='unknown')

do i=1,N(13)

  write(204,*) x(i,13),z(i,13)
  
end do


open(205,file='IWBPath17.dat',status='unknown')

do i=1,N(17)

  write(205,*) x(i,17),z(i,17)
  
end do


! The trianglar area isolines

open(310,file='TriPath1.dat',status='unknown')

do i=1,il

  write(310,*) xtri(1,i),ztri(1,i)
  !print*,N(i)
end do

open(311,file='TriPath20.dat',status='unknown')

do i=1,il

  write(311,*) xtri(20,i),ztri(20,i)
  !print*,N(i)
end do

open(312,file='TriPath40.dat',status='unknown')

do i=1,il

  write(312,*) xtri(40,i),ztri(40,i)
  !print*,N(i)
end do

open(313,file='TriPath60.dat',status='unknown')

do i=1,il

  write(313,*) xtri(60,i),ztri(60,i)
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


open(400,file='TriPathR1.dat',status='unknown')

do i=1,il

  write(400,*) xtrir(1,i),ztrir(1,i)
  !print*,N(i)
end do

open(401,file='TriPathR20.dat',status='unknown')

do i=1,il

  write(401,*) xtrir(20,i),ztrir(20,i)
  !print*,N(i)
end do

open(402,file='TriPathR40.dat',status='unknown')

do i=1,il

  write(402,*) xtrir(40,i),ztrir(40,i)
  !print*,N(i)
end do

open(403,file='TriPathR60.dat',status='unknown')

do i=1,il

  write(403,*) xtrir(60,i),ztrir(60,i)
  !print*,N(i)
end do

!RefLection lines

open(501,file='RefPath1.dat',status='unknown')

do i=1,N(il)

  write(501,*) xref(i,1),zref(i,1)
  
end do


open(502,file='RefPath5.dat',status='unknown')

do i=1,N(il)

  write(502,*) xref(i,5),zref(i,5)
  
end do


open(503,file='RefPath9.dat',status='unknown')

do i=1,N(il)

  write(503,*) xref(i,9),zref(i,9)
  
end do


open(504,file='RefPath13.dat',status='unknown')

do i=1,N(il)

  write(504,*) xref(i,13),zref(i,13)
  
end do


open(505,file='RefPath17.dat',status='unknown')

do i=1,N(il)

  write(505,*) xref(i,17),zref(i,17)
  
end do



print*,"Primary and harmonic parameters"

print*,"wave numbers",kx

print*,"frequencies",w0,Ntri

print*,"z",z(N(1),1)

print*,"h",h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Pycnocline Calculations               !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!this do loop is used to calculate the number of points in pycnocline

ds=0.005


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

print*,"okey pyc"

!print*,Qx
!print*,Qz
! now lets find the interpolation points in pycnocline


allocate(xpyc(Npyc,il))
allocate(zpyc(Npyc,il))






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


!RefLection lines

open(601,file='PycPath1.dat',status='unknown')

do i=1,Npyc

  write(601,*) xpyc(i,1),zpyc(i,1)
  
end do


open(602,file='PycPath5.dat',status='unknown')

do i=1,Npyc

  write(602,*) xpyc(i,5),zpyc(i,5)
  
end do


open(603,file='PycPath9.dat',status='unknown')

do i=1,Npyc

  write(603,*) xpyc(i,9),zpyc(i,9)
  
end do


open(604,file='PycPath13.dat',status='unknown')

do i=1,Npyc

  write(604,*) xpyc(i,13),zpyc(i,13)
  
end do


open(605,file='PycPath17.dat',status='unknown')

do i=1,Npyc

  write(605,*) xpyc(i,17),zpyc(i,17)
  
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                             !
!     the reflection paths in pycnocline      !
!                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! let's find the number of the points in pycnocline 

ds=0.00001

! the number of the reflections 

RefP=4


! path for primary frequency
 
 call PathPoints(Npri,z0,z(N(1),1),nzr,kx,r,dt,ds,N1,w0)

print*,"NumPoint",Npri

allocate(xpri(Npri*RefP,il))
allocate(zpri(Npri*RefP,il))


allocate(Cgxpri(Npri*RefP,il))
allocate(Cgzpri(Npri*RefP,il))



 do i=1,il

  call RayPath(Npri,z0,z(N(1),1),x(N(i),i),kx,r,dt,ds,N1,w0,xpri(:,i),zpri(:,i),RefP,Cgxpri(:,i),Cgzpri(:,i))

 end do 



!print*,"z",z(N(1),1)

open(610,file='RayPri1.dat',status='unknown')

do i=1,Npri*RefP

 write(610,*) (xpri(i,1)-x(N(il),17))/lambdax,zpri(i,1)/lambdax

end do


open(611,file='RayPri5.dat',status='unknown')

do i=1,Npri*RefP

 write(611,*) (xpri(i,5)-x(N(il),17))/lambdax,zpri(i,5)/lambdax

end do


open(612,file='RayPri9.dat',status='unknown')

do i=1,Npri*RefP

 write(612,*) (xpri(i,9)-x(N(il),17))/lambdax,zpri(i,9)/lambdax

end do


open(613,file='RayPri13.dat',status='unknown')

do i=1,Npri*RefP

 write(613,*) (xpri(i,13)-x(N(il),17))/lambdax,zpri(i,13)/lambdax

end do

open(614,file='RayPri17.dat',status='unknown')

do i=1,Npri*RefP

 write(614,*) (xpri(i,il)-x(N(il),17))/lambdax,zpri(i,il)/lambdax

end do
! path for first harmonic




end program EnergyIntegral
