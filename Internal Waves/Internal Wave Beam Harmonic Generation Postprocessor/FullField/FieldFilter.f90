program FieldFilter

    implicit none

  !p1,p2 are the pointers for the nodal distance btwn the center of IWB and pycnocline in z axis
  integer i,j,k,p1,p2,nxr,nzr,Npyc,Nref,t1,coun,ti,tf,ttot
  ! the number of the data files going to be used
  integer NumDat
  ! the number of time steps to be used 
  integer Num
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
  real*8 , allocatable::x1(:,:),z1(:,:),brunt(:)
  ! ds the time integration difference 
  real*8 ds ,dx,dz,dl
   !Nt dummy BV profile
  real*8 Nt
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
   complex , allocatable ::UFour(:,:),VFour(:,:)
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   nxr=768
  
  nzr=533


   ds=0.1

   xcen=0.2

   zcen =0.13


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! initial and final time step for time averaging! 
     ti=200
     tf=316
   ! and total time step tott=tf-ti+1
     ttot=tf-ti+1
     nor=10.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   pi=4.*datan(1.d0)
   
   kx=-2*pi/0.0875

   kz=-1.*kx

   zlen=0.4

  ! the number of the pressure and velocity data files
  NumDat=200

  Num=512

  allocate(time(NumDat))

   allocate(Nums(NumDat))

  allocate(Pfiles(NumDat))

  allocate(Vfiles(NumDat))




allocate(x1(nxr,nzr))
allocate(z1(nxr,nzr))  

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

  


allocate(Press(j,nxr*nzr))

allocate(u(j,nxr*nzr))

allocate(v(j,nxr*nzr))


!print*,"NumDat"
!print*,j

 !  do k=1,j
    
  ! print*,k,time(k),Nums(k)

  ! end do

allocate(PressD(Num,nxr*nzr))

allocate(uD(Num,nxr*nzr))

allocate(vD(Num,nxr*nzr))

! primary , harmonic and frequency domain pressure 

allocate(PressH2(Num,nxr*nzr))

allocate(PressH(Num,nxr*nzr))

allocate(PressP(Num,nxr*nzr))

allocate(PressZ(Num,nxr*nzr))

allocate(PressFour(Num,nxr*nzr))

! primary , harmonic and frequency domain u velocity 

allocate(UH2(Num,nxr*nzr))

allocate(UH(Num,nxr*nzr))

allocate(UP(Num,nxr*nzr))

allocate(UZ(Num,nxr*nzr))

allocate(UFour(Num,nxr*nzr))

! primary , harmonic and frequency domain u velocity 

allocate(VH2(Num,nxr*nzr))

allocate(VH(Num,nxr*nzr))

allocate(VP(Num,nxr*nzr))

allocate(VZ(Num,nxr*nzr))

allocate(VFour(Num,nxr*nzr))


! pressure velocity products and their L2 norms

allocate(PU(Num,nxr*nzr))

allocate(PV(Num,nxr*nzr))

allocate(PVelAbs(Num,nxr*nzr))

! pressure velocity products and their L2 norms of primary wave

allocate(PUP(Num,nxr*nzr))

allocate(PVP(Num,nxr*nzr))

allocate(PVelAbsP(Num,nxr*nzr))

! pressure velocity products and their L2 norms of primary wave

allocate(PUH(Num,nxr*nzr))

allocate(PVH(Num,nxr*nzr))

allocate(PVelAbsH(Num,nxr*nzr))


!Mean Quantities

allocate(PUMean(nxr,nzr))
allocate(PVMean(nxr,nzr))  

allocate(PUZMean(nxr,nzr))
allocate(PVZMean(nxr,nzr))  

allocate(PUPMean(nxr,nzr))
allocate(PVPMean(nxr,nzr))  

allocate(PUH1Mean(nxr,nzr))
allocate(PVH1Mean(nxr,nzr))  

allocate(PUH2Mean(nxr,nzr))
allocate(PVH2Mean(nxr,nzr))  

allocate(U2Mean(nxr,nzr))
allocate(V2Mean(nxr,nzr))
allocate(EMean(nxr,nzr))


   ! here we generate the file names for Pressure
    call pread(Pfiles,NumDat)

   ! here we generate the file names for Velocity
    call vread(Vfiles,NumDat)
   
   ! here I read all data from all files
 
    do i=1,j

    print*,Nums(i)
    call DataRead(pfiles(Nums(i)),vfiles(Nums(i)),Press(i,:),u(i,:),v(i,:),nxr,nzr,x1,z1)

    end do
    print*,v(1,35*nxr+250)
  
    
    ! that routine interpolates data in time with uniform time step size
 
    call TimeInter(PressD,uD,vD,Press,u,v,time,j,Num,nxr,nzr,dt1)
       
    print*,vD(1,35*nxr+250)

  


  w0=-1.*N1*kx/dsqrt(kx**2.+kz**2.)
  
  !dt1=2.*pi/(w0*8.)

  ! dt1=0.1
 
  bw=w0*dt1/2.

   print*,"w,bw,dt1"
   print*,w0,bw,dt1
print*,"okey0"

 open(20,file='Unfiltered.dat',status='unknown')

   do i=1,Num
  
    write(20,*) i,vD(i,35*nxr+250)   

   end do

 
   call BandPassFilter(VP(:,35*nxr+250),vD(:,35*nxr+250),VFour(:,35*nxr+250),w0*dt1/nor,bw/nor,Num)

 open(21,file='filtered.dat',status='unknown')

   do i=1,Num
  
    write(21,*) i,vP(i,35*nxr+250)   

   end do

  open(22,file='fourier.dat',status='unknown')

   do i=1,Num
  
    write(22,*) i,abs(VFour(i,35*nxr+250))   

   end do



do i=1,nxr*nzr
 
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
 
do j=1,nzr

  do i=1,nxr

!  write(19,*) x1(i,j)/lambdax,z1(i,j)/lambdax,abs(PressFour(t1,(j-1)*nxr+i))

   write(20,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressD(t1,(j-1)*nxr+i)/Pnon
  
   write(21,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressP(t1,(j-1)*nxr+i)/Pnon

   write(22,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressH(t1,(j-1)*nxr+i)/Pnon
  
  end do 

 end do


!U velocity output for entire field total primary and harmonic
open(23,file='TUt100.dat',status='unknown')
  
 open(24,file='PUt100.dat',status='unknown') 

   open(25,file='HUt100.dat',status='unknown') 
 
do j=1,nzr

  do i=1,nxr

   write(23,*) x1(i,j)/lambdax,z1(i,j)/lambdax,uD(t1,(j-1)*nxr+i)/u0
  
   write(24,*) x1(i,j)/lambdax,z1(i,j)/lambdax,UP(t1,(j-1)*nxr+i)/u0

   write(25,*) x1(i,j)/lambdax,z1(i,j)/lambdax,UH(t1,(j-1)*nxr+i)/u0
  
  end do 

 end do

!V velocity output for entire field total primary and harmonic
open(26,file='TVt100.dat',status='unknown')
  
 open(27,file='PVt100.dat',status='unknown') 

   open(28,file='HVt100.dat',status='unknown') 
 
do j=1,nzr

  do i=1,nxr

   write(26,*) x1(i,j)/lambdax,z1(i,j)/lambdax,vD(t1,(j-1)*nxr+i)/u0
  
   write(27,*) x1(i,j)/lambdax,z1(i,j)/lambdax,VP(t1,(j-1)*nxr+i)/u0

   write(28,*) x1(i,j)/lambdax,z1(i,j)/lambdax,VH(t1,(j-1)*nxr+i)/u0
  
  end do 

end do
!PU product for total primary and harmonic
open(29,file='TPUt100.dat',status='unknown')
  
 open(30,file='PPUt100.dat',status='unknown') 

   open(31,file='HPUt100.dat',status='unknown') 
 
do j=1,nzr

  do i=1,nxr

 !  write(29,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressD(t1,(j-1)*nxr+i)*uD(t1,(j-1)*nxr+i)/Enon
  
 !  write(30,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressP(t1,(j-1)*nxr+i)*UP(t1,(j-1)*nxr+i)/Enon

 !  write(31,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressH(t1,(j-1)*nxr+i)*UH(t1,(j-1)*nxr+i)/Enon
  
  end do 

 end do

!PV product for total primary and harmonic
open(32,file='TPVt100.dat',status='unknown')
  
 open(33,file='PPVt100.dat',status='unknown') 

   open(34,file='HPVt100.dat',status='unknown') 
 
do j=1,nzr

  do i=1,nxr

 !  write(32,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressD(t1,(j-1)*nxr+i)*vD(t1,(j-1)*nxr+i)/Enon
  
 !  write(33,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressP(t1,(j-1)*nxr+i)*VP(t1,(j-1)*nxr+i)/Enon

 !  write(34,*) x1(i,j)/lambdax,z1(i,j)/lambdax,PressH(t1,(j-1)*nxr+i)*VH(t1,(j-1)*nxr+i)/Enon
  
  end do 

 end do



 
print*,"okey1"
  
 do j=1,nzr

  do i=1,nxr
 
  ! we are sending 107 time steps since it is equla to one wave period
   
  k=(j-1)*nxr+i
  
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


 do j=1,nzr

  do i=1,nxr
 
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

