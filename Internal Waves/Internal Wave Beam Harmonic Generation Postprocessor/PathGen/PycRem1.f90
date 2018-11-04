
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Pycnocline Calculations               !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!this do loop is used to calculate the number of points in pycnocline

ds=0.005


do i=1,1000

if (ds*i > xlen-x(N(il),il)) then

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
 call QuadPoints(Qx,Qz,0.,1.,x(N(il),il),z(N(il),il)+h/2,h/2,il) 

!print*,Qx
!print*,Qz
! now lets find the interpolation points in pycnocline


allocate(xpyc(Npyc,il))
allocate(zpyc(Npyc,il))


! the interpolated data in pycnocline
allocate(PressIntPyc(Num,Npyc,il))

allocate(UIntPyc(Num,Npyc,il))

allocate(VIntPyc(Num,Npyc,il))

! fourier transformed data
allocate(PressIntPycFour(Num,Npyc,il))

allocate(UIntPycFour(Num,Npyc,il))

allocate(VIntPycFour(Num,Npyc,il))

! primary frequency data 
allocate(PressIntPycP(Num,Npyc,il))

allocate(UIntPycP(Num,Npyc,il))

allocate(VIntPycP(Num,Npyc,il))

! first harmonic data
allocate(PressIntPycH(Num,Npyc,il))

allocate(UIntPycH(Num,Npyc,il))

allocate(VIntPycH(Num,Npyc,il))

! second harmonic data
allocate(PressIntPycH2(Num,Npyc,il))

allocate(UIntPycH2(Num,Npyc,il))

allocate(VIntPycH2(Num,Npyc,il))


! group velocities to denote the unit normal in integration
allocate(CgxPyc(Npyc,il))

allocate(CgzPyc(Npyc,il))

! energy in pycnocline
allocate(EnergyPyc(Num,Npyc))

allocate(EnergyPycP(Num,Npyc))

allocate(EnergyPycH(Num,Npyc))

allocate(EnergyPycH2(Num,Npyc))


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
print*,"Okey",Np1(i),dx,dz,h



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
 
do j=1,Num
 
  do i=1,Npyc


   call interpVal(IRPyc(i,:,:),PressIntPyc(j,i,:),Press(j,:),PointPyc(i,:,1),PointPyc(i,:,2),nxr,nzr,Np1(i),il)
 
   call interpVal(IRPyc(i,:,:),UIntPyc(j,i,:),u(j,:),PointPyc(i,:,1),PointPyc(i,:,2),nxr,nzr,Np1(i),il)

   call interpVal(IRPyc(i,:,:),VIntPyc(j,i,:),v(j,:),PointPyc(i,:,1),PointPyc(i,:,2),nxr,nzr,Np1(i),il)

  end do

 end do

 print*,UIntPyc(300,1,:)
! bandpass filtering in pycnocline region 
do i=1,Npyc
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressIntPycP(:,i,j),PressIntPyc(:,i,j),PressIntPycFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(UIntPycP(:,i,j),UIntPyc(:,i,j),UIntPycFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntPycP(:,i,j),VIntPyc(:,i,j),VIntPycFour(:,i,j),w0*dt1/5.,bw/5.,Num)
 
   !  First Harmonic
  call BandPassFilter(PressIntPycH(:,i,j),PressIntPyc(:,i,j),PressIntPycFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(UIntPycH(:,i,j),UIntPyc(:,i,j),UIntPycFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntPycH(:,i,j),VIntPyc(:,i,j),VIntPycFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressIntPycH2(:,i,j),PressIntPyc(:,i,j),PressIntPycFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(UIntPycH2(:,i,j),UIntPyc(:,i,j),UIntPycFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntPycH2(:,i,j),VIntPyc(:,i,j),VIntPycFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

  end do
 end do


print*,"velocity"
print*,UIntPycH2(300,1,:)

do i=1,Npyc

 do j=1,Num 

  
   ! total enegy integration along pycnocline
   call PressVelocityInteg(DumE,EnergyPyc(i,j),PressIntPyc(j,i,:),UIntPyc(j,i,:),VIntPyc(j,i,:),il,0.,1.,h/2.)

   ! Primary Energy integration along pycnocline
   call PressVelocityInteg(DumE,EnergyPycP(i,j),PressIntPycP(j,i,:),UIntPycP(j,i,:),VIntPycP(j,i,:),il,0.,1.,h/2.)

   ! First Harmonic Energy integration along pycnocline
   call PressVelocityInteg(DumE,EnergyPycH(i,j),PressIntPycH(j,i,:),UIntPycH(j,i,:),VIntPycH(j,i,:),il,0.,1.,h/2.)

    ! Second Harmonic Energy integration along pycnocline
   call PressVelocityInteg(DumE,EnergyPycH2(i,j),PressIntPycH2(j,i,:),UIntPycH2(j,i,:),VIntPycH2(j,i,:),il,0.,1.,h/2.)

 end do

end do


open(500,file='TIntEpyc.dat',status='unknown')

do i=1,Npyc

 write(500,*)  (i-1)*ds/lambdax,EnergyPyc(i,Tout)

end do 

open(501,file='PIntEpyc.dat',status='unknown')

do i=1,Npyc

 write(501,*)  (i-1)*ds/lambdax,EnergyPycP(i,Tout)

end do 

open(502,file='HIntEpyc.dat',status='unknown')

do i=1,Npyc

 write(502,*)  (i-1)*ds/lambdax,EnergyPycH(i,Tout)

end do 


open(503,file='H2IntEpyc.dat',status='unknown')

do i=1,Npyc

 write(503,*)  (i-1)*ds/lambdax,EnergyPycH2(i,Tout)

end do 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!   Reflecting Energy From Pycnocline Boundary   !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! The Exact Same lower bound pycnocline points are used 

! those are PressIntPyc(:,1,1),UIntPyc(:,1,1),VIntPyc(:,1,1)

! but this time , rather than U velocity v velocity Should be observed 


open(600,file='RefEnPycT.dat',status='unknown')

do i=1,Npyc

 write(600,*)  (i-1)*ds/lambdax,-1.*PressIntPyc(Tout,1,i)*VIntPyc(Tout,1,i)

end do 

open(601,file='RefEnPycP.dat',status='unknown')

do i=1,Npyc

 write(601,*)  (i-1)*ds/lambdax,-1.*PressIntPycP(Tout,1,i)*VIntPycP(Tout,1,i)

end do 

open(602,file='RefEnPycH.dat',status='unknown')

do i=1,Npyc

 write(602,*)  (i-1)*ds/lambdax,-1.*PressIntPycH(Tout,1,i)*VIntPycH(Tout,1,i)

end do 

open(603,file='RefEnPycH.dat',status='unknown')

do i=1,Npyc

 write(603,*)  (i-1)*ds/lambdax,-1.*PressIntPycH2(Tout,1,i)*VIntPycH2(Tout,1,i)

end do 
