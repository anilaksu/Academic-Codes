!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Reflected Beam Calculations           !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(xref(N(1),il))
allocate(zref(N(1),il))


allocate(CgxRef(N(1),il))
allocate(CgzRef(N(1),il))

allocate(IRRef(N(1),il,nxr*nzr))

!the interpolated data
allocate(PressIntRefFour(Num,N(1),il))

allocate(UIntRefFour(Num,N(1),il))

allocate(VIntRefFour(Num,N(1),il))

!the interpolated data
allocate(PressIntRef(Num,N(1),il))

allocate(UIntRef(Num,N(1),il))

allocate(VIntRef(Num,N(1),il))

! the primary frequency data
allocate(PressIntRefP(Num,N(1),il))

allocate(UIntRefP(Num,N(1),il))

allocate(VIntRefP(Num,N(1),il))

! first harmonic data
allocate(PressIntRefH(Num,N(1),il))

allocate(UIntRefH(Num,N(1),il))

allocate(VIntRefH(Num,N(1),il))

allocate(PointRef(N(1),nxr*nzr,2))


! energy in pycnocline
allocate(EnergyRef(Num,N(1)))

allocate(EnergyRefP(Num,N(1)))

allocate(EnergyRefH(Num,N(1)))



! the dimension of the point array is N(1) as it is the largest 
allocate(Np2(N(1)))

! the width of the reflected beam is determined as the difference between x-coordinate of the IWB in to pycnocline


! zdimr =( x(N(1),1)-x(N(5),5))/2.

 xcenr =( x(N(1),1)+x(N(il),il))/2.

 zcenr = z(N(il),il)

 ! calculation of the reflected beam path 

 ! it is not generated using path command as it does not differentiate positive and negative direction  


 do j=1 ,il

  do i=1 ,N(il)




 
  xref(i,j)=2.*xcenr-x(i,j)

  zref(i,j)=z(i,j)

  ! the group velocities are calculated here it shares the same magnitude with IWB but it propogates in negative z direction

  CgzRef(i,j)=-1.*Cgz(i,j)
 
  CgxRef(i,j)=Cgx(i,j)
 
 

  end do 

 end do 


! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once


do i=1,N(il)

!to determine the distance in which the interpolation data to be extracted 

dz=2.*abs(z(2,1)-z(1,1))

dx=2.*abs(x(2,1)-x(1,1))

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointRef(i,:,:),xref(i,:),zref(i,:),il,nzr,nxr,Np2(i),dx,dz)

print*,"Energy Ref",Np2(i)

!allocate them but it will change in every iteration 
allocate(xi(Np2(i)))

allocate(zi(Np2(i)))


do j=1,Np2(i)

xi(j)=x1(PointRef(i,j,1),PointRef(i,j,2))

zi(j)=z1(PointRef(i,j,1),PointRef(i,j,2))


end do 

 call interpolation2D(IRRef(i,:,:),xref(i,:),xi,zref(i,:),zi,Np2(i),5)

deallocate(xi)

deallocate(zi)

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!    Reflected Pressure ,Velocity,Energy         !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 


do j=1,Num
 
  do i=1,N(il)

    call interpVal(IRRef(i,:,:),PressIntRef(j,i,:),Press(j,:),PointRef(i,:,:1),PointRef(i,:,2),nxr,nzr,Np2(i),il)
 
    call interpVal(IRRef(i,:,:),UIntRef(j,i,:),u(j,:),PointRef(i,:,:1),PointRef(i,:,2),nxr,nzr,Np2(i),il)

    call interpVal(IRRef(i,:,:),VIntRef(j,i,:),v(j,:),PointRef(i,:,:1),PointRef(i,:,2),nxr,nzr,Np2(i),il)

  end do

 end do


! bandpass filtering in reflection region 
do i=1,N(il)
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressIntRefP(:,i,j),PressIntRef(:,i,j),PressIntRefFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(UIntRefP(:,i,j),UIntRef(:,i,j),UIntRefFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntRefP(:,i,j),VIntRef(:,i,j),VIntRefFour(:,i,j),w0*dt1/5.,bw/5.,Num)
 
   !  First Harmonic
  call BandPassFilter(PressIntRefH(:,i,j),PressIntRef(:,i,j),PressIntRefFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(UIntRefH(:,i,j),UIntRef(:,i,j),UIntRefFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntRefH(:,i,j),VIntRef(:,i,j),VIntRefFour(:,i,j),2.*w0*dt/5.,bw/5.,Num)
 

  end do
 end do

print*,"Okey Ref"


! in reflected integration , direction should be corrected  
kx=-1.*kx

do i=1,N(il)

 do j=1,Num 

  
   ! total enegy integration along IWB 
   call PressVelocityInteg(DumE,EnergyRef(i,j),PressIntRef(j,i,:),UIntRef(j,i,:),VIntRef(j,i,:),il,kx,kz,zdim)

   ! Primary Energy integration along IWB
   call PressVelocityInteg(DumE,EnergyRefP(i,j),PressIntRefP(j,i,:),UIntRefP(j,i,:),VIntRefP(j,i,:),il,kx,kz,zdim)

   ! First Harmonic Energy integration along IWB
   call PressVelocityInteg(DumE,EnergyRefH(i,j),PressIntRefH(j,i,:),UIntRefH(j,i,:),VIntRefH(j,i,:),il,kx,kz,zdim)

 end do

end do

! again let's correct ds length 

ds=0.1


open(600,file='TIntERef.dat',status='unknown')

do i=N(il),1,-1

 write(600,*)  (i-1)*ds/lambdax,EnergyRef(i,300)

end do 

open(601,file='PIntERef.dat',status='unknown')

do i=N(il),1,-1

 write(601,*)  (i-1)*ds/lambdax,EnergyRefP(i,300)

end do 

open(602,file='HIntERef.dat',status='unknown')

do i=N(il),1,-1

 write(602,*)  (i-1)*ds/lambdax,EnergyRefH(i,300)

end do 

