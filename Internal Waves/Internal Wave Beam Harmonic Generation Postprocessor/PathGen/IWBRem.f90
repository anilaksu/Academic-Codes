dz=abs(Qz(2)-Qz(1))

dx=abs(Qx(2)-Qx(1))


print*,"okey 4"
do i=1,N(il)

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,Point(i,:,:),x(i,:),z(i,:),il,nzr,nxr,Np(i),dx,dz)

!print*,Np(1)


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
 

      call interpVal(IR(i,:,:),PressInt(j,i,:),Press(j,:),Point(i,:,1),Point(i,:,2),nxr,nzr,Np(i),il)
   
      call interpVal(IR(i,:,:),UInt(j,i,:),u(j,:),Point(i,:,1),Point(i,:,2),nxr,nzr,Np(i),il)

      call interpVal(IR(i,:,:),VInt(j,i,:),v(j,:),Point(i,:,1),Point(i,:,2),nxr,nzr,Np(i),il)

  
  end do

end do


!print*,"Pressure"
!print*,PressInt(300,1,:)
!print*,"u velocity"
!print*,UInt(300,1,:)
!print*,"v velocity"
!print*,VInt(300,1,:)

print*,"okey"


! bandpass filtering in incident wave beam region 

 do i=1,N(il)
   do j=1,il

    !Primary Wave
   call BandPassFilter(PressIntP(:,i,j),PressInt(:,i,j),PressIntFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(UIntP(:,i,j),UInt(:,i,j),UIntFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntP(:,i,j),VInt(:,i,j),VIntFour(:,i,j),w0*dt1/5.,bw/5.,Num)
 
   !  First Harmonic
  call BandPassFilter(PressIntH(:,i,j),PressInt(:,i,j),PressIntFour(:,i,j),2.*w0*dt1/25.,bw/25.,Num)

  call BandPassFilter(UIntH(:,i,j),UInt(:,i,j),UIntFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntH(:,i,j),VInt(:,i,j),VIntFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)
 

  end do
 end do

!print*,"Pressure"
!print*,PressInt(300,1,:)
!print*,"u velocity"
!print*,UInt(300,1,:)
!print*,"v velocity"
!print*,VInt(300,1,:)

do i=1,N(il)

 do j=1,Num 

  
   ! total enegy integration along IWB 
   call PressVelocityInteg(DumE,Energy(i,j),PressInt(j,i,:),UInt(j,i,:),VInt(j,i,:),il,kx,kz,zdim)

   ! Primary Energy integration along IWB
   call PressVelocityInteg(DumE,EnergyP(i,j),PressIntP(j,i,:),UIntP(j,i,:),VIntP(j,i,:),il,kx,kz,zdim)

   ! First Harmonic Energy integration along IWB
   call PressVelocityInteg(DumE,EnergyH(i,j),PressIntH(j,i,:),UIntH(j,i,:),VIntH(j,i,:),il,kx,kz,zdim)

 end do

end do

open(404,file='TPress300.dat',status='unknown')

do j=1,il

 write(404,*) ((j-1)*zdim/(il-1))/lambdax,PressInt(300,1,j)

end do 
print*,"okey"


open(405,file='TEnergyFluxIWB.dat',status='unknown')

do j=1,il

 E1=PressInt(300,1,j)*(dcos(pi/4.)*UInt(300,1,j)+dcos(pi/4.)*VInt(300,1,j))/Enon

 write(405,*) (j*zdim/(il-1))/lambdax,E1

end do 

open(406,file='PEnergyFluxIWB.dat',status='unknown')

do j=1,il

 E1=PressIntP(300,1,j)*(dcos(pi/4.)*UIntP(300,1,j)+dcos(pi/4.)*VIntP(300,1,j))/Enon

 write(406,*)  (j*zdim/(il-1))/lambdax,E1

end do 

open(407,file='HEnergyFluxIWB.dat',status='unknown')

do j=1,il

 E1=PressIntH(300,1,j)*(dcos(pi/4.)*UIntH(300,1,j)+dcos(pi/4.)*VIntH(300,1,j))/Enon

 write(407,*)  (j*zdim/(il-1))/lambdax,E1

end do 

open(400,file='TIntEIWB.dat',status='unknown')

do i=1,N(il)

 write(400,*)  (i-1)*ds/lambdax,Energy(i,300)

end do 

open(401,file='PIntEIWB.dat',status='unknown')

do i=1,N(il)

 write(401,*)  (i-1)*ds/lambdax,EnergyP(i,300)

end do 

open(402,file='HIntEIWB.dat',status='unknown')

do i=1,N(il)

 write(402,*)  (i-1)*ds/lambdax,EnergyH(i,300)

end do 

! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once



!to determine the distance in which the interpolation data to be extracted 



