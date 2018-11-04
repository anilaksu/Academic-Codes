!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!   Incident Beam Pressure ,Velocity,Energy      !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

print*,"Okey"

do j=1,Num
 
  do i=1,2


    call interpVal(IR(i,:,:),PressInt(j,i,:),Press(j,:),Point(i,:,:1),Point(i,:,2),nxr,nzr,Np(i),5)
 
    call interpVal(IR(i,:,:),UInt(j,i,:),u(j,:),Point(i,:,:1),Point(i,:,2),nxr,nzr,Np(i),5)

    call interpVal(IR(i,:,:),VInt(j,i,:),v(j,:),Point(i,:,:1),Point(i,:,2),nxr,nzr,Np(i),5)

    !print*,UInt(i,:)
  end do

 end do
 



! bandpass filtering in incident wave beam region 

 do i=1,2
   do j=1,5

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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Pycnocline Calculations               !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!this do loop is used to calculate the number of points in pycnocline

ds=0.01
do i=1,100

if (ds*i > xlen-x(N(5),5)) then

Npyc=i-1

goto 500
end if




end do

500 continue



! to calculate the quadrature points in pycnocline 

! we start from x(N(5),5) to reduce the computational effort 

! And also our beam enters to pycnocline at that point 
 call QuadPoints(Qx,Qz,1.,0.,x(N(5),5),z(N(5),5)+h/2,h/2,5) 

!print*,Qx
!print*,Qz
! now lets find the interpolation points in pycnocline


allocate(xpyc(Npyc,5))
allocate(zpyc(Npyc,5))


! the interpolated data in pycnocline
allocate(PressIntPyc(Num,Npyc,5))

allocate(UIntPyc(Num,Npyc,5))

allocate(VIntPyc(Num,Npyc,5))

! fourier transformed data
allocate(PressIntPycFour(Num,Npyc,5))

allocate(UIntPycFour(Num,Npyc,5))

allocate(VIntPycFour(Num,Npyc,5))

! primary frequency data 
allocate(PressIntPycP(Num,Npyc,5))

allocate(UIntPycP(Num,Npyc,5))

allocate(VIntPycP(Num,Npyc,5))

! first harmonic data
allocate(PressIntPycH(Num,Npyc,5))

allocate(UIntPycH(Num,Npyc,5))

allocate(VIntPycH(Num,Npyc,5))

! group velocities to denote the unit normal in integration
allocate(CgxPyc(Npyc,5))

allocate(CgzPyc(Npyc,5))

! energy in pycnocline
allocate(EnergyPyc(Num,Npyc))

allocate(EnergyPycP(Num,Npyc))

allocate(EnergyPycH(Num,Npyc))



 CgxPyc=1.
 CgzPyc=0.

do i=1,Npyc 

do j=1,5
! ds increment in x 
xpyc(i,j)=(i-1)*ds+Qx(j)

!z is same along the path 

zpyc(i,j)=Qz(j)

end do


end do




allocate(Np1(Npyc))

allocate(IRPyc(Npyc,5,nxr*nzr))

allocate(PointPyc(Npyc,nxr*nzr,2))



! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once


do i=1,Npyc

dz=4.*abs(Qz(2)-Qz(1))

dx=dz

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointPyc(i,:,:),xpyc(i,:),zpyc(i,:),5,nzr,nxr,Np1(i),dx,dz)
print*,"Okey",Np1(i),dx,dz,h



!allocate them but it will change in every iteration 
allocate(xi(Np1(i)))

allocate(zi(Np1(i)))


do j=1,Np1(i)

 xi(j)=x1(PointPyc(i,j,1),PointPyc(i,j,2))

 zi(j)=z1(PointPyc(i,j,1),PointPyc(i,j,2))

 
end do 

 call interpolation2D(IRPyc(i,:,:),xpyc(i,:),xi,zpyc(i,:),zi,Np1(i),5)

 
deallocate(xi)

deallocate(zi)



end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!   Pycnocline Pressure ,Velocity,Energy         !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
do j=1,Num
 
  do i=6,8


   call interpVal(IRPyc(i,:,:),PressIntPyc(j,i,:),Press(j,:),PointPyc(i,:,:1),PointPyc(i,:,2),nxr,nzr,Np1(i),5)
 
   call interpVal(IRPyc(i,:,:),UIntPyc(j,i,:),u(j,:),PointPyc(i,:,:1),PointPyc(i,:,2),nxr,nzr,Np1(i),5)

   call interpVal(IRPyc(i,:,:),VIntPyc(j,i,:),v(j,:),PointPyc(i,:,:1),PointPyc(i,:,2),nxr,nzr,Np1(i),5)

  end do

 end do
 

! bandpass filtering in pycnocline region 
do i=6,8
   do j=1,5

    !Primary Wave
   call BandPassFilter(PressIntPycP(:,i,j),PressIntPyc(:,i,j),PressIntPycFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(UIntPycP(:,i,j),UIntPyc(:,i,j),UIntPycFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntPycP(:,i,j),VIntPyc(:,i,j),VIntPycFour(:,i,j),w0*dt1/5.,bw/5.,Num)
 
   !  First Harmonic
  call BandPassFilter(PressIntPycH(:,i,j),PressIntPyc(:,i,j),PressIntPycFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(UIntPycH(:,i,j),UIntPyc(:,i,j),UIntPycFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VIntPycH(:,i,j),VIntPyc(:,i,j),VIntPycFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)
 

  end do
 end do

