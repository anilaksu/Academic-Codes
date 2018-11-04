allocate(IR(N(1),il,nxr*nzr))


! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once



!to determine the distance in which the interpolation data to be extracted 

dz=abs(Qz(2)-Qz(1))

dx=abs(Qx(2)-Qx(1))

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,Point(1,:,:),x(1,:),z(1,:),il,nzr,nxr,Np(1),dx,dz)

print*,Np(1)


!allocate them but it will change in every iteration 
allocate(xi(Np(1)))

allocate(zi(Np(1)))


do j=1,Np(1)

xi(j)=x1(Point(1,j,1),Point(1,j,2))

zi(j)=z1(Point(1,j,1),Point(1,j,2))


end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IR(1,:,:),x(1,:),xi,z(1,:),zi,Np(1),il)


deallocate(xi)

deallocate(zi)

do j=1,Num
 

    call interpVal(IR(1,:,:),PressInt(j,1,:),Press(j,:),Point(1,:,:1),Point(1,:,2),nxr,nzr,Np(1),il)
 
    call interpVal(IR(1,:,:),UInt(j,1,:),u(j,:),Point(1,:,:1),Point(1,:,2),nxr,nzr,Np(1),il)

    call interpVal(IR(1,:,:),VInt(j,1,:),v(j,:),Point(1,:,:1),Point(1,:,2),nxr,nzr,Np(1),il)

  
 end do

!print*,UInt(1,1,:)

print*,"okey"


! bandpass filtering in incident wave beam region 

 do i=1,1
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


open(405,file='TEnergyFluxIWB.dat',status='unknown')

do j=1,il

E1=PressInt(300,1,j)*(dcos(pi/4.)*UInt(300,1,j)+dcos(pi/4.)*VInt(300,1,j))/Enon

write(405,*) (zdim+zdim*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1

end do 

open(406,file='PEnergyFluxIWB.dat',status='unknown')

do j=1,il

E1=PressIntP(300,1,j)*(dcos(pi/4.)*UIntP(300,1,j)+dcos(pi/4.)*VIntP(300,1,j))/Enon

write(406,*) (zdim+zdim*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1

end do 

open(407,file='HEnergyFluxIWB.dat',status='unknown')

do j=1,il

E1=PressIntH(300,1,j)*(dcos(pi/4.)*UIntH(300,1,j)+dcos(pi/4.)*VIntH(300,1,j))/Enon

write(407,*) (zdim+zdim*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1

end do 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Pycnocline Calculations               !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!this do loop is used to calculate the number of points in pycnocline

Npyc=3

! to calculate the quadrature points in pycnocline 

! we start from x(N(5),5) to reduce the computational effort 

! And also our beam enters to pycnocline at that point 
 call QuadPoints(Qx,Qz,1.,0.,x(N(20),20),z(N(20),20)+h/2,h/2,5) 

!print*,Qx
!print*,Qz
! now lets find the interpolation points in pycnocline


allocate(xpyc(Npyc,20))
allocate(zpyc(Npyc,20))


! the interpolated data in pycnocline
allocate(PressIntPyc(Num,Npyc,20))

allocate(UIntPyc(Num,Npyc,20))

allocate(VIntPyc(Num,Npyc,20))

! fourier transformed data
allocate(PressIntPycFour(Num,Npyc,220))

allocate(UIntPycFour(Num,Npyc,20))

allocate(VIntPycFour(Num,Npyc,20))

! primary frequency data 
allocate(PressIntPycP(Num,Npyc,20))

allocate(UIntPycP(Num,Npyc,20))

allocate(VIntPycP(Num,Npyc,20))

! first harmonic data
allocate(PressIntPycH(Num,Npyc,20))

allocate(UIntPycH(Num,Npyc,20))

allocate(VIntPycH(Num,Npyc,20))

! group velocities to denote the unit normal in integration
allocate(CgxPyc(Npyc,20))

allocate(CgzPyc(Npyc,20))

! energy in pycnocline
allocate(EnergyPyc(Num,Npyc))

allocate(EnergyPycP(Num,Npyc))

allocate(EnergyPycH(Num,Npyc))



 CgxPyc=1.
 CgzPyc=0.

ds=2.*lambdax
do i=1,Npyc 

do j=1,20
! ds increment in x 
xpyc(i,j)=(i)*ds+Qx(j)

!z is same along the path 

zpyc(i,j)=Qz(j)

end do


end do




allocate(Np1(Npyc))

allocate(IRPyc(Npyc,20,nxr*nzr))

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
 
  do i=1,Npyc


   call interpVal(IRPyc(i,:,:),PressIntPyc(j,i,:),Press(j,:),PointPyc(i,:,:1),PointPyc(i,:,2),nxr,nzr,Np1(i),5)
 
   call interpVal(IRPyc(i,:,:),UIntPyc(j,i,:),u(j,:),PointPyc(i,:,:1),PointPyc(i,:,2),nxr,nzr,Np1(i),5)

   call interpVal(IRPyc(i,:,:),VIntPyc(j,i,:),v(j,:),PointPyc(i,:,:1),PointPyc(i,:,2),nxr,nzr,Np1(i),5)

  end do

 end do
! bandpass filtering in pycnocline region 
do i=1,Npyc
   do j=1,20

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

!first x-section
open(408,file='TEnergyFluxPyc1.dat',status='unknown')

do j=1,20

    write(408,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPyc(300,1,j)*UIntPyc(300,1,j)/Enon

end do 

open(409,file='PEnergyFluxPyc1.dat',status='unknown')

do j=1,20

    write(409,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPycP(300,1,j)*UIntPycP(300,1,j)/Enon

end do 

open(410,file='HEnergyFluxPyc1.dat',status='unknown')

do j=1,20

    write(410,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPycH(300,1,j)*UIntPycH(300,1,j)/Enon

end do 

!second x-section
open(411,file='TEnergyFluxPyc2.dat',status='unknown')

do j=1,20

  write(411,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPyc(300,2,j)*UIntPyc(300,2,j)/Enon

end do 

open(412,file='PEnergyFluxPyc2.dat',status='unknown')

do j=1,20

E1=PressIntPycP(300,2,j)*UIntPycP(300,2,j)/Enon


 write(412,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1
end do 

open(413,file='HEnergyFluxPyc2.dat',status='unknown')

do j=1,20

  write(413,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPycH(300,2,j)*UIntPycH(300,2,j)/Enon

end do 

!second x-section
open(414,file='TEnergyFluxPyc3.dat',status='unknown')

do j=1,20

    write(414,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPyc(300,3,j)*UIntPyc(300,3,j)/Enon

end do 

open(415,file='PEnergyFluxPyc3.dat',status='unknown')

do j=1,20

    write(415,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPycP(300,3,j)*UIntPycP(300,3,j)/Enon

end do 

open(416,file='HEnergyFluxPyc3.dat',status='unknown')

do j=1,20

    write(416,*) (0.5*h+0.5*h*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,PressIntPycH(300,3,j)*UIntPycH(300,3,j)/Enon

end do 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Reflected Beam Calculations           !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(xref(N(1),20))
allocate(zref(N(1),20))


allocate(CgxRef(N(1),20))
allocate(CgzRef(N(1),20))

allocate(IRRef(N(1),20,nxr*nzr))

!the interpolated data
allocate(PressIntRefFour(Num,N(1),20))

allocate(UIntRefFour(Num,N(1),20))

allocate(VIntRefFour(Num,N(1),20))

!the interpolated data
allocate(PressIntRef(Num,N(1),20))

allocate(UIntRef(Num,N(1),20))

allocate(VIntRef(Num,N(1),20))

! the primary frequency data
allocate(PressIntRefP(Num,N(1),20))

allocate(UIntRefP(Num,N(1),20))

allocate(VIntRefP(Num,N(1),20))

! first harmonic data
allocate(PressIntRefH(Num,N(1),20))

allocate(UIntRefH(Num,N(1),20))

allocate(VIntRefH(Num,N(1),20))

allocate(PointRef(N(1),nxr*nzr,2))


! energy in pycnocline
allocate(EnergyRef(Num,N(1)))

allocate(EnergyRefP(Num,N(1)))

allocate(EnergyRefH(Num,N(1)))



! the dimension of the point array is N(1) as it is the largest 
allocate(Np2(N(1)))

! the width of the reflected beam is determined as the difference between x-coordinate of the IWB in to pycnocline


! zdimr =( x(N(1),1)-x(N(5),5))/2.

 xcenr =( x(N(1),1)+x(N(20),20))/2.

 zcenr = z(N(20),20)

 ! calculation of the reflected beam path 

 ! it is not generated using path command as it does not differentiate positive and negative direction  


 do i=1 ,N(1)-1

 do j=1 ,20


 
 xref(i,j)=2.*xcenr-x(i,j)

 zref(i,j)=z(i,j)

 ! the group velocities are calculated here it shares the same magnitude with IWB but it propogates in negative z direction

 CgzRef(i,j)=-1.*Cgz(i,j)
 
 CgxRef(i,j)=Cgx(i,j)
 
 !print*,xref(i,j),zref(i,j)

 end do 

 end do 


! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once

!to determine the distance in which the interpolation data to be extracted 

dz=2.*abs(z(2,1)-z(1,1))

dx=2.*abs(x(2,1)-x(1,1))

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointRef(1,:,:),xref(1,:),zref(1,:),20,nzr,nxr,Np2(1),dx,dz)

print*,"Energy Ref",Np2(1),dx,dz

!allocate them but it will change in every iteration 
allocate(xi(Np2(1)))

allocate(zi(Np2(1)))


do j=1,Np2(1)

 xi(j)=x1(PointRef(1,j,1),PointRef(1,j,2))

 zi(j)=z1(PointRef(1,j,1),PointRef(1,j,2))

print*,j,xi(j),zi(j),PointRef(1,j,1),PointRef(1,j,2)
print*,x1(PointRef(1,j,1),PointRef(1,j,2)),z1(PointRef(1,j,1),PointRef(1,j,2))
end do 

 call interpolation2D(IRRef(1,:,:),xref(1,:),xi,zref(1,:),zi,Np2(1),20)

deallocate(xi)

deallocate(zi)


print*,"okeyref"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!    Reflected Pressure ,Velocity,Energy         !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 


do j=1,Num
 
print*,"okey ref ",j

    call interpVal(IRRef(1,:,:),PressIntRef(j,1,:),Press(j,:),PointRef(1,:,1),PointRef(1,:,2),nxr,nzr,Np2(1),20)
 if (j==1 ) then 
print*,PressIntRef(j,1,:)
end if 
 
print*,"okey ref ",j

    call interpVal(IRRef(1,:,:),UIntRef(j,1,:),u(j,:),PointRef(1,:,1),PointRef(1,:,2),nxr,nzr,Np2(1),20)

print*,"okey ref ",j

    call interpVal(IRRef(1,:,:),VIntRef(j,1,:),v(j,:),PointRef(1,:,1),PointRef(1,:,2),nxr,nzr,Np2(1),20)

print*,"okey ref ",j

 end do

print*,"okey ref "

! bandpass filtering in reflection region 
do i=1,1
   do j=1,20

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

open(417,file='TEnergyFluxRef.dat',status='unknown')

do j=1,20

E1=PressIntRef(300,1,j)*(dcos(pi/4.)*UIntRef(300,1,j)-dsin(pi/4.)*VIntRef(300,1,j))/Enon

write(417,*) (zdim+zdim*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1

end do 

open(418,file='PEnergyFluxRef.dat',status='unknown')

do j=1,20

E1=PressIntRefP(300,1,j)*(dcos(pi/4.)*UIntRefP(300,1,j)-dsin(pi/4.)*VIntRefP(300,1,j))/Enon

write(418,*) (zdim+zdim*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1

end do 

open(419,file='HEnergyFluxRef.dat',status='unknown')


do j=1,20

E1=PressIntRefH(300,1,j)*(dcos(pi/4.)*UIntRefH(300,1,j)-dsin(pi/4.)*VIntRefH(300,1,j))/Enon

write(419,*) (zdim+zdim*dcos((2.*j-1.)*pi/(2.*20)))/lambdax,E1

end do 

!print*,UIntRef(1,1,:)
!open(321,file='RefPoint.dat',status='unknown')

!do i=1,20 

 !write(321,*) xref(1,i),zref(1,i)

