!Finally Energy Integration of wave 

do j=1,Num
 do i=1,N(1)-1

   call IntegLine(Energy(i,j),PressInt(j,i,:),UInt(j,i,:),VInt(j,i,:),Cgx(i,:),Cgz(i,:),zdim,5)

 end do 
end do



!Finally Energy Integration of primary wave 

do j=1,Num
 do i=1,N(1)-1

   call IntegLine(EnergyP(i,j),PressIntP(j,i,:),UIntP(j,i,:),VIntP(j,i,:),Cgx(i,:),Cgz(i,:),zdim,5)

 end do 
end do


!Finally Energy Integration of harmonic wave 

do j=1,Num
 do i=1,N(1)-1

   call IntegLine(EnergyH(i,j),PressIntH(j,i,:),UIntH(j,i,:),VIntH(j,i,:),Cgx(i,:),Cgz(i,:),zdim,5)

 end do 
end do


! L2 norm of the pressure velocity product 


open(310,file='L2PVt2.dat',status='unknown')

do i=1,N(5)
  do j=1,5
    write(310,*) x(i,j),z(i,j),abs(PressInt(2,i,j))*dsqrt(UInt(2,i,j)**2.+VInt(2,i,j)**2.)
  
  end do   
   
end do 

open(311,file='L2PVt100.dat',status='unknown')

do i=1,N(5)
  do j=1,5
    write(311,*) x(i,j),z(i,j),abs(PressInt(100,i,j))*dsqrt(UInt(100,i,j)**2.+VInt(100,i,j)**2.)
  
  end do   
   
end do 


open(312,file='L2PVt300.dat',status='unknown')

do i=1,N(5)
  do j=1,5
    write(31,*) x(i,j),z(i,j),abs(PressInt(300,i,j))*dsqrt(UInt(300,i,j)**2.+VInt(300,i,j)**2.)
  
  end do   
   
end do 



! total energy of wave along the path at sample time steps
open(316,file='EnAlPathTt2.dat',status='unknown')

do i=1,N(5)

    write(316,*) path1(i)/lambdax,Energy(i,2)/Enon
     
   
end do 


open(317,file='EnAlPathTt100.dat',status='unknown')

do i=1,N(5)

    write(317,*) path1(i)/lambdax,Energy(i,100)/Enon

end do 

open(318,file='EnAlPathTt300.dat',status='unknown')

do i=1,N(5)

    write(318,*) path1(i)/lambdax,Energy(i,300)/Enon

end do 


! energy of the primary wave along the path at sample time steps
open(310,file='EnAlPatht2.dat',status='unknown')

do i=1,N(5)

    write(310,*) path1(i)/lambdax,EnergyP(i,2)/Enon
     
   
end do 

open(311,file='EnAlPatht100.dat',status='unknown')

do i=1,N(5)

    write(311,*) path1(i)/lambdax,EnergyP(i,100)/Enon

end do 

open(312,file='EnAlPatht300.dat',status='unknown')

do i=1,N(5)

    write(312,*) path1(i)/lambdax,EnergyP(i,300)/Enon

end do 

! energy of the harmonic wave along the path at sample time steps
open(313,file='EnAlPathHt2.dat',status='unknown')

do i=1,N(5)

    write(313,*) path1(i)/lambdax,EnergyH(i,2)/Enon

end do 

open(314,file='EnAlPathHt100.dat',status='unknown')

do i=1,N(5)

    write(314,*) path1(i)/lambdax,EnergyH(i,100)/Enon

end do 

open(315,file='EnAlPathHt300.dat',status='unknown')

do i=1,N(5)

   write(315,*) path1(i)/lambdax,EnergyH(i,300)/Enon

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

dz=50.*abs(Qz(2)-Qz(1))

dx=dz

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointPyc(i,:,:),xpyc(i,:),zpyc(i,:),5,nzr,nxr,Np1(i),dx,dz)
!print*,"Okey",Np1(i),dx,dz,h



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


! output the time series 
open(200,file='VelocityIso1.dat',status='unknown')
   


do j=1,Num

 write(200,*) (j-1)*ds ,UIntPyc(j,1,1)

end do

open(201,file='VelocityFour.dat',status='unknown')
   


do j=1,Num

 write(201,*) real(UIntPycFour(j,1,1))

end do

open(202,file='VelocityFilter.dat',status='unknown')
   

do j=1,Num

 write(202,*) UIntPycP(j,1,1)

end do

!Finally Energy Integration of wave 

do j=1,Num
 do i=1,Npyc

   call IntegLine(EnergyPyc(i,j),PressIntPyc(j,i,:),UIntPyc(j,i,:),VIntPyc(j,i,:),CgxPyc(i,:),CgzPyc(i,:),h/2.,5)

 end do 
end do

!Finally Energy Integration of primary frequency  

do j=1,Num
 do i=1,Npyc

   call IntegLine(EnergyPycP(i,j),PressIntPycP(j,i,:),UIntPycP(j,i,:),VIntPycP(j,i,:),CgxPyc(i,:),CgzPyc(i,:),h/2.,5)

 end do 
end do

!Finally Energy Integration of first harmonic

do j=1,Num
 do i=1,Npyc

   call IntegLine(EnergyPycH(i,j),PressIntPycH(j,i,:),UIntPycH(j,i,:),VIntPycH(j,i,:),CgxPyc(i,:),CgzPyc(i,:),h/2.,5)

 end do 
end do
! L2 norm of the pressure velocity product 


open(397,file='PL2PVt2.dat',status='unknown')

do i=1,Npyc
  do j=1,5
    write(397,*) xpyc(i,j),zpyc(i,j),abs(PressIntPyc(2,i,j))*dsqrt(UIntPyc(2,i,j)**2.+VIntPyc(2,i,j)**2.)
  
  end do   
   
end do 

open(398,file='PL2PVt100.dat',status='unknown')

do i=1,Npyc
  do j=1,5
    write(398,*) xpyc(i,j),zpyc(i,j),abs(PressIntPyc(100,i,j))*dsqrt(UIntPyc(100,i,j)**2.+VIntPyc(100,i,j)**2.)
  
  end do   
   
end do 


open(399,file='PL2PVt300.dat',status='unknown')

do i=1,Npyc
  do j=1,5
    write(399,*) xpyc(i,j),zpyc(i,j),abs(PressIntPyc(300,i,j))*dsqrt(UIntPyc(300,i,j)**2.+VIntPyc(300,i,j)**2.)
  
  end do   
   
end do 


! total energy of wave along the path at sample time steps
open(400,file='PycEnAlPathTt2.dat',status='unknown')

do i=1,Npyc

    write(400,*) (i-1)*ds/lambdax,EnergyPyc(i,2)/Enon
     
   
end do 

open(401,file='PycEnAlPathTt100.dat',status='unknown')

do i=1,Npyc

    write(401,*) (i-1)*ds/lambdax,EnergyPyc(i,100)/Enon

end do 

open(402,file='PycEnAlPathTt300.dat',status='unknown')

do i=1,Npyc

    write(402,*) (i-1)*ds/lambdax,EnergyPyc(i,300)/Enon

end do 

! energy of primary frequency along the path at sample time steps

open(403,file='PycEnAlPatht2.dat',status='unknown')

do i=1,Npyc

    write(403,*) (i-1)*ds/lambdax,EnergyPycP(i,2)/Enon
     
   
end do 

open(404,file='PycEnAlPatht100.dat',status='unknown')

do i=1,Npyc

    write(404,*) (i-1)*ds/lambdax,EnergyPycP(i,100)/Enon

end do 

open(405,file='PycEnAlPatht300.dat',status='unknown')

do i=1,Npyc

    write(405,*) (i-1)*ds/lambdax,EnergyPycP(i,300)/Enon

end do 


! energy of first harmonic along the path at sample time steps

open(403,file='PycEnAlPathH2.dat',status='unknown')

do i=1,Npyc

    write(403,*) (i-1)*ds/lambdax,EnergyPycH(i,2)/Enon
     
   
end do 

open(404,file='PycEnAlPathH100.dat',status='unknown')

do i=1,Npyc

    write(404,*) (i-1)*ds/lambdax,EnergyPycH(i,100)/Enon

end do 

open(405,file='PycEnAlPathH300.dat',status='unknown')

do i=1,Npyc

    write(405,*) (i-1)*ds/lambdax,EnergyPycH(i,300)/Enon

end do 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!	   Reflected Beam Calculations           !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(xref(N(1),5))
allocate(zref(N(1),5))


allocate(CgxRef(N(1),5))
allocate(CgzRef(N(1),5))

allocate(IRRef(N(1),5,nxr*nzr))

!the interpolated data
allocate(PressIntRefFour(Num,N(1),5))

allocate(UIntRefFour(Num,N(1),5))

allocate(VIntRefFour(Num,N(1),5))

!the interpolated data
allocate(PressIntRef(Num,N(1),5))

allocate(UIntRef(Num,N(1),5))

allocate(VIntRef(Num,N(1),5))

! the primary frequency data
allocate(PressIntRefP(Num,N(1),5))

allocate(UIntRefP(Num,N(1),5))

allocate(VIntRefP(Num,N(1),5))

! first harmonic data
allocate(PressIntRefH(Num,N(1),5))

allocate(UIntRefH(Num,N(1),5))

allocate(VIntRefH(Num,N(1),5))

allocate(PointRef(N(1),nxr*nzr,2))


! energy in pycnocline
allocate(EnergyRef(Num,N(1)))

allocate(EnergyRefP(Num,N(1)))

allocate(EnergyRefH(Num,N(1)))



! the dimension of the point array is N(1) as it is the largest 
allocate(Np2(N(1)))

! the width of the reflected beam is determined as the difference between x-coordinate of the IWB in to pycnocline


! zdimr =( x(N(1),1)-x(N(5),5))/2.

 xcenr =( x(N(1),1)+x(N(5),5))/2.

 zcenr = z(N(5),5)

 ! calculation of the reflected beam path 

 ! it is not generated using path command as it does not differentiate positive and negative direction  


 do i=1 ,N(1)-1

 do j=1 ,5


 
 xref(i,j)=2.*xcenr-x(i,j)

 zref(i,j)=z(i,j)

 ! the group velocities are calculated here it shares the same magnitude with IWB but it propogates in negative z direction

 CgzRef(i,j)=-1.*Cgz(i,j)
 
 CgxRef(i,j)=Cgz(i,j)
 
 

 end do 

 end do 


! in that do loop I generate the interpolation matrix and store them for each time step

! as the same interpolation matrix is used at each time step I am generating it once


do i=1,N(1)-1

!to determine the distance in which the interpolation data to be extracted 

dz=abs(z(2,1)-z(1,1))/1.5

dx=abs(x(2,1)-x(1,1))/1.5

!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointRef(i,:,:),xref(i,:),zref(i,:),5,nzr,nxr,Np2(i),dx,dz)

!print*,"Energy Ref",Np2(i),dz,dx

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
 
  do i=1,N(1)-1

    call interpVal(IRRef(i,:,:),PressIntRef(j,i,:),Press(j,:),PointRef(i,:,:1),PointRef(i,:,2),nxr,nzr,Np2(i),5)
 
    call interpVal(IRRef(i,:,:),UIntRef(j,i,:),u(j,:),PointRef(i,:,:1),PointRef(i,:,2),nxr,nzr,Np2(i),5)

    call interpVal(IRRef(i,:,:),VIntRef(j,i,:),v(j,:),PointRef(i,:,:1),PointRef(i,:,2),nxr,nzr,Np2(i),5)

  end do

 end do


! bandpass filtering in reflection region 
do i=1,N(1)-1
   do j=1,5

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

!Finally Energy Integration of total wave in reflection

do j=1,Num
 do i=1,N(1)-1

   call IntegLine(EnergyRef(i,j),PressIntRef(j,i,:),UIntRef(j,i,:),VIntRef(j,i,:),CgxRef(i,:),CgzRef(i,:),zdim,5)

 end do 
end do

!Finally Energy Integration of primary wave in reflection
do j=1,Num
 do i=1,N(1)-1


   call IntegLine(EnergyRefP(i,j),PressIntRefP(j,i,:),UIntRefP(j,i,:),VIntRefP(j,i,:),CgxRef(i,:),CgzRef(i,:),zdim,5)

 end do 
end do

!Finally Energy Integration of primary wave in reflection
do j=1,Num
 do i=1,N(1)-1

   call IntegLine(EnergyRefH(i,j),PressIntRefH(j,i,:),UIntRefH(j,i,:),VIntRefP(j,i,:),CgxRef(i,:),CgzRef(i,:),zdim,5)

 end do 
end do
! L2 norm of the pressure velocity product 


open(497,file='RL2PVt2.dat',status='unknown')

do i=1,N(1)-1
  do j=1,5
    write(497,*) xref(i,j),zref(i,j),PressIntRef(2,i,j)*dsqrt(UIntRef(2,i,j)**2.+VIntRef(2,i,j)**2.)
  
  end do   
   
end do 

open(498,file='RL2PVt100.dat',status='unknown')

do i=1,N(1)-1
  do j=1,5
    write(498,*) xref(i,j),zref(i,j),PressIntRef(100,i,j)*dsqrt(UIntRef(100,i,j)**2.+VIntRef(100,i,j)**2.)
  
  end do   
   
end do 


open(499,file='RL2PVt300.dat',status='unknown')

do i=1,N(1)-1
  do j=1,5
    write(499,*) xref(i,j),zref(i,j),PressIntRef(300,i,j)*dsqrt(UIntRef(300,i,j)**2.+VIntRef(300,i,j)**2.)
  
  end do   
   
end do 

! total energy of wave along the path at sample time steps
open(503,file='RefEnAlPathTt2.dat',status='unknown')

do i=1,N(1)-1

    write(503,*) path1(i)/lambdax,EnergyRef(i,2)/Enon
     
   
end do 

open(504,file='RefEnAlPathTt100.dat',status='unknown')

do i=1,N(1)-1

    write(504,*) path1(i)/lambdax,EnergyRef(i,100)/Enon

end do 

open(505,file='RefEnAlPathTt200.dat',status='unknown')

do i=1,N(1)-1

    write(505,*) path1(i)/lambdax,EnergyRef(i,300)/Enon

end do 


! energy of primary frequency along the path at sample time steps
open(506,file='RefEnAlPatht2.dat',status='unknown')

do i=1,N(1)-1

    write(506,*) path1(i)/lambdax,EnergyRefP(i,2)/Enon
     
   
end do 

open(507,file='RefEnAlPatht100.dat',status='unknown')

do i=1,N(1)-1

    write(507,*) path1(i)/lambdax,EnergyRefP(i,100)/Enon

end do 

open(508,file='RefEnAlPatht300.dat',status='unknown')

do i=1,N(1)-1

    write(508,*) path1(i)/lambdax,EnergyRefP(i,300)/Enon

end do 

! energy of first harmonic along the path at sample time steps
open(509,file='RefEnAlPathHt2.dat',status='unknown')

do i=1,N(1)-1

    write(509,*) path1(i)/lambdax,EnergyRefH(i,2)/Enon
     
   
end do 

open(510,file='RefEnAlPathHt100.dat',status='unknown')

do i=1,N(1)-1

    write(510,*) path1(i)/lambdax,EnergyRefH(i,100)/Enon

end do 

open(511,file='RefEnAlPathHt300.dat',status='unknown')

do i=1,N(1)-1

    write(511,*) path1(i)/lambdax,EnergyRefH(i,300)/Enon

end do 


