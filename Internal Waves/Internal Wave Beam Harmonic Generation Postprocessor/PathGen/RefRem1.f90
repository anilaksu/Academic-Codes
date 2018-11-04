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

dx=0.5*abs(xref(1,1)-xref(1,2))

dz=0.5*abs(zref(1,1)-zref(1,2))

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
   call BandPassFilter(PressRefP(:,i,j),PressRef(:,i,j),PressRefFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(URefP(:,i,j),URef(:,i,j),URefFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VRefP(:,i,j),VRef(:,i,j),VRefFour(:,i,j),w0*dt1/5.,bw/5.,Num)
 
   !  First Harmonic
  call BandPassFilter(PressRefH1(:,i,j),PressRef(:,i,j),PressRefFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(URefH1(:,i,j),URef(:,i,j),URefFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(VRefH1(:,i,j),VRef(:,i,j),VRefFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressRefH2(:,i,j),PressRef(:,i,j),PressRefFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(URefH2(:,i,j),URef(:,i,j),URefFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VRefH2(:,i,j),VRef(:,i,j),VRefFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)
 
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

open(427,file='PriPU1Ref300out.dat',status='unknown')

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





