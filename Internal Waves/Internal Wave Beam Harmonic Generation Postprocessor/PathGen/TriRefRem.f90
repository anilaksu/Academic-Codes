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

do i=1,il

  write(400,*) xtrir(40,i),ztrir(40,i)
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

dx=7.*abs(xtrir(1,1)-xtrir(1,2))

dz=7.*abs(ztrir(1,1)-ztrir(1,2))

print*,"dx and dz",dx,dz

do i=1,Ntri


!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointTriR(i,:,:),xtrir(i,:),ztrir(i,:),il,nzr,nxr,Nptrir(i),dx,dz)

print*,Nptrir(i)


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

do i=1,Nri

  do j=1,Num
 

      call interpVal(IRTriR(i,:,:),PressTriR(j,i,:),Press(j,:),PointTriR(i,:,1),PointTriR(i,:,2),nxr,nzr,Nptrir(i),il)
   
      call interpVal(IRTriR(i,:,:),UTriR(j,i,:),u(j,:),PointTriR(i,:,1),PointTriR(i,:,2),nxr,nzr,Nptrir(i),il)

      call interpVal(IRTriR(i,:,:),VTriR(j,i,:),v(j,:),PointTriR(i,:,1),PointTriR(i,:,2),nxr,nzr,Nptrir(i),il)

  
  end do

end do



! the arclength starting from the entry of triangular region 

arctrir=0.


do i=2,Ntri
 
 arctrir(i)=arctrir(i-1)+dsqrt((xtrir(i,1)-xtrir(i-1,1))**2.+(ztrir(i,1)-ztrir(i-1,1))**2.)

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
   call BandPassFilter(PressTriRP(:,i,j),PressTriR(:,i,j),PressTriRFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(UTriRP(:,i,j),UTriR(:,i,j),UTriRFour(:,i,j),w0*dt1/5.,bw/5.,Num)

   call BandPassFilter(VTriRP(:,i,j),VTriR(:,i,j),VTriRFour(:,i,j),w0*dt1/5.,bw/5.,Num)
 
   !  First Harmonic
  call BandPassFilter(PressTriRH1(:,i,j),PressTriR(:,i,j),PressTriRFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(UTriRH1(:,i,j),UTriR(:,i,j),UTriRFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(VTriRH1(:,i,j),VTriR(:,i,j),VTriRFour(:,i,j),2.*w0*dt1/5.,bw/5.,Num)
 
   !  Second Harmonic
  call BandPassFilter(PressTriRH2(:,i,j),PressTriR(:,i,j),PressTriRFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(UTriRH2(:,i,j),UTriR(:,i,j),UTriRFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)

  call BandPassFilter(VTriRH2(:,i,j),VTriR(:,i,j),VTriRFour(:,i,j),3.*w0*dt1/5.,bw/5.,Num)
 
  end do
 end do

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

open(322,file='TotPU1Trir400.dat',status='unknown')

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

  write(600,*) arcref(i)/lambdax,EnergyRef(i)/Enon
 
end do

open(601,file='PEnergyTriR200.dat',status='unknown')

do i=1,Ntri

  write(601,*) arcref(i)/lambdax,EnergyRefP(i)/Enon
 
end do


open(602,file='1HEnergyTriR200.dat',status='unknown')

do i=1,Ntri
  write(602,*) arcref(i)/lambdax,EnergyRefH1(i)/Enon
 
end do


open(603,file='2HEnergyTriR200.dat',status='unknown')

do i=1,Ntri

  write(603,*) arcref(i)/lambdax,EnergyRefH2(i)/Enon
 
end do




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

  write(610,*) arcref(i)/lambdax,EnergyRef(i)/Enon
 
end do

open(611,file='PEnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(611,*) arcref(i)/lambdax,EnergyRefP(i)/Enon
 
end do


open(612,file='1HEnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(612,*) arcref(i)/lambdax,EnergyRefH1(i)/Enon
 
end do


open(613,file='2HEnergyTriR300.dat',status='unknown')

do i=1,Ntri

  write(613,*) arcref(i)/lambdax,EnergyRefH2(i)/Enon
 
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

  write(620,*) arcref(i)/lambdax,EnergyRef(i)/Enon
 
end do

open(621,file='PEnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(621,*) arcref(i)/lambdax,EnergyRefP(i)/Enon
 
end do


open(622,file='1HEnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(622,*) arcref(i)/lambdax,EnergyRefH1(i)/Enon
 
end do


open(623,file='2HEnergyTriR400.dat',status='unknown')

do i=1,Ntri

  write(623,*) arcref(i)/lambdax,EnergyRefH2(i)/Enon
 
end do



