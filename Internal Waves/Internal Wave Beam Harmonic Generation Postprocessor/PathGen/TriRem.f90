

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

print*," okey 3"
ds=0.0001

do i=1,Ntri

 zdimr(i)=0.

 !x(N(il)+i,1)

 zdum=z(N(il)+i,1)


 ! in this loop , we find the length of the integration path 
  do j=1,2000 
 
   zdum=zdum+ds

   call BVprofile(BV,N1,r,zdum,z0,dt)

    if(zdum>z(N(il),il)+h .or. BV<N1) then 

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


   !if (zdum < ztri(i,il)) then 

    ! zdimr(i)=zdimr(i)-eps 

     !goto 600     

   !end if 

 end do

print*,"starting point",x(N(il),1),z(N(il),1)

open(300,file='TriPath.dat',status='unknown')

do i=1,il

  write(300,*) xtri(1,i),ztri(1,i)
  !print*,N(i)
end do

open(301,file='FPointLineP.dat',status='unknown')

do i=1,Ntri

  write(301,*) xtri(i,1),ztri(i,1)
  !print*,N(i)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
! Interpolation in Triangular region for primary frequency !
!                                                          !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(Nptri(Ntri))

allocate(IRTri(Ntri,il,nxr*nzr))

allocate(PointTri(Ntri,nxr*nzr,2))

allocate(PressIntTri(Num,Ntri,il))

allocate(UIntTri(Num,Ntri,il))

allocate(VIntTri(Num,Ntri,il))

allocate(EnergyTri(Ntri))

!allocate(EnergyTotTri(Ntri))

dx=abs(xtri(1,1)-xtri(1,2))

dz=abs(ztri(1,1)-ztri(1,2))

print*,"dx and dz",dx,dz
do i=1,Ntri


!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointTri(i,:,:),xtri(i,:),ztri(i,:),il,nzr,nxr,Nptri(i),dx,dz)

print*,Nptri(i)


!allocate them but it will change in every iteration 
allocate(xi(Nptri(i)))

allocate(zi(Nptri(i)))


   do j=1,Nptri(i)

     xi(j)=x1(PointTri(i,j,1),PointTri(i,j,2))

     zi(j)=z1(PointTri(i,j,1),PointTri(i,j,2))
 ! print*,xi(j),zi(j)

   end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IRTri(i,:,:),xtri(i,:),xi,ztri(i,:),zi,Nptri(i),il)


  deallocate(xi)

  deallocate(zi)

end do

!print*,"okey tri"
!print*,IRTri(1,:,1)

do i=1,Ntri

  do j=1,1
 

      call interpVal(IRTri(i,:,:),PressIntTri(j,i,:),Press(j,:),PointTri(i,:,1),PointTri(i,:,2),nxr,nzr,Nptri(i),il)
   
      call interpVal(IRTri(i,:,:),UIntTri(j,i,:),u(j,:),PointTri(i,:,1),PointTri(i,:,2),nxr,nzr,Nptri(i),il)

      call interpVal(IRTri(i,:,:),VIntTri(j,i,:),v(j,:),PointTri(i,:,1),PointTri(i,:,2),nxr,nzr,Nptri(i),il)

  
  end do

end do


open(302,file='TriPress.dat',status='unknown')

do i=1,il

  write(302,*) (i-1)*zdimr(1)/il,PressIntTri(1,5,i)
  !print*,N(i)
end do

open(303,file='TriU.dat',status='unknown')

do i=1,il

  write(303,*) (i-1)*zdimr(1)/il,UIntTri(1,5,i)
  !print*,N(i)
end do

open(304,file='TriV.dat',status='unknown')

do i=1,il

  write(304,*) (i-1)*zdimr(1)/il,VIntTri(1,5,i)
  !print*,N(i)
end do

open(305,file='GroupVel.dat',status='unknown')

do i=1,il

  write(305,*) Cgxtri(5,i),Cgztri(5,i)
  !print*,N(i)
end do

! the arclength starting from the entry of triangular region 

arctri=0.

do i=2,Ntri
 
 arctri(i)=arctri(i-1)+dsqrt((xtri(i,1)-xtri(i-1,1))**2.+(ztri(i,1)-ztri(i-1,1))**2.)

end do 

do i=1,Ntri

 call PressVelocityInteg(EnergyTri(i),PressIntTri(1,i,:),UIntTri(1,i,:),VIntTri(1,i,:),il,Cgxtri(i,:),Cgztri(i,:),zdimr(i))


end do 

open(306,file='TriEnergy.dat',status='unknown')

do i=1,Ntri

  write(306,*) arctri(i),EnergyTri(i)
  !print*,N(i)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!  Triangular Reflection Region First Harmonics  !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Ntrih is the number of points of first harmonic in triangular region 

ds=0.005


do i=1,1000

 if (ds*i > x(N(1),1)-x(N(il),il)) then

  Ntrih=i-1

  goto 800
  end if


end do

800 continue

print*,"Ntrih",Ntrih


allocate(xtrih(Ntrih,il))

allocate(ztrih(Ntrih,il))


allocate(zdimrh(Ntrih))

allocate(Cgxtrih(Ntrih,il))

allocate(Cgztrih(Ntrih,il))

allocate(arctrih(Ntrih))

ds=0.0001


! here we are  changing  kx keep that in mind
kxh1=2.*kx

wh1=2.*w0

do i=1,Ntrih

 zdimrh(i)=0.

 !x(N(il)+i,1)

 zdum=z(N(il),il)


 ! in this loop , we find the length of the integration path 
  do j=1,2000 
 
   zdum=zdum+ds

   if (i==1 .and. j==1) then 

     print*,zdum

   end if 
   call BVprofile(BV,N1,r,zdum,z0,dt)

    if(zdum>z(N(il),il)+h .or. BV<2.*N1) then 

    goto 900

    end if 

    
     
     
      kzh1=dsqrt((BV**2.)*(kxh1**2.)/(wh1**2.)-kxh1**2.)

   ! print*,"kzh1",kzh1,N1,BV,wh1
     zdimrh(i)=zdimrh(i)+ds*dsqrt(1.+(kxh1/kzh1)**2.)
   !if(i==1) then
   !print*,zdimr(i),j
   !end if 


   end do

   900 continue 
    
    print*,zdimrh(i),z(N(il),il)+h,zdum
   
    xtrih(i,1)=x(N(il),il)+i*ds

    ztrih(i,1)=z(N(il),il)
  
    call GroupVel(Cgxtrih(i,1),Cgztrih(i,1),N1,ztrih(i,1),z0,dt,r,kxh1,wh1)

   do j=1,il-1

      call BVprofile(BV,N1,r,ztrih(i,j),z0,dt)

      if( BV < wh1) then 

       !print*,j,BV 

        BV=wh1
       
      end if 

      kzh1=dsqrt((BV**2.)*(kxh1**2.)/(wh1**2.)-kxh1**2.)

      xtrih(i,j+1)=xtrih(i,j)+(zdimrh(i)/(il-1))*kxh1/dsqrt(kxh1**2.+kzh1**2.)

      ztrih(i,j+1)=ztrih(i,j)+(zdimrh(i)/(il-1))*kzh1/dsqrt(kxh1**2.+kzh1**2.)

      call GroupVel(Cgxtrih(i,j+1),Cgztrih(i,j+1),N1,ztrih(i,j+1),z0,dt,r,kx,wh1)

   end do 


 end do

open(400,file='HTriPath.dat',status='unknown')

do i=1,il

  write(400,*) xtrih(1,i),ztrih(1,i)
  !print*,N(i)
end do


open(405,file='GroupVelH.dat',status='unknown')

do i=1,il

  write(405,*) Cgxtrih(5,i),Cgztrih(5,i)
  !print*,N(i)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
! Interpolation in Triangular region for first harmonics   !
!                                                          !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(Nptrih(Ntrih))

allocate(IRTrih(Ntrih,il,nxr*nzr))

allocate(PointTrih(Ntrih,nxr*nzr,2))

allocate(PressIntTrih(Num,Ntrih,il))

allocate(UIntTrih(Num,Ntrih,il))

allocate(VIntTrih(Num,Ntrih,il))

allocate(EnergyTrih(Ntrih))

!allocate(EnergyTotTrih(Ntri))

dx=10.*abs(xtrih(1,1)-xtrih(1,2))

dz=10.*abs(ztrih(1,1)-ztrih(1,2))

print*,"dx and dz",dx,dz
do i=1,Ntrih


!here we find the index of the points going to be used in interpolation
 call Points(x1,z1,PointTrih(i,:,:),xtrih(i,:),ztrih(i,:),il,nzr,nxr,Nptrih(i),dx,dz)

print*,Nptrih(i)


!allocate them but it will change in every iteration 
allocate(xi(Nptrih(i)))

allocate(zi(Nptrih(i)))


   do j=1,Nptrih(i)

     xi(j)=x1(PointTrih(i,j,1),PointTrih(i,j,2))

     zi(j)=z1(PointTrih(i,j,1),PointTrih(i,j,2))
 ! print*,xi(j),zi(j)

   end do 

 !print*,"Okey" ,Np(i)   

 call interpolation2D(IRTrih(i,:,:),xtrih(i,:),xi,ztrih(i,:),zi,Nptrih(i),il)


  deallocate(xi)

  deallocate(zi)

end do

print*,"okey trih"
!print*,IRTri(1,:,1)

do i=1,Ntrih

  do j=1,1
 

      call interpVal(IRTrih(i,:,:),PressIntTrih(j,i,:),Press(j,:),PointTrih(i,:,1),PointTrih(i,:,2),nxr,nzr,Nptrih(i),il)
   
      call interpVal(IRTrih(i,:,:),UIntTrih(j,i,:),u(j,:),PointTrih(i,:,1),PointTrih(i,:,2),nxr,nzr,Nptrih(i),il)

      call interpVal(IRTrih(i,:,:),VIntTrih(j,i,:),v(j,:),PointTrih(i,:,1),PointTrih(i,:,2),nxr,nzr,Nptrih(i),il)

  
  end do

end do


open(402,file='TriPressH.dat',status='unknown')

do i=1,il

  write(402,*) (i-1)*zdimrh(1)/il,PressIntTrih(1,5,i)
  !print*,N(i)
end do

open(403,file='TriUH.dat',status='unknown')

do i=1,il

  write(403,*) (i-1)*zdimrh(1)/il,UIntTrih(1,5,i)
  !print*,N(i)
end do

open(404,file='TriVH.dat',status='unknown')

do i=1,il

  write(404,*) (i-1)*zdimrh(1)/il,VIntTrih(1,5,i)
  !print*,N(i)
end do

open(405,file='GroupVelH.dat',status='unknown')

do i=1,il

  write(405,*) Cgxtrih(5,i),Cgztrih(5,i)
  !print*,N(i)
end do

print*,"okey trih 1"

! the arclength starting from the entry of triangular region 

arctrih=0.

do i=2,Ntrih
 
 arctrih(i)=arctrih(i-1)+dsqrt((xtrih(i,1)-xtrih(i-1,1))**2.+(ztrih(i,1)-ztrih(i-1,1))**2.)

end do 

do i=1,Ntrih

 call PressVelocityInteg(EnergyTrih(i),PressIntTrih(1,i,:),UIntTrih(1,i,:),VIntTrih(1,i,:),il,Cgxtrih(i,:),Cgztrih(i,:),zdimrh(i))


end do 

open(406,file='TriEnergyH.dat',status='unknown')

do i=1,Ntrih

  write(406,*) arctrih(i),EnergyTrih(i)
  !print*,N(i)
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						 !
!  Triangular Reflection Region Second Harmonic  !
!						 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Ntrih is the number of points of first harmonic in triangular region 

ds=0.005


do i=1,1000

 if (ds*i > x(N(1),1)-x(N(il),il)) then

  Ntrih2=i-1

  goto 901
  end if


end do

901 continue

print*,"The Number of Points for Second Harmonics",Ntrih2


 allocate(xtrih2(Ntrih2,il))

 allocate(ztrih2(Ntrih2,il))


 allocate(zdimrh2(Ntrih2))

 allocate(Cgxtrih2(Ntrih,il))

 allocate(Cgztrih2(Ntrih,il))

 ds=0.0001


! here we are  changing  kx keep that in mind
kxh2=3.*kx

wh2=3.*w0

do i=1,Ntrih2

 zdimrh2(i)=0.

 !x(N(il)+i,1)

 zdum=z(N(il),il)


 ! in this loop , we find the length of the integration path 
  do j=1,2000 
 
   zdum=zdum+ds

   if (i==1 .and. j==1) then 

     print*,zdum

   end if 
   call BVprofile(BV,N1,r,zdum,z0,dt)

    if(zdum>z(N(il),il)+h .or. BV/wh2<1.015) then 

    goto 902

    end if 

    
     
     
      kzh2=dsqrt((BV**2.)*(kxh2**2.)/(wh2**2.)-kxh2**2.)

    !print*,"kzh2",kzh2,N1,BV,wh2
     zdimrh2(i)=zdimrh2(i)+ds*dsqrt(1.+(kxh2/kzh2)**2.)
   !if(i==1) then
   !print*,zdimr(i),j
   !end if 


   end do

   902 continue 
    
    print*,zdimrh2(i),z(N(il),il)+h,zdum
   
    xtrih2(i,1)=x(N(il),il)+i*ds

    ztrih2(i,1)=z(N(il),il)
  
    call GroupVel(Cgxtrih2(i,1),Cgztrih2(i,1),N1,ztrih2(i,1),z0,dt,r,kxh2,wh2)

   do j=1,il-1

      call BVprofile(BV,N1,r,ztrih2(i,j),z0,dt)

      if( BV < wh2) then 

       !print*,j,BV 

        BV=wh2
       
      end if 

      kzh2=dsqrt((BV**2.)*(kxh2**2.)/(wh2**2.)-kxh2**2.)

      xtrih2(i,j+1)=xtrih2(i,j)+(zdimrh2(i)/(il-1))*kxh2/dsqrt(kxh2**2.+kzh2**2.)

      ztrih2(i,j+1)=ztrih2(i,j)+(zdimrh2(i)/(il-1))*kzh2/dsqrt(kxh2**2.+kzh2**2.)

      call GroupVel(Cgxtrih2(i,j+1),Cgztrih2(i,j+1),N1,ztrih2(i,j+1),z0,dt,r,kxh2,wh2)

   end do 


 end do

open(500,file='H2TriPath.dat',status='unknown')

j=1

do i=1,il

  write(500,*) xtrih2(j,i),ztrih2(j,i)
  !print*,N(i)
end do

open(600,file='FPointline.dat',status='unknown')


do i=1,Ntrih2

  write(600,*) xtrih2(i,1),ztrih2(i,1)
  !print*,N(i)
end do

open(505,file='GroupVelH2.dat',status='unknown')

do i=1,il

  write(505,*) Cgxtrih2(5,i),Cgztrih2(5,i)
  !print*,N(i)
end do



