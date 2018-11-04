
subroutine NumPointsPyc(N,z0,zcen,nzr,h,kx,kz,r,dt,ds,N1,w0)

!this routine is used to find the number of points to pycnocline 
integer i,p1
!all z points 
integer ,intent(in)::nzr
!N number of the data points out !
integer , intent(inout)::N
!the coordinate of the starting point ,pycnocline center,pycnocline thickness!
real*8,intent(in)::zcen,z0,h
!the z coordinate of the integration path
real*8 z1
! the error margin
real*8 eps
! the group velocities
real*8 Cgx,Cgz
! the pycnocline parameters the measure of pycnocline ,r the strength ratio,Nmin
real*8,intent(in) ::dt,r,N1
! ds integration time difference 
real*8,intent(in)::ds
! the wave frequency
real*8 ,intent(in)::w0
eps=10E-3

z1=zcen



do i=1 ,nzr
 

 call GroupVel(Cgx,Cgz,N1,z,z0,dt,r,kx,w0)


z1=z1+Cgz*ds

  if (abs(z0-z1-h/2)<eps) then
   p1=i
  end if
  
  
  end do 
  
   N=p1
  
end subroutine NumPointsPyc

subroutine Path(Cgx,Cgz,x,z,xcen,zcen,N1,N,kx,ds,r,dt,z0,w0)

 integer i
!dz is the distance between two points in z axis 
real*8 dz
!N number of the data points out !
integer , intent(in)::N
!the coordinate of the center of the internal wave !
real*8,intent(in)::xcen,zcen
!the wave numbers of the initial wave
real*8,intent(in)::kx
!The coordinates of the BeamPath!
real*8 ,dimension(N),intent(inout)::x,z
!The time integration step size
real*8 ,intent(in)::ds
!the group velocities
real*8 ,dimension(N),intent(inout)::Cgx,Cgz
!the bv profile parameters
real*8 ,intent(in)::r,dt,N1,z0
!the frequency of the wave
real*8 ,intent(in)::w0



  



x(1)=xcen

z(1)=zcen 



!lets calculate group velocities to calculate the integration path 

call GroupVel(Cgx(1),Cgz(1),N1,z(1),z0,dt,r,kx,w0)

do i=2,N
!print*,Cgx(i-1),Cgz(i-1),z(i-1)

z(i)=z(i-1)+ds*Cgz(i-1)


x(i)=x(i-1)+ds*Cgx(i-1)


call GroupVel(Cgx(i),Cgz(i),N1,z(i),z0,dt,r,kx,w0)

 !print*,i,Cgx(i),Cgz(i),x(i),z(i)
end do 

return 
end subroutine Path





subroutine bisection(x,r)

 ! this is bisection code to find the pycnocline thickness 
    implicit none
   
   integer i,numit
   real*8 xl,xr,xm,xm1,fl,fr,fm
   real*8 ,intent(inout) :: x,r 
   ! this is left limit for your profile and it is supposed to be negative due to how we defined coordinate system for BV Profile 
   xl=-10.

   !this is right limit 
   xr=0.
   !Number of iteration 
   numit=1000

   xm1=xl

  
   do i=1,numit
    xm=(xl+xr)/2.

   call funcbisec (fl,xl,r)

   call funcbisec (fr,xr,r)

   call funcbisec (fm,xm,r)
   
 ! print*,fl,fr,fm,xl,xr,xm

   if(fl*fr==0)then
   print*,xm1
   exit
   end if
  

    if(fm*fr<0.)then
   xl=xm
   goto 100
   end if
   if(fl*fm<0.)then
   xr=xm
   goto 100
   end if
    if((xm-xm1)<1E-10)then
    ! xm1=xm
   !goto 200
   end if
  if((xm1-xm)<0.000000000001)then
   
   exit
   end if
   100 continue
   xm1=xm
  
   end do

  ! Bingo you find your root
  ! 200 continue
     x=xm1 
   ! print*,"Number of iteration",i,abs(xm-xm1)
  
 end subroutine bisection
    
 
    subroutine funcbisec (f,x,r)
   real*8 :: f
   real*8 :: a,b
   real*8 ,intent(in)::x,r
    a=1.-0.5/(r**2)
   b=1.-a 
  
  f=a*(1./(cosh(x)**2.))+b*(1.-tanh(x))-((r+1.)/(2.*r))**2.
   !a=r**2-0.5
  !f=a*(1/(cosh(x)**2))+0.5*(1-tanh(x))-((r+1)/2)**2
  end subroutine funcbisec

   subroutine BVprofile(N,N1,r,z,z0,dt)

  
  real*8 ,intent(inout) :: N,N1,r,z,z0,dt 
  real*8 ::a,arg1,arg2
  !r = Nmax/N1
  !z0=pycnocline center
  !dt measure of pycnocline thickness

   a=r**2-0.5

   arg1=a*(1/(cosh((4./0.45)*(z-z0)/dt)**2))

   arg2=0.5*(1-tanh((4./0.45)*(z-z0)/dt))
  N=N1*sqrt(arg1+arg2)
  


  end subroutine BVprofile 

   subroutine GroupVel(Cgx,Cgz,N1,z,z0,dt,r,kx,w0)

  
  real*8 ,intent(inout) :: N1,r,z,z0,kx
  real*8 ,intent(inout) :: dt,Cgx,Cgz,w0
  real*8 ::a,kz
  real*8 ::N
  !r = Nmax/N1
  !z0=pycnocline center
  !dt measure of pycnocline thickness
  ! N is BV frequency at that location

  call BVprofile(N,N1,r,z,z0,dt)

   if (N/w0< 1.) then
     print*,"ration" ,N/w0

   end if 
   kz=-1.*dsqrt((N/w0)**2.-1)*kx
! print*,N,N1,r,z,z0,dt
  !group velocity in x direction 
  Cgx=N/dsqrt(kx**2.+kz**2.)-N*(kx**2.)/((kx**2.+kz**2.)**1.5)
  !group velocity in y direction
  Cgz=-N*kx*kz/((kx**2.+kz**2.)**1.5)
 


  end subroutine GroupVel



