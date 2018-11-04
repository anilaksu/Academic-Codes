

subroutine RayPath(Np,z0,zs,xs,kx,r,dt,ds,N1,w0,x,z,Nref,Cgx,Cgz)

!this routine is used to find the number of points to pycnocline 
integer i,j,k
!N number of the data points out and Number of reflections!
integer , intent(inout)::Np,Nref
!the coordinate of the starting point ,pycnocline center,pycnocline thickness!
real*8,intent(in)::zs,xs,z0
!the z coordinate of the integration path
real*8 z1,BV
! the error margin and pi number
real*8 eps,pi
!the group velocities
real*8 ,dimension(Np*Nref),intent(inout)::Cgx,Cgz,x,z
! the pycnocline parameters the measure of pycnocline ,r the strength ratio,Nmin
real*8,intent(in) ::dt,r,N1
! ds integration time difference 
real*8,intent(in)::ds
! the wave frequency
real*8 ,intent(in)::w0
eps=10E-3

pi=4.*datan(1.d0)

z1=zcen


x(1)=xs

z(1)=zs

do j=1,Nref

k=(j-1)*Np

 do i=1,Np-1
 

 call GroupVel(Cgx(i+k),Cgz(i+k),N1,z(i+k),z0,dt,r,kx,w0)
 
! the signum funcion reverses the direction of the ray path 

  Cgz(i+k)=dcos(pi*(j-1))*Cgz(i+k)

  x(i+k+1)=x(i+k)+Cgx(i+k)*ds/abs(Cgx(i+k)**2.+Cgz(i+k)**2.)

  z(i+k+1)=z(i+k)+Cgz(i+k)*ds/abs(Cgx(i+k)**2.+Cgz(i+k)**2.)

 end do 
  
  call GroupVel(Cgx(j*Np),Cgz(j*Np),N1,z(j*Np),z0,dt,r,kx,w0)
 
! the signum funcion reverses the direction of the ray path 

  Cgz(j*Np)=(-1.**(j+1))*Cgz(j*Np)
! to set the reflection point as the last point of the previous ray path 

 x(j*Np+1)=x(j*Np)

 z(j*Np+1)=z(j*Np)

end do    
  
!ds=0.001

!do i=1,1000

 !call BVprofile(BV,N1,r,ztri(i,j),z0,dt)

 !call GroupVel(Cgxtri(i,j+1),Cgztri(i,j+1),N1,ztri(i,j+1),z0,dt,r,kx,w0)


!end do


end subroutine RayPath



subroutine PathPoints(Np,z0,zcen,nzr,kx,r,dt,ds,N1,w0)

!this routine is used to find the number of points to pycnocline 
integer i,p1
!all z points 
integer ,intent(in)::nzr
!N number of the data points out !
integer , intent(inout)::Np
!the coordinate of the starting point ,pycnocline starting point
real*8,intent(in)::zcen,z0
!the z coordinate of the integration path and BV frequency
real*8 z1,BV
! the error margin
real*8 eps
! the group velocities
real*8 Cgx,Cgz
! the pycnocline parameters the measure of pycnocline ,r the strength ratio,Nmin
real*8,intent(in) ::dt,r,N1
! ds integration path increment
real*8,intent(in)::ds
! the wave frequency
real*8 ,intent(in)::w0
! the wave number 
real*8 ,intent(in):: kx
eps=1E-1

z1=zcen




do i=1 ,1000
 
 
 call BVprofile(BV,N1,r,z1,z0,dt)

 !print*,"BV freq",BV
 if(BV/N1-1.<eps) then

   goto 300

 end if

 call GroupVel(Cgx,Cgz,N1,z1,z0,dt,r,kx,w0)

 z1=z1+ds*Cgz/abs(Cgx**2.+Cgz**2.)

! print*,"z1",z1

end do

300 continue 


Np=i-3




end subroutine PathPoints

