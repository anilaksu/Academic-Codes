subroutine MeanPU(PU,P,U,V,Cgx,Cgz,N)
!5 Point Gauss Quadrature so N = 5
 integer i
! N is used to allocate the number of points in integration
integer , intent(in)::N
!the group velocity at that point !
real*8,intent(in)::Cgx,Cgz
! the averaged pressure veloctiy product 
real*8,intent(inout)::PU
!pressure and velocity arrays
real*8 ,dimension(N),intent(inout)::P,U,V
!dummy variable
real*8  dummy,dummy1


dummy=0. 

dummy1=dsqrt(Cgx**2.+Cgz**2.)


do i=1,N

  dummy=dummy+P(i)*(U(i)*Cgx+V(i)*Cgz)/dummy1


end do 


PU=dummy/N

return 
end subroutine MeanPU

subroutine QuadPoints(Qx,Qz,kx,kz,xcen,zcen,zdim,N)
!5 Point Gauss Quadrature so N = 5
 integer i
! N is used to allocate the number of points in integration
integer , intent(inout)::N
!the coordinate of the center of the internal wave !
real*8,intent(inout)::xcen,zcen
!the wave numbers of the initial wave and the length of the integration  
real*8,intent(inout)::kx,kz,zdim
!The coordinates of Integration points!
real*8 ,dimension(N),intent(inout)::Qx,Qz
!the angle of the beam , the step size between two adjacent points
real*8  theta ,dl
!pi number
real*8 pi

dl=2.*zdim/(N-1)

pi=4.*datan(1.d0)

if( kx==0.) then
theta=2*datan(1.d0)
goto 100
end if

theta= datan(-kz/kx)

100 continue



do i=1,N

Qx(i)=xcen+(zdim-(i-1)*dl)*dcos(theta)

Qz(i)=zcen-(zdim-(i-1)*dl)*dsin(theta)



end do 


return 
end subroutine QuadPoints
subroutine PressVelocityInteg(Etot,P,U,V,N,Cgx,Cgz,zdim)
!this is Chebyshev point integrator
 integer i,j
! N is used to allocate the number of points in integration
integer , intent(in)::N
! funtion itself and its derivatives 
real*8 ,dimension(N),intent(inout)::P,U,V
! the wave lenghth 
real*8,dimension(N),intent(in)::Cgx,Cgz
!total integrated energy and the integration length
real*8 ,intent(inout)::Etot,zdim
!the group velocity angle and pi number
real*8  theta ,pi
!the magnitude of the group velocity vector
real*8 Cg
! funtion itself and its derivatives 
real*8 ,dimension(N)::E
 

!pi=4.*datan(1.d0)

!if( kx==0.) then
!theta=2*datan(1.d0)
!goto 100
!end if

!theta= datan(-kx/kz)

!100 continue

do i=1,N

  Cg=dsqrt(Cgx(i)**2.+Cgz(i)**2.)

  E(i)=P(i)*(U(i)*Cgx(i)+V(i)*Cgz(i))/Cg

end do 

 call RomBerg(E,Etot,zdim,N)

end subroutine PressVelocityInteg
subroutine RomBerg(f,f1,zdim,N)
!this is Chebyshev point integrator
 integer i,j
! N is used to allocate the number of points in integration
integer , intent(in)::N
! funtion itself and its derivatives 
real*8 ,dimension(N),intent(inout)::f
! the resultant integral 
real*8 ,intent(inout)::f1
!pi number and dummy variable
real*8 pi,dummy
! the dimension
! the normalized coordinates of the chebyshev points
real*8,allocatable:: x(:) 
! step sizes 
real*8 dl1,dl2,dl3
! the idea is to invert derivative matrix to integrate the function 
real*8 fi1,fi2,fi3
! the Number of steps in integration
integer N1,N2,N3 

pi=4.*datan(1.d0)


dl3=8.*zdim/(N-1)

dl2=4.*zdim/(N-1)

dl1=2.*zdim/(N-1)

N1=N-1

N2=N1/2

N3=N2/2

dummy=0.


! Romberg integration steps 
do i=1,N3

  fi3=dummy+0.5*dl3*(f(1+4*i)+f(4*i-3))

  dummy=fi3

end do

dummy=0.

do i=1,N2

  fi2=dummy+0.5*dl2*(f(1+2*i)+f(2*i-1))

  dummy=fi2

end do

dummy=0.

do i=1,N1

  fi1=dummy+0.5*dl1*(f(1+i)+f(i))

  dummy=fi1

end do

fi2=(4.*fi2-fi3)/3.

f1=(4.*fi1-fi2)/3.

return 
end subroutine RomBerg

subroutine Points(x,z,Point,xq,zq,N,nzr,nxr,N1,dx,dz)

!this routine is used to find the interpolation points
integer i,j,k
!Note that l parameter is used to keep the dimension of the points
integer l
!the dimensions of the all the points in the domain
integer ,intent(in)::nzr,nxr
!N number of the quadrature points !
integer , intent(inout)::N
!N1 the number of the grid points will be used in interpolation 
integer ,intent(inout)::N1
! dz ,dx those are the distances determining the interpolation points  
real*8,intent(in)::dx,dz
! the grid points 
real*8 ,dimension(nxr,nzr),intent(inout)::x,z
! the quadrature points
real*8 ,dimension(N),intent(inout)::xq,zq
! the array keeping the index of the points
integer ,dimension(nxr*nzr,2),intent(inout)::Point

! to initiate the l from 1
l=1




do i=1 ,N
 
do j=1,nxr

do k=1,nzr
 
 
  if (abs(xq(i)-x(j,k))<dx .and. abs(zq(i)-z(j,k))<dz ) then

   Point(l,1)=j
   Point(l,2)=k
   l=l+1
  
  end if
  
  
end do 
  
end do

end do  



! now we may the same point more than one , let's clear them
  
 call eliminate(Point,l,nxr,nzr,x,z)

N1=l



end subroutine Points

subroutine eliminate(Point,l,nxr,nzr,x,z)

integer i,j,k,m,n
!Note that l parameter is used to keep the dimension of the points
integer ,intent(inout)::l
! the maximum possible dimensions
integer,intent(in)::nxr,nzr
! the grid points 
real*8 ,dimension(nxr,nzr),intent(inout)::x,z
!the dimensions of the all the points in the domain
integer ,dimension(nxr*nzr,2),intent(inout)::Point
! the error margin for elimination 
real*8 err
! the L2 norm of the distance
real*8 L2norm
err=10E-6
m=1


do n=1,l

 do i=1,l-1
 
   do j=i+1,l-m
L2norm=dsqrt((x(Point(i,1),Point(i,2))-x(Point(j,1),Point(j,2)))**2.+(z(Point(i,1),Point(i,2))-z(Point(j,1),Point(j,2)))**2.)
 if (L2norm<err) then
 if(j<l-m+1) then
 
 m=m+1
 
 end if 
   do k=j,l-1
    Point(k,1)=Point(k+1,1)
    Point(k,2)=Point(k+1,2)
   end do 
  end if 

  end do

 end do 

end do
l=l-m-1
end subroutine eliminate


