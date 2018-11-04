subroutine DisInteg(visloss,visave,Jac,Nx,Nz,dx,dz)

 ! this routine generated to integrate the viscous dissipation
integer i,j,k

!the total number of points in solution domain
integer,intent(in)::Nx,Nz

! step size
real*8,intent(inout)::dx,dz

! viscous loss
real*8,dimension(Nx),intent(inout)::visloss

! the shifted coordinate
real*8,dimension(Nx,Nz),intent(inout)::Jac,visave

! viscous loss integrated along and average viscous loss divided by Jacobian at that point
real*8,allocatable::visz(:),visJ(:,:)

! dummy variable length to be used in Romberg integration
real*8 dummy,zdiml

! since our coordinate transformation we unit cells 

zdiml=(Nz-1)*0.5*dz

dx=1.

allocate(visz(Nx))

allocate(visJ(Nx,Nz))

 ! let's generate visJ

  do i=1,Nx
   do j=1,Nz
  
    visJ(i,j)=visave(i,j)/Jac(i,j)

   end do
  end do
  

! let's integrate it along z 
 
 
   do i=1,Nx
 
     call RomBerg(visJ(i,:),visz(i),zdiml,Nz)
 
   end do

! now let's integrate along x

  visloss=0.
  do i=1,Nx-1

   visloss(i+1)=visloss(i)+0.5*dx*(visz(i)+visz(i+1))    

  end do
  print*,"okey Integ"

end subroutine DisInteg

subroutine VisAver(visave,vis,Num,Nx,Nz)

 ! this routine generated to take average of the viscous loss
integer i,j,k

!the total number of time steps
integer,intent(in)::Num,Nx,Nz

! viscous loss time series
real*8,dimension(Num,Nx,Nz),intent(inout)::vis


! viscous loss time series
real*8,dimension(Nx,Nz),intent(inout)::visave

! dummy variable length to be used in Romberg integration
real*8 dummy




  open(268,file='VisAve.dat',status='unknown')


 ! let's calculate time average

  do i=1,Nx

   do j=1,Nz
   
      dummy=0.

     do k=1,Num
   
      dummy=dummy+vis(k,i,j)

     end do

      visave(i,j)=dummy/Num

   write(268,*)  i,j,visave(i,j)
  
   end do

  end do

print*,"okey aver"
  
end subroutine VisAver



subroutine VisDis(vis,Dx,Dz,u,v,Nx,Nz,xnu)

 ! this routine generated to calculate the viscous dissiantip

integer i,j,k

!the total number of points in solution domain
integer,intent(in)::Nx,Nz

! kinematic viscosity 
real*8,intent(in)::xnu

! the shifted coordinate
real*8,dimension(Nx,Nz),intent(inout)::u,v,vis

! the shifted coordinate
real*8,dimension(Nx*Nz,Nx*Nz),intent(in)::Dx,Dz

! 1D vector for x1 and z1 and their derivative with respect to x and z
real*8,allocatable::u1(:),v1(:),vis1(:),ux(:),uz(:),vx(:),vz(:)


! total number of points
integer Nt

Nt=Nx*Nz

allocate(u1(Nt))

allocate(v1(Nt))

allocate(vis1(Nt))

allocate(ux(Nt))

allocate(uz(Nt))

allocate(vx(Nt))

allocate(vz(Nt))

 do j=1,Nz
   do i=1,Nx

    k=(j-1)*Nx+i

    u1(k)=u(i,j)

    v1(k)=v(i,j)
  

   end do
 end do

 ! let's evaluate dx1/dx
 call matvect(Dx,u1,ux,Nt)

  ! let's evaluate dx1/dz
 call matvect(Dz,u1,uz,Nt)

  ! let's evaluate dz1/dx
 call matvect(Dx,v1,vx,Nt)

  ! let's evaluate dz1/dz
 call matvect(Dz,v1,vz,Nt)

 
  do i=1,Nt

   ! let's calculate viscous dissipation
   vis1(i)=2.*xnu*(ux(i)**2.+0.5*(uz(i)+vx(i))**2.+vz(i)**2.)
   !  print*,Jac1(i)
  
  end do

  !let's return it in two dimensional array 
   do j=1,Nz
     do i=1,Nx

     k=(j-1)*Nx+i

     vis(i,j)=vis1(k)
  
     end do
   end do

print*," okey VisDis"

end subroutine VisDis


subroutine Jacobian(Jac,Dx,Dz,x1,z1,Nx,Nz)

integer i,j,k

!the total number of points in solution domain
integer,intent(in)::Nx,Nz


! the shifted coordinate
real*8,dimension(Nx,Nz),intent(inout)::x1,z1,Jac

! the shifted coordinate
real*8,dimension(Nx*Nz,Nx*Nz),intent(in)::Dx,Dz

! 1D vector for x1 and z1 and their derivative with respect to x and z
real*8,allocatable::x11(:),z11(:),Jac1(:),x1x(:),x1z(:),z1x(:),z1z(:)


! total number of points
integer Nt

Nt=Nx*Nz

allocate(x11(Nt))

allocate(z11(Nt))

allocate(Jac1(Nt))

allocate(x1x(Nt))

allocate(x1z(Nt))

allocate(z1x(Nt))

allocate(z1z(Nt))

 do j=1,Nz
   do i=1,Nx

    k=(j-1)*Nx+i

    x11(k)=x1(i,j)

    z11(k)=z1(i,j)
  

   end do
 end do

 ! let's evaluate dx1/dx
 call matvect(Dx,x11,x1x,Nt)

  ! let's evaluate dx1/dz
 call matvect(Dz,x11,x1z,Nt)

  ! let's evaluate dz1/dx
 call matvect(Dx,z11,z1x,Nt)

  ! let's evaluate dz1/dz
 call matvect(Dz,z11,z1z,Nt)

 
  do i=1,Nt

   ! let's calculate jacobian
   Jac1(i)=dabs(x1x(i)*z1z(i)-x1z(i)*z1x(i))
   !  print*,Jac1(i)
  
  end do
open(267,file='Jaco.dat',status='unknown')

   

  !let's return it in two dimensional array 
   do j=1,Nz
     do i=1,Nx

     k=(j-1)*Nx+i

     Jac(i,j)=Jac1(k)
   
     write(267,*)  i,j,Jac(i,j)

     end do
   end do

end subroutine Jacobian

subroutine CoorShift(x,z,x1,z1,Nx,Nz)

integer i,j,k

!the total number of points in solution domain
integer,intent(in)::Nx,Nz


! the original coordinate
real*8,dimension(Nx,Nz),intent(in)::x,z

! the shifted coordinate
real*8,dimension(nx,Nz),intent(inout)::x1,z1

x1=0.

z1=0.


 do i=1,Nx
   do j=1,Nz

     x1(i,j)=(i-1)*1.
     
     z1(i,j)=(j-1)*1.
     !if (j>1) then
     ! j indice keeps track of the new z1
     !z1(i,j)=z1(i,j-1)+dsqrt((x(i,j)-x(i,j-1))**2.+(z(i,j)-z(i,j-1))**2.)    
     ! end if

     ! if(i>1) then
     ! i indice keeps track of the new x1
     !x1(i,j)=x1(i-1,j)+dsqrt((x(i-1,j)-x(i,j))**2.+(z(i,j)-z(i-1,j))**2.)
     ! end if
     print*,x1(i,j),z1(i,j)
   end do
 end do


end subroutine CoorShift

subroutine dfdx(Afx,f,fx,Nx,Nz)

integer i,j,k

!the total number of points in solution domain
integer,intent(in)::Nx,Nz


! the original coordinate
real*8,dimension(Nx,Nz),intent(inout)::f,fx


! the shifted coordinate
real*8,dimension(Nx*Nz,Nx*Nz),intent(inout)::Afx

! total number of points
 integer Nt

! 1D vector
real*8,allocatable::f1(:),fx1(:)


Nt=Nx*Nz

 allocate(f1(Nt))

 allocate(fx1(Nt))


 do j=1,Nz

   do i=1,Nx

    k=(j-1)*Nx+i

     f1(k)=f(i,j)      
 
   end do

 end do

 call matvect(Afx,f1,fx1,Nt)

 ! now let's write the derivative again in two dimension

 do j=1,Nz

   do i=1,Nx

    k=(j-1)*Nx+i

     fx(i,j)=fx1(k)      
 
   end do

 end do

end subroutine dfdx

 subroutine Diffx(Afx,x,z,Nx,Nz)

integer i,j,k

!the total number of points in solution domain
integer,intent(in)::Nx,Nz


! the original coordinate
real*8,dimension(Nx,Nz),intent(in)::x,z

! the shifted coordinate
real*8,dimension(Nx*Nz,Nx*Nz),intent(inout)::Afx

! total number of points
 integer Nt

! 1D vector
real*8,allocatable::x1(:),z1(:)

! interpolation matrix and its inverse and diff matrix
real*8,allocatable::A(:,:),A1(:,:),Ax(:,:)
Nt=Nx*Nz
!  now let's store all values in 1D vector to use RBF

 allocate(x1(Nt))

 allocate(z1(Nt))

allocate(A(Nt,Nt))

allocate(A1(Nt,Nt))

allocate(Ax(Nt,Nt))


 do j=1,Nz

   do i=1,Nx

    k=(j-1)*Nx+i

     x1(k)=x(i,j)    
  
     z1(k)=z(i,j)    
 
   end do

 end do


 ! let's generate the interpolation matrix
 Call matgen(A,x1,z1,Nt)
 
  print*,Nt
 ! let's invert it
 Call invertSVD(Nt,A,A1)

 ! let's generate differentiation matrix
 Call matgenx(Ax,x1,z1,Nt)  
  
  ! the resultant matrix 
 CALL matmultnon(Ax,A1,Afx,Nt,Nt,Nt)

 
end subroutine Diffx




subroutine matgenx(A,x,y,N)


integer , intent(in)::N

real*8 ,dimension(N,N),intent(inout)::A

real*8 ,dimension(N),intent(inout)::x,y

integer i,j,k

real*8 c,mqx

 c=2E-10

do i=1,N
do j=1,N
! the famous rbf ! 
 ! each entry of the matrix is the values of the rbf !
 ! by summing the effect of the each entry multiplied by corresponding alpha coefficent , the resultant interploation is obtained  
 call multix(mqx,x(i),x(j),y(i),y(j),c)

A(i,j)=mqx



end do 
end do

return 
end subroutine matgenx


 subroutine multix(mqx,xi,xj,yi,yj,c)
 
  ! the famous radial basis functions derivative with respec to x !  
   real*8 :: mqx
   real*8 ,intent(in)::xi,xj,yi,yj,c
  
  mqx=(xi-xj)/sqrt((xi-xj)**2.+(yi-yj)**2.+c*c)
  !r=dsqrt((xi-xj)**2.+(yi-yj)**2.)
  
 !  mq=dlog(r)*(r**2.)

  end subroutine multix
 

 
