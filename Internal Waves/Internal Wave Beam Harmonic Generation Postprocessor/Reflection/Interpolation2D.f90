subroutine TimeInter(Press1,UVel1,VVel1,Press,UVel,VVel,time,Np,Num,nxr,nzr,dt)

integer i,j,k

!the total number of points in solution domain
integer,intent(in)::nxr,nzr

! the number of data point in time 
integer,intent(in)::Np

! time array
real*8,dimension(Np),intent(in)::time

! the original pressure and velocity data 
real*8,dimension(Np,nxr*nzr),intent(in)::Press,UVel,VVel

! the original pressure and velocity data 
real*8,dimension(Num,nxr*nzr),intent(inout)::Press1,UVel1,VVel1

! incremental time difference
real*8 ,intent(inout)::dt

! the interpolation time array and dummy time used to imitate y axis in interpolation
real*8 ,allocatable::IntTime(:),timdum1(:),timdum2(:)

! total time difference
real*8 ttime



! the matrces to be used in interpolation
real*8,allocatable::A(:,:),A1(:,:),Ai(:,:),Ar(:,:)

allocate(IntTime(Num))

allocate(timdum1(Np))

allocate(timdum2(Num))

allocate(A(Np,Np))

allocate(A1(Np,Np))

allocate(Ai(Num,Np))

allocate(Ar(Num,Np))

ttime=time(Np)-time(1)

dt=ttime/(Num-1)


do i=1,Num 

 IntTime(i)=time(1)+dt*(i-1)
! print*,IntTime(i)
end do


timdum1=0.

timdum2=0.

!this function generates A matrix using x coordinates
 CALL matgen(A,time,timdum1,Np)


 !this function generates Ai matrix using xi coordinates
 CALL matgenxi(Ai,IntTime,time,timdum2,timdum1,Num,Np)

  

 ! the inversion of the system matrix !
 CALL invertSVD(Np,A,A1)

 

  ! the resultant matrix 
 CALL matmultnon(Ai,A1,Ar,Num,Np,Np)
 
 !print*,Ar(1,2)
 do i=1,nxr*nzr

   call matvectnon(Ar,Press(:,i),Press1(:,i),Num,Np)

   call matvectnon(Ar,UVel(:,i),UVel1(:,i),Num,Np)

   call matvectnon(Ar,VVel(:,i),VVel1(:,i),Num,Np)
 
 end do

 
return 
end subroutine TimeInter






subroutine interpVal(A,IntVal,DatVal,Pointx,Pointz,nxr,nzr,Np,N)

!the total number of points in solution domain
integer,intent(in)::nxr,nzr

! the number of points usedin interpolation
integer,intent(in)::Np


! the interpolation matrix 
real*8,dimension(N,Np),intent(inout)::A

! the Interpolated value 
real*8,dimension(N),intent(inout):: IntVal

! the data used to generate the interpolation
real*8 ,dimension(nxr*nzr),intent(inout):: DatVal

! the array keeps track of the data point used in interpolation 
integer ,dimension(Np),intent(in)::Pointx,Pointz


integer i,j,k



real*8 ,allocatable ::DatIn(:)

!Data to be used in interpolation
allocate(DatIn(Np))

 
 !print*,Pointx(1),Pointz(1)

 !print*,A(:,2)

 
do i=1,Np



DatIn(i)=DatVal((Pointz(i)-1)*nxr+Pointx(i))

end do



 call matvectnon(A,DatIn,IntVal,N,Np)

return 
end subroutine interpVal


subroutine interpolation2D(Ar,xi,x,yi,y,N1,N2)

! this routine includes all required routines for the interpolation !

! all you need to do is to give the interpolation data f , x and the coordinates of the resultant interpolation xi !

! you will obtain the resultant interploation as fi !

!N1 is the size of actual data!

!N2 is the size of interpolation data!
integer , intent(in)::N1,N2

!f vector is the value of the actual data !

!x is the coordinate vector of the data !
real*8 ,dimension(N1),intent(inout)::x,y

!fi vector is the value of the interpolated data !

!xi is the coordinate vector of the interpolated data !
real*8 ,dimension(N2),intent(inout)::yi,xi

! Ar the interpolation matrix
real*8,dimension(N2,N1),intent(inout)::Ar

integer i,j,k



real*8 ,allocatable ::A(:,:),A1(:,:),Ai(:,:)



!A matrix is the matrix used to form interpolation coeffiecnts
allocate(A(N1,N1))

!A1 matrix is the inverse of the A matrix by multiply it with f vector , interpolation coeffiecients obtained !
allocate(A1(N1,N1))

!Ai matrix is used to obtain the interpolated data using interpolation coefficients at interpolation locations !
allocate(Ai(N2,N1))


 

 !this function generates A matrix using x coordinates
 CALL matgen(A,x,y,N1)

 !this function generates Ai matrix using xi coordinates
 CALL matgenxi(Ai,xi,x,yi,y,N2,N1)

  

 ! the inversion of the system matrix !
 CALL invertSVD(N1,A,A1)

  ! the resultant matrix 
 CALL matmultnon(Ai,A1,Ar,N2,N1,N1)
  

return 
end subroutine interpolation2D



subroutine matgenxi(A,xi,x,yi,y,N1,N2)


!N1 is the number of the columns !
!N2 is the number of the rows!

integer , intent(in)::N1,N2

real*8 ,dimension(N1,N2),intent(inout)::A

real*8 ,dimension(N1),intent(inout)::x,y

real*8 ,dimension(N2),intent(inout)::xi,yi

integer i,j,k

real*8 c,mq

 c=2E-10

do i=1,N1
do j=1,N2
! the famous rbf ! 
 ! each entry of the matrix is the values of the rbf !
 ! by summing the effect of the each entry multiplied by corresponding alpha coefficent , the resultant interploation is obtained  
 call multi(mq,xi(i),x(j),yi(i),y(j),c)

A(i,j)=mq

end do 
end do

return 
end subroutine matgenxi


subroutine matgen(A,x,y,N)


integer , intent(in)::N

real*8 ,dimension(N,N),intent(inout)::A

real*8 ,dimension(N),intent(inout)::x,y

integer i,j,k

real*8 c,mq

 c=2E-10

do i=1,N
do j=1,N
! the famous rbf ! 
 ! each entry of the matrix is the values of the rbf !
 ! by summing the effect of the each entry multiplied by corresponding alpha coefficent , the resultant interploation is obtained  
 call multi(mq,x(i),x(j),y(i),y(j),c)

A(i,j)=mq



end do 
end do

return 
end subroutine matgen 

subroutine matvect(A,x,b,N)
! simple matrix vector multiplication for square matrix 

integer , intent(in)::N

real*8 ,dimension(N,N),intent(inout)::A

real*8 ,dimension(N),intent(inout)::x,b

integer i,j,k

real*8 c

!c is a dummy variable!

b=0.

do i=1,N
 c=0.

do k=1,N
 c=c+A(i,k)*x(k)
end do 

b(i)=c

end do

return 
end subroutine matvect



subroutine matvectnon(L,U,A,N1,N2)
! simple matrix vector multiplication for non-square matrix 

!For non-square matrix 

integer , intent(in)::N1,N2

real*8 ,dimension(N1,N2),intent(inout)::L

real*8 ,dimension(N2),intent(inout)::U

real*8 ,dimension(N1),intent(inout)::A

integer i,j,k

real*8 c

!c is a dummy variable!

A=0.

do i=1,N1

 c=0.

do k=1,N2
 c=c+L(i,k)*U(k)
end do 

A(i)=c



end do

return 
end subroutine matvectnon

subroutine matmultnon(L,U,A,N1,N2,N3)
! simple matrix vector multiplication for non-square matrix 

!For non-square matrix 

integer , intent(in)::N1,N2,N3

real*8 ,dimension(N1,N2),intent(inout)::L

real*8 ,dimension(N2,N3),intent(inout)::U

real*8 ,dimension(N1,N3),intent(inout)::A

integer i,j,k

real*8 c
!the vector where the nonrectangular matrix vector products stored
real*8,allocatable::V(:)

allocate(V(N1))
!c is a dummy variable!

A=0.

do i=1,N3

 CALL matvectnon(L,U(:,i),V,N1,N2)

do k=1,N1

 A(k,i)=V(k)

end do 


end do

return 
end subroutine matmultnon

              
   subroutine multi(mq,xi,xj,yi,yj,c)
 
  ! the famous radial basis function !  
   real*8 :: mq
   real*8 ,intent(in)::xi,xj,yi,yj,c
  
  mq=sqrt((xi-xj)**2.+(yi-yj)**2.+c*c)
  !r=dsqrt((xi-xj)**2.+(yi-yj)**2.)
  
 !  mq=dlog(r)*(r**2.)

  end subroutine multi
 
   
 

subroutine invertSVD(n,A,Ainv)
! matrix inversion via singular value decompostion !
  integer :: n
  real*8  :: A(n,n),Ainv(n,n)
  real*8,allocatable :: U(:,:),S(:),Vt(:,:),work(:),dum(:,:)
  integer :: info,i,j,k
  allocate(U(n,n),S(n),Vt(n,n),dum(n,n),work(5*n))
  dum = A
  call dgesvd('A' ,'A'   , n, n, dum,   n, S, U,   n, Vt,    n, work,  5*n ,info)
  do i = 1,n
    do j = 1,n
      Ainv(i,j) = 0.d0
      do k = 1,n
        Ainv(i,j) = Ainv(i,j)+Vt(k,i)*U(j,k)/S(k)
      end do
    end do
  end do
  deallocate(U,S,Vt,work,dum)
end subroutine invertSVD
