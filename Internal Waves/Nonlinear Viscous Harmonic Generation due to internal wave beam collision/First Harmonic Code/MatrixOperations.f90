!!! This routine is generated to perform required matrix operations
!!! New operations can be added
subroutine matvect(A,x,y,N)
	! simple matrix vector multiplication for square matrix 
	integer , intent(in)::N
	real*8 ,dimension(N,N),intent(inout)::A
	! [A]x=b
	real*8 ,dimension(N),intent(inout)::x,y

	integer i,j,k
	!c is a dummy variable!
	real*8 c
	! y has to start as zero vector
	y=0.

	do i=1,N
		c=0.
		do k=1,N
			c=c+A(i,k)*x(k)
		end do 
		y(i)=c
	end do

	return 
end subroutine matvect

subroutine cmatvect(A,x,y,N)
	! simple matrix vector multiplication for square matrix 
	integer , intent(in)::N
	complex*8,dimension(N,N),intent(inout)::A
	! [A]x=b
	complex*8 ,dimension(N),intent(inout)::x,y

	integer i,j,k
	!c is a dummy variable!
	complex c
	! y has to start as zero vector
	y=0.

	do i=1,N
		c=0.
		do k=1,N
			c=c+A(i,k)*x(k)
		end do 
		y(i)=c
	end do

	return 
end subroutine cmatvect

subroutine cvectvect(uf,u,f,N)
	! this routine performs vector vector dot product in complex domain
	integer , intent(in)::N
	complex,dimension(N),intent(inout)::u,f,uf
	integer i
	do i=1,N 
		uf(i)=u(i)*f(i)
	end do
	return 
end subroutine cvectvect

subroutine matvectnon(L,x,y,N1,N2)
	! simple matrix vector multiplication for non-square matrix 
	!For non-square matrix 
	! N1 is the number of rows in Matrix L
	! N2 is the number of columns in Matrix L
	integer , intent(in)::N1,N2
	! [L]x=y 
	real*8 ,dimension(N1,N2),intent(inout)::L
	! x is the vector multiplied with matrix L
	real*8 ,dimension(N2),intent(inout)::x
	! y is the resultant vector
	real*8 ,dimension(N1),intent(inout)::y

	integer i,k
	!c is a dummy variable!
	real*8 c
	! y has to start as zero vector
	y=0.

	do i=1,N1
		c=0.
		do k=1,N2
			c=c+L(i,k)*x(k)
		end do 
		y(i)=c
	end do

	return 
end subroutine matvectnon

subroutine matmultnon(L,U,A,N1,N2,N3)
	! simple matrix vector multiplication for non-square matrix 
	!For non-square matrix 
	![L][U]=[A]
	integer , intent(in)::N1,N2,N3
	! N1 is the number of rows in Matrix L
	! N2 is the number of columns in Matrix L
	real*8 ,dimension(N1,N2),intent(inout)::L
	! N2 is the number of rows in Matrix U
	! N3 is the number of columns in Matrix U
	real*8 ,dimension(N2,N3),intent(inout)::U
	! N1 is the number of rows in Matrix A
	! N3 is the number of columns in Matrix A
	real*8 ,dimension(N1,N3),intent(inout)::A

	integer i,j,k
	!c is a dummy variable!
	real*8 c
	!the vector where the nonrectangular matrix vector products stored
	real*8,allocatable::V(:)
	allocate(V(N1))
	! A has to start as zero matrix
	A=0.

	do i=1,N3
		CALL matvectnon(L,U(:,i),V,N1,N2)
		do k=1,N1
			A(k,i)=V(k)
		end do 
	end do

	return 
end subroutine matmultnon

subroutine IdMatrix(Id,N)
	! this routine generates identity matrix
	integer , intent(in)::N
	! N1 is the number of rows in Matrix L
	! N2 is the number of columns in Matrix L
	real*8 ,dimension(N,N),intent(inout)::Id
	integer i
	! first let's set it to zero
	Id=0.
	do i=1,N
		Id(i,i)=1.
	end do

	return 
end subroutine IdMatrix
