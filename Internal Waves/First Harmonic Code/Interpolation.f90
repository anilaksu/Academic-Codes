!!! this file contains all interpolation routines 
!!! First One is my favourite one Radial Basis Function Interpolation
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

subroutine matgen(A,x,y,N)
	! This subroutine generates required RBF solution matrix for interpolation
	integer , intent(in)::N
	! the solution matrix A
	real*8 ,dimension(N,N),intent(inout)::A
	! [A]x=y
	real*8 ,dimension(N),intent(inout)::x,y

	integer i,j
	real*8 c,mq
	! c coefficient used in RBF mq function f
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

subroutine matgenxi(A,xi,x,yi,y,N1,N2)
! this routine is generated to interpolate the data on desired points given as x,y
integer , intent(in)::N1,N2
	!N1 is the number of the columns !
	!N2 is the number of the rows!
	real*8 ,dimension(N1,N2),intent(inout)::A
	! Interpolation points
	real*8 ,dimension(N1),intent(inout)::x,y
	! Data Points coordinates xi, yi
	real*8 ,dimension(N2),intent(inout)::xi,yi

	integer i,j,k
	real*8 c,mq
	! c coefficient used in RBF mq function f
	c=2E-10

	do i=1,N1
		do j=1,N2
			! the famous rbf ! 
			! each entry of the matrix is the values of the rbf !
			! by summing the effect of the each entry multiplied by corresponding alpha coefficent , the resultant interploation is obtained  
			call multi(mq,xi(j),x(i),yi(j),y(i),c)
			A(i,j)=mq
		end do 
	end do

	return 
end subroutine matgenxi

subroutine multi(mq,xi,xj,yi,yj,c) 
	! the famous radial basis function !  
	real*8 :: mq
	real*8 ,intent(in)::xi,xj,yi,yj,c
	mq=sqrt((xi-xj)**2.+(yi-yj)**2.+c*c)
end subroutine multi



subroutine getInterpolationPoints(Indices,counter,x_int,x_data,r,N_int,N_data)
	! This subroutine finds the indices of data points 
	! to be used in interpolation within a given radius
	integer , intent(in)::N_int,N_data
	! Indice array
	integer ,dimension(N_data),intent(inout)::Indices
	! data points
	real*8 ,dimension(N_data,2),intent(inout)::x_data
	! the interpolation points
	real*8 ,dimension(N_int,2),intent(inout)::x_int
	! the radius
	real*8 r,r_dist 
	! the number of points for interpolation
	integer, intent(inout)::counter
	! logical statement to keep track of the indice
	logical ::True
	
	integer i,j
	! the counter to keep track of the indices
	counter=0
	
	do i=1,N_int
		do j=1,N_data
			! the radial distance between data point and the interpolation point
			r_dist=dsqrt((x_data(j,1)-x_int(i,1))**2.+(x_data(j,2)-x_int(i,2))**2.)
			if (r_dist<r) then
				! if it is the first time to fall within a given radius, there is no need to check
				! if it is taken previously
				if (counter==0) then
					counter=1
					Indices(counter)=j
				! now it is time to check if it is taken before
				else
					call checkPoints(Indices,j,counter,True)
					if(True) then
						!print*,counter,j
						counter=counter+1
						Indices(counter)=j
					end if
				end if
			end if
			
		end do
	end do

	return 
end subroutine getInterpolationPoints

subroutine checkPoints(IndiceArray,Indice,N_data,True)
	! This subroutine checks if the incide already taken to indice array or not
	integer , intent(in)::N_data
	! Indice array 
	Integer,dimension(N_data),intent(inout)::IndiceArray
	! the new indice
	integer, intent(inout)::Indice
	! logical statement to keep track of the indice
	logical, intent(inout)::True
	
	integer i
	! the counter to keep track of the indices
	counter=0
	! first assume false, loop will determine it
	True=.true.
	do i=1,N_data
		if(IndiceArray(i)==Indice) then 
			True=.false.
		end if
	end do

	return 
end subroutine checkPoints

subroutine getDataForInterpolation(IndiceArray,totalData,intData,x_data,x_intData,counter,N_data)
	! This subroutine checks if the incide already taken to indice array or not
	integer , intent(in)::N_data,counter
	! Indice array 
	Integer,dimension(counter),intent(in)::IndiceArray
	! total coordinates
	real*8, dimension(N_data,2),intent(in)::x_data
	! the interpolation coordinates
	real*8 ,dimension(counter,2),intent(inout)::x_intData
	! total data
	real*8, dimension(N_data),intent(in)::totalData
	! data to be used in the interpolation 
	real*8 ,dimension(counter),intent(inout)::intData
	
	integer i
	
	do i=1,counter
		x_intData(i,:)=x_data(IndiceArray(i),:)
		intData(i)=totalData(IndiceArray(i))
	end do

	return 
end subroutine getDataForInterpolation
