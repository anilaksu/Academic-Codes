!!! this subroutine includes all associated integral transforms
!!! These include fourier transform and Legendre Transform
!!! It will be updated


subroutine ComplexFFT(Odata,FourData,N,ssign)
	! This function performs Forward and Backward Fourier Transform
	! If ssign=1, Forward Fourier Transform
	! If ssign=-1, Backward Fourier Transform
	! the size of the data 
	integer ,intent(in):: N,ssign
	! the original data and Filtered Data
	complex*8 ,dimension(N),intent(inout)::FourData,OData
	!! All arrays to use fftmod.f routines
	complex ,allocatable :: A(:,:),WORK(:,:)
	INTEGER IFAX1(13)
	complex ,allocatable :: TRIGS1(:)

	allocate(A(N+2,1))
	allocate(WORK(N,1))
	allocate(TRIGS1(2*N))

	CALL CFTFAX(N, IFAX1, TRIGS1) 


	! let's store the data into FFT matrix
	do i=1,N
	if (ssign==1) then				
			A(i,1)=OData(i)
		else 	
			A(i,1)=FourData(i)
		end if		
	end do

	! here we obtaine the fourier transform of the data
	 CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, N+2, N, 1, -1*ssign)

	! here we store the fourier transform of the data
	do i=1,N
		if (ssign==1) then				
			FourData(i)=A(i,1)
		else 	
			OData(i)=A(i,1)/N
		end if
	end do

end subroutine ComplexFFT


subroutine RealFFT(Odata,FourData,N,ssign)
	! This function performs Forward and Backward Fourier Transform
	! If ssign=1, Forward Fourier Transform
	! If ssign=-1, Backward Fourier Transform
	! the size of the data 
	integer ,intent(in):: N,ssign
	! the original data 
	real*8 ,dimension(N),intent(inout)::OData
	! the fourier transformed data
	complex*8 ,dimension(N),intent(inout)::FourData
	!! All arrays to use fftmod.f routines
	complex ,allocatable :: A(:,:),WORK(:,:)
	INTEGER IFAX1(13)
	complex ,allocatable :: TRIGS1(:)

	allocate(A(N+2,1))
	allocate(WORK(N,1))
	allocate(TRIGS1(2*N))

	CALL CFTFAX(N, IFAX1, TRIGS1) 


	! let's store the data into FFT matrix
	do i=1,N
	if (ssign==1) then				
			A(i,1)=CMPLX(OData(i),0)
		else 	
			A(i,1)=FourData(i)
		end if		
	end do
	! here we obtaine the fourier transform of the data
	 CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, N+2, N, 1, -1*ssign)

	! here we store the fourier transform of the data
	do i=1,N
	! Note that in real fourier transform first half of the transfrom coefficient are used
		if (ssign==1 ) then
			if (i>1 .and. i<=N/2) then
				FourData(i)=2.*A(i,1)
			else if(i==1 .or. i==(N/2+1) ) then
					FourData(i)=A(i,1)
				else 
					FourData(i)=0.
			end if
		else 	
			OData(i)=real(A(i,1)/N)
		end if
	end do
end subroutine RealFFT

subroutine LegendreTransform(Odata,LegenData,GLLpoints,w,N)
	! This function performs Forward Legendre Transform
	! the size of the data 
	integer ,intent(in):: N
	! the original data and transformed Data
	real*8 ,dimension(N),intent(inout)::OData
	! the original data and transformed Data
	real*8 ,dimension(N+1),intent(inout)::LegenData
	! the grid points and corresponding weights
	real*8 ,dimension(N),intent(in)::GLLpoints,w
	integer i,j
    !! Legendre Polynomial associated arrays
	real*8 ,allocatable :: Ln(:,:),Lpn(:,:),q(:,:),qp(:,:)
	
	allocate(Ln(N+1,N))
	allocate(Lpn(N+1,N))
	allocate(q(N+1,N))
	allocate(qp(N+1,N))
	
	! let's generate all required Legendre polynomials for Legendre Transform
	do i=1,(N+1)
		do j=1,N
			call LegendrePolynomials(Ln(i,j),Lpn(i,j),q(i,j),qp(i,j),GLLpoints(j),i-1)
		end do
	end do
	
	! let's calculate legendre transform
	do i=1,(N+1)
		LegenData(i)=0.
		do j=1,N
			LegenData(i)=LegenData(i)+w(j)*Ln(i,j)*OData(j)
		end do
		LegenData(i)=(2.*(i-1)+1.)*LegenData(i)/2.
	end do
	

end subroutine LegendreTransform

subroutine LegendreBackTransform(Odata,LegenData,GLLpoints,N)
	! This function performs Forward Legendre Transform
	! the size of the data 
	integer ,intent(in):: N
	! the original Data
	real*8 ,dimension(N),intent(inout)::OData
	! the transformed Data
	real*8 ,dimension(N+1),intent(inout)::LegenData
	! the grid points and corresponding weights
	real*8 ,dimension(N),intent(in)::GLLpoints
	integer i,j
    !! Legendre Polynomial associated arrays
	real*8 ,allocatable :: Ln(:,:),Lpn(:,:),q(:,:),qp(:,:)
	
	allocate(Ln(N+1,N))
	allocate(Lpn(N+1,N))
	allocate(q(N+1,N))
	allocate(qp(N+1,N))
	
	! let's generate all required Legendre polynomials for Legendre Transform
	do i=1,(N+1)
		do j=1,N
			call LegendrePolynomials(Ln(i,j),Lpn(i,j),q(i,j),qp(i,j),GLLpoints(j),i-1)
		end do
	end do
	
	! let's calculate legendre transform
	do i=1,N
		OData(i)=0.
		do j=1,N+1
			OData(i)=OData(i)+Ln(j,i)*LegenData(j)
		end do
	end do
	
end subroutine LegendreBackTransform