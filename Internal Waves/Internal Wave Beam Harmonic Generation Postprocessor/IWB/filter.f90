subroutine BandPassFilter(Fdata,OData,FourData,wo,bw,N)

! the size of the data 
integer N
! the original data and Filtered Data
real*8 ,dimension(N),intent(inout)::Fdata,OData
 !the fourier transromed data 
complex ,dimension(N),intent(inout):: FourData
! the frequency where the filter centered about and band width 
real*8 ,intent(in)::wo ,bw
! the applicable versio of the frequency and bandwith ,pi number
real*8 w1,bw1,pi
!the filer vector 
complex ,allocatable:: f(:)
! the filtered data in fourier domain
complex,allocatable:: FourFData(:)

complex ,allocatable :: A(:,:),WORK(:,:)

INTEGER IFAX1(13)
complex ,allocatable :: TRIGS1(:)

pi=4.*datan(1.d0)

w1=N*0.5*wo/pi

bw1=bw*N*0.5/pi

!print*,w1,N,wo,bw1

allocate(A(N+2,1))
allocate(WORK(N,1))
allocate(TRIGS1(2*N))

allocate(f(N))
allocate(FourFData(N))

 CALL CFTFAX(N, IFAX1, TRIGS1) 


! let's store the data into FFT matrix
do i=1,N

A(i,1)=OData(i)

end do
 
! here we obtaine the fourier transform of the data
 CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, N+2, N, 1, -1)
 
! here we store the fourier transform of the data
do i=1,N

FourData(i)=A(i,1)

end do


! let's generate the bandpass filter in that do loop 
do i=1,N

!f(i)=1./N
if ( i > w1-bw1 .and. i< w1+bw1 ) then 
f(i)=1./N
end if

if ( i > N-w1-bw1 .and. i< N-w1+bw1 ) then 
f(i)=1./N
end if


end do 


! in that do loop we filter the data in fourier domain
 call vectvect(FourFData,FourData,f,N)

! now let's inverse fourier transform the filtered data
 
! let's store the data into FFT matrix
do i=1,N

A(i,1)=FourFData(i)

end do

 CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, N+2, N, 1, 1)

do i=1,N

FData(i)=Real(A(i,1))

end do


end subroutine BandPassFilter


subroutine vectvect(uf,u,f,N)


integer , intent(in)::N

complex,dimension(N),intent(inout)::u,f,uf

integer i




do i=1,N 

uf(i)=u(i)*f(i)

end do

return 
end subroutine vectvect

            

