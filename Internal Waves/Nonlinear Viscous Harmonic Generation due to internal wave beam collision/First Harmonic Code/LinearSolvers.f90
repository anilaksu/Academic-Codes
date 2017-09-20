!! This routine includes all linear solvers 
!! This routine will be updated as new linear solvers added
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