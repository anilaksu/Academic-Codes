!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                         !
!            Written by Anil A. Aksu on 04.07.2020 update on 16.05.2021    			      !
!  		Analytical solution of Euler beam under cantilevered boundary conditions	      !
!                                                                                         !
!                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getRootFunction(f,dummy)

	real*8   arg_1, arg_2, arg_3, arg_4		                        ! dummy summation argument 
	real*8, intent(in) :: dummy										! input argument
	real*8, intent(inout) :: f                                      ! output function
	

	arg_1 = dcos(dummy)*dsinh(dummy)
		
	arg_2 = dsin(dummy)*dcosh(dummy)

	
	f = arg_1 + arg_2  					 ! Reduced eigenrelation
	!f = (4.*dummy*(arg_1*arg_2 - arg_3*arg_4))
	!/((L**2.)*(dcos(dummy)+dsin(dummy))*(dcos(dummy)-dsin(dummy)))
end subroutine getRootFunction

subroutine getEigenRoot(x_root, x_ls, x_rs)

	integer   i,j				                        ! dummy indices
	real*8, intent(in) :: x_ls, x_rs					! left and right start points 
	real*8, intent(inout) :: x_root 					! root
	
	real*8  x_l, x_r, x_m								! input arguments
	real*8  f_l, f_r, f_m                               ! left, right and middle function
	real*8 Cutoff										! error cut off
	
	x_l = x_ls
	x_r = x_rs 
	Cutoff	= 10.**(-1.)

	f_m = 1.                                            ! to initialize while loop here we set it to random large number
	!do while (dabs(f_m) > Cutoff)
	do i=1,100
		!print*, "Okey I am here"
		call getRootFunction(f_l,x_l)                   ! left value
		call getRootFunction(f_r,x_r)					! right value
		
		x_m = 0.5*(x_l+x_r)                             ! middle point
		call getRootFunction(f_m,x_m)					! middle value
		!print*, f_m
		if (f_m * f_r < 0.) then 						! let's replace with left point
			x_l = x_m									
		else											! let's replace with right point
			x_r = x_m
		end if

		x_root = x_m
	end do
	

end subroutine getEigenRoot

subroutine getBeamEigenvalues(Nfreq,Evalues,N,L,alpha)
	! this function calculates eigenvalues and corresponding natural frequencies for beam modes
	integer i, j, k                                                 ! dummy integer
	integer, intent(in):: N        	                            	! Number of eigenvalues
	real*8 dummy, prev                                             ! dummy and previous dummy argument
	real*8 f_cur, f_prev, f_test
	real*8, intent(in):: L, alpha                                   ! domain length
	real*8  pi                                                      ! pi constant
	real*8, intent(inout),dimension(N) :: Evalues, Nfreq    		! Eigenvalue and Natural frequency arrays\
	real*8, dimension(3) :: Beta_n 			   						! Eigencoefficients

    pi = 4.*datan(1.d0) 
	Evalues = 0.
	i = 0
	j = 0
	f_prev = 0.           ! Initial start value
	prev = 0.			  ! Initical start value
	! let's get the eigenvalues
	do while (j < N) 
	!do j=1, N
		dummy = i/10000.
		
		call getRootFunction(f_cur ,dummy)

		if (f_cur*f_prev < 0) then
			j = j+1 
			call getEigenRoot(Evalues(j), prev, dummy)
			call getRootFunction(f_test, Evalues(j))
			Evalues(j) = Evalues(j)/(.5*L)
			Nfreq(j) = alpha*(Evalues(j)**2.)
			!print*, j, dummy, Evalues(j), f_cur, f_prev, f_test
			!, j
		end if 
		i = i + 1	      ! let's update i value
		! here we set these to previuos values to check the sign change
		prev = dummy 
		f_prev = f_cur
	end do
	
	do i = 1, N
		call getRootFunction(f_test, 2.*Evalues(i)/L)
		!print*, "Eigenvalue "
		!print*, i, Evalues(i),f_test
	end do
	

	
end subroutine getBeamEigenvalues

subroutine getBeamModeCoeffs(Beta_n,Evalues,N,L)
	! this function calculates corresponding beam modes with their normalized basis function coefficients 
	integer i		                                                ! dummy integer
	integer, intent(in):: N        	                            	! Number of eigenvalues
	real*8, intent(in):: L		                                    ! domain length
	real*8, intent(in),dimension(N) :: Evalues    				! Eigenvalue and Natural frequency arrays
	real*8, intent(inout), dimension(N,2) :: Beta_n 			    ! Eigencoefficients

	do i = 1, N
		Beta_n(i,1) = 1.
		Beta_n(i,2) = -dcos(.5*Evalues(i)*L)/dcosh(.5*Evalues(i)*L)
	end do
	
end subroutine getBeamModeCoeffs

subroutine getCoefficientMatrix(M_nm, Evalue_n,Evalue_m, w_x, x_grid, Nx,L)
	! this function calculates coefficient matrix
	integer i, j, k                                             ! dummy integer
	integer, intent(in)::Nx         	                       	! Number of grid points
	real*8, intent(in):: Evalue_n,Evalue_m                      ! Two corresponding eigenvalues
	real*8, intent(in):: L		                                ! domain length
	real*8, intent(in),dimension(Nx) :: w_x,x_grid 			   	! LGL weights and grid points
	real*8, intent(inout),dimension(2,2) :: M_nm			   	! Corresponding coefficient matrix
	real*8  dummy
	M_nm = 0.
	! First element of first row
	dummy = 0.
	do i=1,Nx
		dummy = dummy + w_x(i)*dcos(Evalue_n*x_grid(i))*dcos(Evalue_m*x_grid(i))
	end do

	M_nm(1,1) = 0.5*dummy*L
	!	print*, dummy, L, M_nm(1,1)
	! Second element of first row
	dummy = 0.
	do i=1,Nx
		dummy = dummy + w_x(i)*dcos(Evalue_n*x_grid(i))*dcosh(Evalue_m*x_grid(i))
	end do
	
	M_nm(1,2) = 0.5*dummy*L
	
	! First element of second row
	dummy = 0.
	do i=1,Nx
		dummy = dummy + w_x(i)*dcosh(Evalue_n*x_grid(i))*dcos(Evalue_m*x_grid(i))
	end do
	
	M_nm(2,1) = 0.5*dummy*L
	! Second element of second row
	dummy = 0.
	do i=1,Nx
		dummy = dummy + w_x(i)*dcosh(Evalue_n*x_grid(i))*dcosh(Evalue_m*x_grid(i))
	end do
	
	M_nm(2,2) = 0.5*dummy*L


end subroutine getCoefficientMatrix

subroutine getCoefficientVector(V, Beta_m, Evalue_m, w_x, q, x_grid, Nx,L)
	! this function calculates coefficient matrix
	integer i, j, k                                             ! dummy integer
	integer, intent(in)::Nx         	                       	! Number of grid points
	real*8, intent(in):: Evalue_m  			                    ! Corresponding eigenvalue
	real*8, intent(in):: q  				                    ! Uniform distirbuted load
	real*8, intent(in):: L		                                ! domain length
	real*8, intent(in),dimension(Nx) :: w_x,x_grid 			   	! LGL weights and grid points
	real*8, intent(in),dimension(2) :: Beta_m			       	! Corresponding coefficient matrix
	real*8, intent(inout) :: V	            		        	! Corresponding coefficient matrix
	
	real*8  dummy
	
	V = 0.
	! First element 
	dummy = 0.
	do i=1,Nx
		dummy = dummy + q*w_x(i)*dcos(Evalue_m*x_grid(i))
	end do

	V = V + Beta_m(1)* 0.5*dummy*L

	! Second element 
	dummy = 0.
	do i=1,Nx
		dummy = dummy + q*w_x(i)*dcosh(Evalue_m*x_grid(i))
	end do
	
	V = V + Beta_m(2)* 0.5*dummy*L
	
	
end subroutine getCoefficientVector

subroutine getBeamModes(M,Evalues,N, Nx, w_x, x_grid, L)
	! this function calculates corresponding beam modes with their normalized basis function coefficients 
	integer i, j, k                                                 ! dummy integer
	integer, intent(in):: N, Nx        	                            ! Number of eigenvalues and grid points
	real*8, intent(in),dimension(N) :: x_grid				     	! Grid points
	real*8, intent(in) :: L 								     	! Domain length
	real*8, intent(in),dimension(N) :: Evalues				     	! Eigenvalues
	real*8, intent(inout),dimension(N,Nx) :: M                      ! System level forcing eigenvalue relation
	real*8,dimension(N,2) :: Beta_n									! Eigenrelations between modes for each eigenvalue

	call getBeamModeCoeffs(Beta_n,Evalues,N,L)					    ! Here we calculate eigenrelation between basis function for each eigenvalue
	
	!call getCoefficientMatrix(M_nm, Evalue_n,Evalue_m, w_x, x_grid, Nx,L)
	do i=1,N
		do j=1,Nx
			M(i,j) = Beta_n(i,1)*dcos(Evalues(i)*x_grid(j))+Beta_n(i,2)*dcosh(Evalues(i)*x_grid(j)) 
		end do
	end do
	
end subroutine getBeamModes

subroutine getBeamAnalyticalSolution(w_disp, EI, rho, q, t, N, Nx, w, x_grid, L)
	! this function generates the analytocal solution 
	integer i, j, k, m                                                 ! dummy integer
	integer, intent(in):: N, Nx        	                            ! Number of eigenvalues and grid points
	real*8, intent(in),dimension(Nx) :: x_grid				     	! Grid points
	real*8, intent(in),dimension(Nx) :: w					     	! weights
	real*8, intent(inout),dimension(Nx) :: w_disp				     	! displacement
	real*8, intent(in) :: L 								     	! Domain length
	real*8, intent(in) :: EI	 								    ! Beam stiffness
	real*8, intent(in) :: rho	 								    ! Beam density
	real*8, intent(in) :: q     								    ! Uniform distribute load
	real*8, intent(in) :: t		 								    ! time
	real*8, dimension(N) :: Evalues, Nfreq	         				! Eigenvalues
	real*8, dimension(N,Nx) :: Modes                     		    ! System level forcing eigenvalue relation
	real*8, dimension(N,2) :: Beta_n								! Eigenrelations between modes for each eigenvalue
	real*8, dimension(2,2) :: M_nm						     		! Element-wise coefficient matrix
	real*8, dimension(2) :: dummy_vect						     	! Dummy vector
	real*8, dimension(N,N) :: M_sys, M_inv				     		! System matrix and its inverse
	real*8, dimension(N) :: F_sys						     		! System forcing
	real*8, dimension(N) :: d_n						     		    ! Mode coefficients
	real*8 dummy
	
	call getBeamEigenvalues(Nfreq,Evalues,N,L, dsqrt(EI/rho))              ! Here we calculate eigenvalues for the Euler beam equation 
	call getBeamModeCoeffs(Beta_n,Evalues,N,L)					    ! Here we calculate eigenrelation between basis function for each eigenvalue
	call getBeamModes(Modes,Evalues,N, Nx, w, x_grid, L)
	k = 2 ! Just to call DEGMM function
	m = 1 ! Just to call DEGMM function
	
	do i =1,N
		do j=1,N
			call getCoefficientMatrix(M_nm, Evalues(j),Evalues(i), w, x_grid, Nx,L)	   ! Here we generate element-wise coefficient matrix	
			call matvect(M_nm,Beta_n(i,:),dummy_vect,k)
			call vectvect(dummy_vect,Beta_n(j,:),dummy,k)
			
			M_sys(i,j) = dummy*(Evalues(j)**4.)
			
		end do 
			call getCoefficientVector(dummy, Beta_n(i,:), Evalues(i), w, q, x_grid, Nx,L) ! Forcing vector calculation
			F_sys(i) = -dummy/EI	
	end do 
	
	call invertSVD(N,M_sys,M_inv)    ! Let's calculate the inverse of the system matrix
	call DGEMM("N","N",N,m,N,1.d0,M_inv,N,F_sys,N,0.d0,d_n,N)
	
	! Let's generate the full solution
	w_disp = 0.
	do i=1, N
		w_disp = w_disp + d_n(i)*(1.-dcos(Nfreq(i)*t))*Modes(i,:)
		!w_disp = w_disp + d_n(i)*Modes(i,:)
	end do

	!do i=1,  N
	!	print*,i, d_n(i)           ! Here we output corresponding eigencoefficient
	!end do	
	
	do i=1,  Nx
		!print*,i , x_grid(i), w_disp(i),Modes(1,i)
	end do	
end subroutine getBeamAnalyticalSolution