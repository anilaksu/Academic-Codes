!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                         !
!                     Written by Anil A. Aksu on 04.07.2020                               !
!   Analytical and Numerical subroutines of Thermoelastic beam solver
!                                                                                         !
!                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get1DHeatWeakForm(DwD, D, w, N)
	! this function generates the weak form of a heat diffusion operator
	integer i, j                                                    ! dummy integer
	integer, intent(in):: N    	                                    ! the number of points on grid 
	real*8, intent(in),dimension(N,N) :: D                          ! Differentiaton Matrix
    real*8, intent(in),dimension(N) :: w                            ! weights
	real*8, intent(inout),dimension(N,N) :: DwD                     ! weak heat diffusion operator on LGL grid points
	real*8, dimension(N,N) :: wD                                   	! Product of the diagonal weight matrix multipled with the differentian matrix
	real*8, dimension(N,N) :: D_trans                               ! Transpose of the differentian matrix
	
	! let's get  weak heat diffusion operator 
	do i=1,N
		wD(i,:) = w(i)*D(i,:)
	end do
	
	call getTranspose(D,D_trans,N,N)
	
	Call DGEMM("N","N",N,N,N,1.d0,D_trans,N,wD,N,0.d0,DwD,N)
	
end subroutine get1DHeatWeakForm

subroutine getHeatSystemXdirection(Kx_inv, Kx_right, DwDx, DiffX,LeftBC, RightBC,  rho, cp, kappa, Lx, dt, w_x, Nx)
	! this function generates the weak form of a heat diffusion operator in x direction
	integer i, j                                                    ! dummy integers
	integer, intent(in):: Nx        	                            ! the number of points on grid in x direction
	logical, intent(in):: LeftBC, RightBC                           ! Left and right boundary conditions if BC == True, it is Neumann condition, else Dirichlet condition
	real*8, intent(in),dimension(Nx) :: w_x                     	! weights in x direction
	real*8, intent(in) :: rho                       		        ! Density
	real*8, intent(in) :: cp	                     		        ! specific heat
	real*8, intent(in) :: Lx                     		        	! Length
	real*8, intent(in) :: kappa                      		        ! Elasticity
	real*8, intent(in) :: dt                      		            ! time step size
	real*8, intent(in),dimension(Nx,Nx) :: DwDx, DiffX              ! weak beam matrix and differentiaon matrices
	real*8, intent(inout),dimension(Nx,Nx) :: Kx_inv, Kx_right      ! Stiffness matrix of the beam in state space form

	! First, we generate the system matrix for the implicit step
	Kx_right = 0.
	! Left boundary condition
	if (LeftBC) then
		Kx_right(1,:) = -DiffX(1,:)
	else
		Kx_right(1,1) = 1.
	end if
	
	! Right boundary condition
	if (RightBC) then
		Kx_right(Nx,:) = DiffX(Nx,:)
	else
		Kx_right(Nx,Nx) = 1.
	end if
	
	do i = 2, Nx-1
		Kx_right(i,:) = 0.5*dt*(kappa/(rho*cp))*DwDx(i,:)                     ! Here we store heat operator normalized with material properties
		Kx_right(i,i) = Kx_right(i,i) +   w_x(i)  
	end do 
	
	call invertSVD(Nx,Kx_right,Kx_inv)
	
	! Second, we generate the system matrix for the explicit
	Kx_right = 0.
	
	do i = 2, Nx-1
		Kx_right(i,:) = -0.5*dt*(kappa/(rho*cp))*DwDx(i,:)                     ! Here we store heat operator normalized with material properties
		Kx_right(i,i) = Kx_right(i,i) +    w_x(i)  
	end do 
	
end subroutine getHeatSystemXdirection

subroutine getCompHeatSystemYdirection(Ky_inv,Ky_right,DwDy,DiffY,TopBC,BotBC,rho1,rho2,cp1,cp2,kappa1,kappa2,g,Ly1,Ly2,dt,w_y,Ny)
	! this function generates the weak form of a heat diffusion operator in y direction for composite beam
	integer i, j                                                    ! dummy integers
	integer, intent(in):: Ny        	                            ! the number of points on grid in y direction
	logical, intent(in):: TopBC, BotBC                              ! top and bottom boundary conditions if BC == True, it is Neumann condition, else Dirichlet condition
	real*8, intent(in),dimension(Nx) :: w_y                     	! weights in y direction
	real*8, intent(in) :: rho1, rho2                       		    ! Density
	real*8, intent(in) :: cp1, cp2	                     		    ! specific heat
	real*8, intent(in) :: Ly1, Ly2                     		        ! Heights in each domain
	real*8, intent(in) :: kappa1, kappa2, g                         ! Conductivity in each domain and interface conductance
	real*8, intent(in) :: dt                      		            ! time step size
	real*8, intent(in),dimension(Ny,Ny) :: DwDy, DiffY              ! weak beam matrix and differentiaon matrices
	real*8, intent(inout),dimension(2*Ny,2*Ny) :: Ky_inv, Ky_right      ! Stiffness matrix of the beam in state space form

	! First, we generate the system matrix for the implicit step
	
	! Bottom boundary condition
	if (BotBC) then
		Ky_right(1,1:Ny) = -DiffY(1,:)
	else
		Ky_right(1,1) = 1.
	end if
	
	! Top boundary condition
	if (TopBC) then
		Ky_right(2*Ny,Ny+1:2*Ny) = DiffY(Ny,:)
	else
		Ky_right(2*Ny,2*Ny) = 1.
	end if
	
	! First domain
	do i = 2, Ny-1
		Ky_right(i,1:Ny) = 0.5*dt*(kappa1/(rho1*cp1))*DwDy(i,:)                     ! Here we store heat operator normalized with material properties in the first domain
		Ky_right(i,i) = Ky_right(i,i) +  w_y(i)  
	end do 
	
	! Second domain
	do i = Ny+2, 2*Ny-1
		Ky_right(i,Ny+1:2*Ny) = 0.5*dt*(kappa2/(rho2*cp2))*DwDy(i-Ny,:)                ! Here we store heat operator normalized with material properties in the second domain
		Ky_right(i,i) = Ky_right(i,i) +   w_y(i-Ny)  
	end do 
	
	! Interface condition
	Ky_right(Ny,1:Ny) = -kappa1*DiffY(Ny,:)              ! Derivative with respect to y
	Ky_right(Ny,Ny+1) = -g                                ! Conductance coefficient for neighbour point
	Ky_right(Ny,Ny) = Ky_right(Ny,Ny) + g                ! Conductance coefficient for the current point
	
	Ky_right(Ny+1,Ny+1:2*Ny) = kappa2*DiffY(1,:)        ! Derivative with respect to y
	Ky_right(Ny+1,Ny) = -g                                ! Conductance coefficient for neighbour point
	Ky_right(Ny+1,Ny+1) = Ky_right(Ny+1,Ny+1) + g        ! Conductance coefficient for current point
	
	call invertSVD(2*Ny,Ky_right,Ky_inv)
	!Ky_inv = Ky_right
	
	! Second, we generate the system matrix for the explicit
	Ky_right = 0.
	
	! First domain
	do i = 2, Ny-1
		Ky_right(i,:) = -0.5*dt*(kappa1/(rho1*cp1))*DwDy(i,:)                     ! Here we store heat operator normalized with material properties
		Ky_right(i,i) = Ky_right(i,i) +    w_y(i)  
	end do 
	! Second domain
	do i = Ny+2, 2*Ny-1
		Ky_right(i,Ny+1:2*Ny) = -0.5*dt*(kappa2/(rho2*cp2))*DwDy(i-Ny,:)                     ! Here we store heat operator normalized with material properties
		Ky_right(i,i) = Ky_right(i,i) +    w_y(i-Ny)  
	end do 
	
end subroutine getCompHeatSystemYdirection

subroutine getCompHeatSystem(K_inv,K_right,DwDx,DiffX,DwDy,DiffY,rho1,rho2,cp1,cp2,kappa1,kappa2,g,dt,Lx,Ly1,Ly2,w_x,w_y,Nx,Ny)
	! this function generates the weak form of a heat diffusion operator in x and y directions 
	integer i, j, k, jstart, jend                                   ! dummy integers
	integer, intent(in):: Nx, Ny        	                        ! the number of points on grid in x and y directions
	real*8, intent(in) :: rho1, rho2, cp1, cp2, kappa1, kappa2! Thermal properties in both domains
	real*8, intent(in) :: g, dt 									! Conductance across domains andd time step
	real*8, intent(in) :: Lx, Ly1, Ly2								! Domain length in x and y directions
	real*8, intent(in),dimension(Nx) :: w_x		                    ! weights in x and y directions
	real*8, intent(in),dimension(Ny) :: w_y      	                ! weights in x and y directions
	real*8, intent(in),dimension(Nx,Nx) :: DwDx, DiffX		        ! Differentiatin matrix and heat matrix in weak form in x direction
	real*8, intent(inout),dimension(Ny,Ny) :: DwDy, DiffY       		! Differentiatin matrix and heat matrix in weak form in y direction
	real*8, intent(inout),dimension(2*Nx*Ny,2*Nx*Ny) :: K_inv, K_right    ! Full heat matrices
	real*8, dimension(Nx*Ny,Nx*Ny) :: K1_inv, K2_inv,K1_right, K2_right   ! Left Handside and right handside  for both domains
	real*8, dimension(Nx*Ny,Nx*Ny) :: K1, K2		  ! Left Handside matrices in x for both domains
	logical	LeftBC, RightBC, BotBC, TopBC             ! boundary conditions if BC == True, it is Neumann condition, else Dirichlet condition

	LeftBC  = .TRUE.                   ! It states Neumann boundary condition on the left boundary
	RightBC = .FALSE.				   ! It states Dirichlet boundary condition on the right boundary
	
	TopBC  = .TRUE.                    ! It states Neumann boundary condition on the top boundary
	BotBC  = .TRUE.	                   ! It states Neumann boundary condition on the bottom boundary
	
	call  get2DHeatSystem(K1_inv,K1_right,DwDx, DiffX, DwDy,DiffY,LeftBC,RightBC, BotBC,TopBC,rho1, cp1, kappa1,dt,w_x,w_y,Nx,Ny)
	! Let's modify the differentiation matrix and the weak form of the heat system in y direction 
	DiffY = DiffY * (Ly1 / Ly2)           
	DwDy  = DwDy  * ((Ly1 / Ly2)**2.)
	call  get2DHeatSystem(K2_inv,K2_right,DwDx, DiffX, DwDy,DiffY,LeftBC,RightBC, BotBC,TopBC,rho2, cp2, kappa2,dt,w_x,w_y,Nx,Ny)
	! Let's recorrect them for upcoming operations
	DiffY = DiffY * (Ly2 / Ly1)
	DwDy  = DwDy  * ((Ly2 / Ly1)**2.)
	
	call invertSVD(Nx * Ny, K1_inv, K1)
	call invertSVD(Nx * Ny, K2_inv, K2)
	
	! Left handside matrix
	K_right = 0.                 ! Note that it will be inverted to form the left handside matrix
	! First part of the domain
	K_right(1:Nx*Ny,1:Nx*Ny) = K1                	  
	! Second part of the domain
	K_right(Nx*Ny+1:2*Nx*Ny,Nx*Ny+1:2*Nx*Ny) = K2               
	
	! Interface conditions between these domains
	do i = 2, Nx
		jstart = (Ny-1) * Nx
		jend = Ny * Nx
		do j = 1, Ny
			! Interface condition in the first domain derivative part
			K_right(jstart+i,i+(j-1)*Nx) = -kappa1*DiffY(Ny,j)              	 ! Derivative with respect to y
			! Interface condition in the second domain derivative part
			K_right(jend+i,jend+i+(j-1)*Nx) = kappa2*DiffY(1,j) * (Ly1 / Ly2)    ! Derivative with respect to y
		end do 
		! Interface condition in the first domain temperature difference part
		K_right(jstart+i,jend+i) = -g                                			! Conductance coefficient for neighbour point
		K_right(jstart+i,jstart+i) = K_right(jstart+i,jstart+i) + g             ! Conductance coefficient for the current point
		! Interface condition in the second domain temperature difference part
		K_right(jend+i,jstart+i) = -g                               		    ! Conductance coefficient for neighbour point
		K_right(jend+i,jend+i) = K_right(jend+i,jend+i) + g    					! Conductance coefficient for current point
	end do
	
	call invertSVD(2 * Nx * Ny, K_right, K_inv)
	! Right handside matrix
	K_right = 0.
	! First part of the domain           
	K_right(1:Nx*Ny,1:Nx*Ny) = K1_right			  
	! Second part of the domain            
	K_right(Nx*Ny+1:2*Nx*Ny,Nx*Ny+1:2*Nx*Ny) = K2_right	
	
end subroutine getCompHeatSystem

subroutine get2DHeatSystem(K_inv,K_right,DwDx, DiffX, DwDy,DiffY,LeftBC,RightBC, BotBC,TopBC,rho, cp, kappa,dt,w_x,w_y,Nx,Ny)
	! this function generates the weak form of a heat diffusion operator in x and y directions 
	integer i, j, k, jstart, jend                                   ! dummy integers
	integer, intent(in):: Nx, Ny        	                        ! the number of points on grid in x and y directions
    logical, intent(in):: LeftBC, RightBC, BotBC, TopBC             ! boundary conditions if BC == True, it is Neumann condition, else Dirichlet condition
	real*8, intent(in),dimension(Nx) :: w_x		                    ! weights in x direction
	real*8, intent(in),dimension(Nx) :: w_y     	                ! weights in y direction
	real*8, intent(in),dimension(Nx,Nx) :: DwDx, DiffX              ! weak heat matrix and differentation matrices in x direction
	real*8, intent(in),dimension(Ny,Ny) :: DwDy, DiffY              ! weak heat matrix and differentation matrices in y direction
	real*8, intent(in) :: dt             	          		        ! Time step size
	real*8, intent(in) :: rho                       		        ! Density
	real*8, intent(in) :: cp	                     		        ! specific heat
	real*8, intent(in) :: kappa                      		        ! Diffusivity
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny) :: K_inv, K_right  ! Full heat matrices
	
	! First, we generate the system matrix for the implicit step
	K_right = 0.
	
	do j = 1, Ny
		jstart = (j-1) * Nx+1     ! Start index
		jend = j * Nx			  ! End index
		
		! Left boundary condition
		if (LeftBC) then
			K_right(jstart,jstart:jend) = -DiffX(1,:)
		else
			K_right(jstart,jstart) = 1.
		end if

		! Right boundary condition
		if (RightBC) then
			K_right(jend, jstart:jend) = DiffX(Nx,:)
		else
			K_right(jend ,jend) = 1.
		end if

		do i = 2, Nx-1
			K_right(jstart+i-1,jstart:jend) = 0.5*dt*w_y(j) * (kappa/(rho*cp))*DwDx(i,:)                     ! Here we store heat operator normalized with material properties
			K_right(jstart+i-1,jstart+i-1) = K_right(jstart+i-1,jstart+i-1) +  w_y(j)*w_x(i)                 ! Addition of mass matrix
		end do 
	end do 
	
	! Then, we add heat matrices in y direction to the left handside system matrix 
	
	do i =1, Nx
		do j=1,Ny 
			if( i == 1 .OR.  i == Nx ) then
				continue
			else if ( j == 1 ) then               ! Boundary condition 
					K_right(i+(j-1)*Nx,:) = 0.
					! Bottom boundary condition
					do k=1,Ny
						if (BotBC) then
							K_right(i+(j-1)*Nx,i+(k-1)*Nx) = -DiffY(1,k)
						else
							if(k == 1) then
								K_right(i+(j-1)*Nx,i+(k-1)*Nx) = 1.
							else
								K_right(i+(j-1)*Nx,i+(k-1)*Nx) = 0.
							end if
						end if
					end do
				 else if ( j == Ny ) then               ! Boundary condition 
					K_right(i+(j-1)*Nx,:) = 0.
					! Top boundary condition
					do k=1,Ny
						if (TopBC) then
							K_right(i+(j-1)*Nx,i+(k-1)*Nx) = DiffY(Ny,k)
						else
							if(k == Ny) then
								K_right(i+(j-1)*Nx,i+(k-1)*Nx) = 1.
							else
								K_right(i+(j-1)*Nx,i+(k-1)*Nx) = 0.
							end if
						end if	
					end do
			else	! In domain
				do k=1,Ny
					K_right(i+(j-1)*Nx,i+(k-1)*Nx) = K_right(i+(j-1)*Nx,i+(k-1)*Nx)+0.5*dt*w_x(i)*(kappa/(rho*cp))*DwDy(j,k)
				end do
			end if
		end do
	end do
	
	call invertSVD(Nx * Ny, K_right, K_inv)
	!K_inv = K_right
	! Second, we generate the system matrix for the explicit
	K_right = 0.
	do j = 1, Ny
		jstart = (j-1) * Nx+1     ! Start index
		jend = j * Nx			  ! End index
		do i = 2, Nx-1
			K_right(jstart+i-1,jstart:jend) = -0.5*dt*w_y(j) * (kappa/(rho*cp))*DwDx(i,:)                     ! Here we store heat operator normalized with material properties
			K_right(jstart+i-1,jstart+i-1) = K_right(jstart+i-1,jstart+i-1) +  w_y(j)*w_x(i)                  ! Addition of mass matrix
		end do 
	end do 
	
	! Then, we add heat matrices in y direction to the right handside system matrix 	
	do i =1, Nx
		do j=1,Ny 
			do k=1,Ny
				if( i == 1 .OR.  i == Nx ) then
					continue
				else if ( j == 1 ) then               ! Boundary condition 
						K_right(i+(j-1)*Nx,:) = 0.    ! Bottom boundary condition
				else if ( j == Ny ) then         	  ! Boundary condition 
						K_right(i+(j-1)*Nx,:) = 0.    ! Top boundary condition
				else	! In domain
						K_right(i+(j-1)*Nx,i+(k-1)*Nx) = K_right(i+(j-1)*Nx,i+(k-1)*Nx)-0.5*dt*w_x(i)*(kappa/(rho*cp))*DwDy(j,k)
				end if
			end do
		end do
	end do
	
end subroutine get2DHeatSystem


subroutine getTotalCompBoundarForcing(F_boundary,LeftBC, RightBC, TopBC, BotBC, uL,uR,uT,uB,kappa1,kappa2,Nx,Ny)
	! this function generates the forcing due to boundary conditions
	integer i, j                                                    ! dummy integers
	integer, intent(in):: Nx, Ny        	                        ! the number of points on grid in x and y direction
	logical, intent(in):: LeftBC, RightBC, TopBC, BotBC             ! All boundary conditions if BC == True, it is Neumann condition, else Dirichlet condition
	real*8, intent(in),dimension(2*Ny) :: uL, uR                   	! boundary forcing in x direction
	real*8, intent(in),dimension(Nx) :: uT, uB                   	! boundary forcing in x direction
	real*8, intent(in) :: kappa1, kappa2                      		! conductivities in domain 1 and domain 2
	real*8, intent(inout),dimension(2*Nx*Ny) :: F_boundary	        ! Forcing vector at boundaries

	! Boundary forcing in x direction
	do i = 1, 2*Ny
		if ( i <= Ny) then        ! First domain
			! Left boundary condition
			if (LeftBC) then
				F_boundary(Nx*(i-1)+1) = uL(i)/kappa1
			else
				F_boundary(Nx*(i-1)+1) = uL(i)
			end if

			! Right boundary condition
			if (RightBC) then
				F_boundary(Nx*i) = uR(i)/kappa1
			else
				F_boundary(Nx*i) = uR(i)
			end if
		else
			! Left boundary condition
			if (LeftBC) then
				F_boundary(Nx*(i-1)+1) = uL(i)/kappa2
			else
				F_boundary(Nx*(i-1)+1) = uL(i)
			end if

			! Right boundary condition
			if (RightBC) then
				F_boundary(Nx*i) = uR(i)/kappa2
			else
				F_boundary(Nx*i) = uR(i)
			end if
		end if 
	end do 
	
	! Boundary forcing in y direction
	do i = 1, Nx-1
		! Bottom boundary condition
		if (BotBC) then
			F_boundary(i) = F_boundary(i)+uB(i)/kappa1
		else
			F_boundary(i)= uB(i)
		end if

		! Top boundary condition
		if (TopBC) then
			F_boundary((2*Ny-1)*Nx+i) = F_boundary((2*Ny-1)*Nx+i)+uT(i)/kappa2
		else
			F_boundary((2*Ny-1)*Nx+i) = uT(i)
		end if
	end do 
	
	!print*, "Conductivity in second domain: ", kappa2 

end subroutine getTotalCompBoundarForcing

subroutine getTotalBoundarForcing(F_boundary,LeftBC, RightBC, TopBC, BotBC, uL,uR,uT,uB,kappa,Nx,Ny)
	! this function generates the forcing due to boundary conditions in single material 
	integer i, j                                                    ! dummy integers
	integer, intent(in):: Nx, Ny        	                        ! the number of points on grid in x and y direction
	logical, intent(in):: LeftBC, RightBC, TopBC, BotBC             ! All boundary conditions if BC == True, it is Neumann condition, else Dirichlet condition
	real*8, intent(in),dimension(Ny) :: uL, uR                   	! boundary forcing in x direction
	real*8, intent(in),dimension(Nx) :: uT, uB                   	! boundary forcing in x direction
	real*8, intent(in) :: kappa			                     		! conductivity
	real*8, intent(inout),dimension(Nx*Ny) :: F_boundary	        ! Forcing vector at boundaries

	 F_boundary = 0.                         ! Let's reset the boundary forcing vector
	! Boundary forcing in x direction
	do i = 1, Ny
		! Left boundary condition
		if (LeftBC) then
			F_boundary(Nx*(i-1)+1) = uL(i)/kappa
		else
			F_boundary(Nx*(i-1)+1) = uL(i)
		end if

		! Right boundary condition
		if (RightBC) then
			F_boundary(Nx*i) = uR(i)/kappa
		else
			F_boundary(Nx*i) = uR(i)
		end if
	end do 
	
	! Boundary forcing in y direction
	do i = 2, Nx-1
	! Bottom boundary condition
		if (BotBC) then
			F_boundary(i) = F_boundary(i)+uB(i)/kappa
		else
			F_boundary(i)= uB(i)
		end if

		! Top boundary condition
		if (TopBC) then
			F_boundary((2*Ny-1)*Nx+i) = F_boundary((2*Ny-1)*Nx+i)+uT(i)/kappa
		else
			F_boundary((2*Ny-1)*Nx+i) = uT(i)
		end if
	end do 

end subroutine getTotalBoundarForcing

subroutine getThermalMoment(M_t, Temp, alpha, E, Ly, w_y, Nx, Ny)
	! this function computes thermal moment at each time step
	integer i, j                                                    ! dummy integer
	integer, intent(in):: Nx, Ny    	                            ! the number of points on grid in x and y direction
	real*8, intent(in) :: Ly                      		            ! Height
	real*8, intent(in) :: E                      		            ! Elasticity
	real*8, intent(in) :: alpha                    		            ! Thermal expansion coefficient
	real*8, intent(in),dimension(Ny) :: w_y                      	! weights in y direction
	real*8, intent(in),dimension(Nx*Ny) ::Temp		                ! Temperature distribution
	real*8, intent(inout),dimension(Ny) :: M_t                      ! Thermal Moment
	! boundary forcing in x direction
	do i=1,Nx
		M_t(i) = 0. 
		do j = 1,Ny
			M_t(i) = M_t(i) + 0.5*Ly*alpha*E*w_y(j)*Temp(i+Nx*(j-1))
		end do 
	end do

	
end subroutine  getThermalMoment