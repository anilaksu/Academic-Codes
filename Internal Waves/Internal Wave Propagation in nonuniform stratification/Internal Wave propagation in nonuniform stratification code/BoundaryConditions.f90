!! this routine generates all the associated boundary condition matrices
subroutine OneDBC(BCMatrix,LB,RB,N)
	!! This function generates boundary conditions 
	!! if LB==1, it is neumann boundary condtion 
	!! if LB==0, it is dirichlet boundary condition
	integer i,j,k 
	! the number of grid points and left boundary condition LB and right boundary condtion RB
	integer, intent(in):: N, LB, RB
	!the differentiation matrix
	real*8, intent(inout),dimension(N,N):: BCMatrix
	!! Differentiation Matrix 
	real*8,allocatable:: D1(:,:)

	! the allocation of first derivative matrix 
	allocate(D1(N,N))

	! let's first generate differentiation matrix
	call FirstDiff(D1,N)
	BCMatrix=0.

	!! dirichlet-dirichlet boundary condition
	if(LB==0 .and. RB==0) then
		BCMatrix(1,1)=1.
		BCMatrix(N,N)=1.
		 ! dirichlet-neumann boundary condition
	else if (LB==0 .and. RB==1) then
			BCMatrix(1,1)=1.
			BCMatrix(N,:)=D1(N,:)
		  ! neumann-dirichlet boundary condition
		 else if (LB==1 .and. RB==0) then
			 	BCMatrix(1,:)=D1(1,:)
				BCMatrix(N,N)=1.
	 		   ! neumann-neumann boundary condition
			  else if (LB==1 .and. RB==1) then
				 BCMatrix(1,:)=D1(1,:)
				 BCMatrix(N,:)=D1(N,:)
	end if 

	
end subroutine OneDBC

subroutine OneDBCApply(BCMatrix,GEMatrix,SysMatrix,N)
	!! This function generates system matrix
	integer i
	! the number of grid points 
	integer, intent(in):: N
	!the boundary condition matrix, governing equation matrix and system matrix
	real*8, intent(inout),dimension(N,N):: BCMatrix,GEMatrix,SysMatrix
	!real*8, intent(inout),dimension(N-2,N-2):: SysMatrix
	
	do i=1,N
		if(i==1 .or. i==N) then
			!SysMatrix(i,:)=BCMatrix(i,:)
			SysMatrix(i,:)=BCMatrix(i,:)
		else
			SysMatrix(i,:)=-1.*GEMatrix(i,:)
		end if
	end do
	
	!do i=1,N-2
	!	SysMatrix(i,:)=GEMatrix(i+1,2:(N-1))
	!end do
	
end subroutine OneDBCApply

subroutine OneDSEMBC(SEMBCMatrix,Nsub,Ngrid)
	!! this function calculates the discontinous laplace operator in 1D with Nsub subdomain with 
	!! Ngrid grid points in each domain
	integer i,j
	!! the start and end indices
	integer jstart, jend
	! the number of grid points
	integer, intent(in):: Nsub,Ngrid
	!the boundary condition matrix
	real*8, intent(inout),dimension(Nsub*Ngrid,Nsub*Ngrid):: SEMBCMatrix	
	!! Differentiation Matrix 
	real*8,allocatable:: D1(:,:)
	! the allocation of first derivative matrix 
	allocate(D1(Ngrid,Ngrid))
	! let's first generate differentiation matrix
	call FirstDiff(D1,Ngrid)
	BCMatrix=0.
	
	! the start and end points of the right boundary 
	jstart=Ngrid*(Nsub-1)+1
	jend=Ngrid*Nsub
	
	!! dirichlet-dirichlet boundary condition
	if(LB==0 .and. RB==0) then
		SEMBCMatrix(1,1)=1.
		SEMBCMatrix(jend,jend)=1.
		 ! dirichlet-neumann boundary condition
	else if (LB==0 .and. RB==1) then
			SEMBCMatrix(1,1)=1.
			SEMBCMatrix(jend,jstart:jend)=D1(Ngrid,:)
		  ! neumann-dirichlet boundary condition
		 else if (LB==1 .and. RB==0) then
			 	SEMBCMatrix(1,1:Ngrid)=D1(1,:)
				SEMBCMatrix(jend,jend)=1.
	 		   ! neumann-neumann boundary condition
			  else if (LB==1 .and. RB==1) then
				 SEMBCMatrix(1,1:Ngrid)=D1(1,:)
				 SEMBCMatrix(jend,jstart:jend)=D1(Ngrid,:)
	end if 

	
end subroutine OneDSEMBC

subroutine OneDSEMBCApply(BCMatrix,GEMatrix,SysMatrix,Nsub,Ngrid)
	!! This function generates system matrix
	integer i
	!! the start and end indices
	integer jend
	! the number of subintervals and the number of grid points 
	integer, intent(in):: Nsub,Ngrid
	!the boundary condition matrix, governing equation matrix and system matrix
	real*8, intent(inout),dimension(Nsub*Ngrid,Nsub*Ngrid):: BCMatrix,GEMatrix,SysMatrix
	
	! the start and end points of the right boundary 
	jend=Ngrid*Nsub
	
	do i=1,jend
		if(i==1 .or. i==jend) then
			SysMatrix(i,:)=BCMatrix(i,:)
		else
			SysMatrix(i,:)=-1.*GEMatrix(i,:)
		end if
	end do
	
	!do i=1,N-2
	!	SysMatrix(i,:)=GEMatrix(i+1,2:(N-1))
	!end do
	
end subroutine OneDSEMBCApply

subroutine TwoDBCFirst(BCMatrix,LB,RB,BB,TB,Nx,Ny)
	!! This function generates boundary conditions 
	!! if LB==1, there is no condition
	!! if LB==0, it is dirichlet boundary condition
	integer i,j,k 
	!! the start and end indices
	integer jstart, jend
	! the number of grid points in x and y direction 
	! and left boundary condition LB and right boundary condtion RB
	! the top boundary condition and the bottom boundary condition
	integer, intent(in):: Nx, Ny, LB, RB, TB, BB
	!the differentiation matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: BCMatrix
	
	BCMatrix=0.
	
	! let's impose the boundary conditions in x coordinate
	do i=1,Ny
		jstart=(i-1)*Nx+1
		jend=i*Nx
		!! dirichlet-dirichlet boundary condition
		if(LB==0 .and. RB==0) then
			BCMatrix(jstart,jstart)=1.
			BCMatrix(jend,jend)=1.
			 ! dirichlet on left boundary
		else if (LB==0 .and. RB==1) then
				BCMatrix(jstart,jstart)=1.
			  ! neumann on right boundary condition
			 else if (LB==1 .and. RB==0) then
					BCMatrix(jend,jend)=1.
		end if 
	end do
	
	! let's impose the boundary conditions in y coordinate
	
		jstart=1
		jend=Ny
		!! dirichlet-dirichlet boundary condition
		if(TB==0 .and. BB==0) then
			do j=1,Ny
				do i=2,Nx
					! the bottom boundary 
					if(j==1) then 
						BCMatrix(jstart+i-1,jstart+i-1)=1.+BCMatrix(jstart+i-1,jstart+i-1)
						! the top boundary
					else if(j==Ny) then 					
						BCMatrix((jend-1)*Nx+i,(jend-1)*Nx+i)=1.+BCMatrix((jend-1)*Nx+i,(jend-1)*Nx+i)
					end if 
				end do
			end do
			 ! dirichlet on bottom boundary condition
		else if (TB==1 .and. BB==0) then
				do j=1,Ny
					do i=2,Nx
						! the bottom boundary 
						if(j==1) then 
							BCMatrix(jstart+i-1,jstart+i-1)=1.+BCMatrix(jstart+i-1,jstart+i-1)
						end if
							
					end do
				end do
			  ! dirichlet on top boundary condition
			 else if (TB==0 .and. BB==1) then
					 	do j=1,Ny
							do i=2,Nx
								! the top boundary
								if(j==Ny) then 					
									BCMatrix((jend-1)*Nx+i,(jend-1)*Nx+i)=1.+BCMatrix((jend-1)*Nx+i,(jend-1)*Nx+i)
								end if 
							end do
						end do
					 
		end if 
	
	
end subroutine TwoDBCFirst

subroutine TwoDBCApplyFirst(BCMatrix,GEMatrix,SysMatrix,Nx,Ny)
	!! This function applies boundary condition for first order partial differential equation 
	integer i,j
	!! the start and end indices
	integer jstart,jend
	! the number of subintervals and the number of grid points 
	integer, intent(in):: Nx,Ny
	!the boundary condition matrix, governing equation matrix and system matrix
	real*8, intent(inout),dimension(Nx*Ny,Nx*Ny):: BCMatrix,GEMatrix,SysMatrix
	
	! the start and end points of the right boundary 
	SysMatrix=0.
	

	! the boundary conditions added
	do i=1,Nx*Ny
		! if there is no boundary condition
		if(sum(BCMatrix(i,:))<10.**-2.) then	
			SysMatrix(i,:)=GEMatrix(i,:)
		else		
			SysMatrix(i,:)=BCMatrix(i,:)
		end if
	end do


	
end subroutine TwoDBCApplyFirst
