!!! This file contains all grid generation operations including:
!!! 1) 1D GLL and Fourier Points
!!! 2) 2D GLL-GLL, GLL-Fourier, Fourier-Fourier points
!!! 3) Coordinate Mapping functions like Numerical Jacobian Matrix calculations
subroutine OneDTransform(x_trans,x_surf,x_grid,Nsurf,Ngrid)
	integer i
	! the number of points on grid and surface
	integer, intent(in):: Nsurf,Ngrid
	!the grid coordinate array
	real*8, intent(inout),dimension(Ngrid,2):: x_grid,x_trans
	!the grid coordinate array
	real*8, intent(inout),dimension(Nsurf,2):: x_surf
	! the length of domain
	real*8 domLength
	! the grid points array and corresponding weight
	real*8, allocatable:: ds_surf(:)
	! the transformed surface 
	real*8, allocatable:: x_surtrans(:,:)
	! all required matrices for mapping
	real*8, allocatable:: A(:,:),Ainv(:,:),Amap(:,:),Apro(:,:)
	allocate(ds_surf(Nsurf))
	allocate(x_surtrans(Nsurf,2))
	domLength=0.
	ds_surf(1)=0.
	
	! let's calculate the arclength at each point
	do i=2,Nsurf
		ds_surf(i)=sqrt((x_surf(i-1,1)-x_surf(i,1))**2.+(x_surf(i-1,2)-x_surf(i,2))**2.)
		domLength=domLength+ds_surf(i)
	!	print*, ds_surf(i),domLength
	end do
	!do i=1, Nsurf
	!	print*, x_surf(i,:)
	!end do
	
	! Now let's map them to mother interval -1 to 1
	do i=2,Nsurf
		if(x_grid(1,1)==x_grid(2,1)) then
			x_surtrans(1,2)=-1
			x_surtrans(i,2)=x_surtrans(i-1,2)+2.*ds_surf(i)/domLength
			x_surtrans(:,1)=x_grid(1,1)
		else
				x_surtrans(1,1)=-1
				x_surtrans(i,1)=x_surtrans(i-1,1)+2.*ds_surf(i)/domLength
				x_surtrans(:,2)=x_grid(1,2)
			end if
	end do

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                                                     !
	!             Coordinate Transformation               !
	!                                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	allocate(A(Nsurf,Nsurf))
	allocate(Ainv(Nsurf,Nsurf))
	allocate(Amap(Ngrid,Nsurf))
	allocate(Apro(Ngrid,Nsurf))
	
	! let's generate the coefficient matrix
	call matgen(A,x_surtrans(:,1),x_surtrans(:,2),Nsurf)
	call invertSVD(Nsurf,A,Ainv)	
	call matgenxi(Amap,x_surtrans(:,1),x_grid(:,1),x_surtrans(:,2),x_grid(:,2),Ngrid,Nsurf)
	call matmultnon(Amap,Ainv,Apro,Ngrid,Nsurf,Nsurf)
	
	! let's perform mapping
	! for x coordinate
	call matvectnon(Apro,x_surf(:,1),x_trans(:,1),Ngrid,Nsurf)
	! for y coordinate
	call matvectnon(Apro,x_surf(:,2),x_trans(:,2),Ngrid,Nsurf)
	
	
end subroutine OneDTransform

subroutine TwoDMap(x_grid,w,Ngrid)
	! this function generates 2D map with Ngrid x Ngrid Points
	integer i
	! the number of points on grid  
	integer, intent(in):: Ngrid
	!the grid coordinate array
	real*8, intent(inout),dimension(Ngrid*Ngrid,2):: x_grid
	!the weights
	real*8, intent(inout),dimension(Ngrid):: w
	! 1D GLL points and corresponding weight w
	real*8, allocatable::x_1D(:)
	! all required matrices for mapping
	! 1D GLL points and corresponding weight
	allocate(x_1D(Ngrid))
	
	call GLLPoints(x_1D,w,Ngrid)
	
	! let's generate the grid
	do i=1,Ngrid		
		if(i==1)then
			kstart=1
		else
			kstart=(i-1)*Ngrid+1
		end if
		kend=Ngrid*i
		! the x coordinate of the grid
		x_grid(kstart:kend,1)=x_1D
		! the y coordinate of the grid
		x_grid(kstart:kend,2)=-1.*x_1D(i)
	end do
		
end subroutine TwoDMap

subroutine TwoDMapXY(x_grid,wx,wy,Nx,Ny)
	integer i,j
	! the number of points on grid in x nad y direction  
	integer, intent(in):: Nx, Ny
	!the grid coordinate array
	real*8, intent(inout),dimension(Ny*Ny,2):: x_grid
	!the weights in x direction
	real*8, intent(inout),dimension(Nx):: wx
	!the weights in y direction
	real*8, intent(inout),dimension(Ny):: wy
	! 1D GLL points and corresponding weight w
	real*8, allocatable::x_1Dx(:),x_1Dy(:)
	! all required matrices for mapping
	! 1D GLL points and corresponding weight
	allocate(x_1Dx(Nx))
	allocate(x_1Dy(Ny))
	
	call GLLPoints(x_1Dx,wx,Nx)
	call GLLPoints(x_1Dy,wy,Ny)
	
	! let's generate the grid
	do i=1,Nx
		do j=1,Ny
			! the x coordinate of the grid
			x_grid((j-1)*Nx+i,1)=x_1Dx(i)
			! the y coordinate of the grid
			x_grid((j-1)*Nx+i,2)=x_1Dy(j)
		end do
	end do
		
end subroutine TwoDMapXY

subroutine TwoDMeshTransform(x_trans,x_surf,x_grid,Nsurf,Ngrid)
	! this routine maps and generates arbitrary surface to 2D square domain
	! x_surf is the surface elements of arbitrary domain, it is given in counter-clockwise direction
	! starting from the bottom right corner
	integer i,j,k,jstart,jend,kstart,kend
	! the number of points on surface
	integer, intent(in),dimension(4):: Nsurf
	! the number of points on grid in one edge
	integer, intent(in)::Ngrid
	!the grid coordinate array
	real*8, intent(inout),dimension(Ngrid*Ngrid,2):: x_grid,x_trans
	!the grid coordinate array
	real*8, intent(inout),dimension(sum(Nsurf),2):: x_surf
	
	! the transformed surface and grid
	real*8, allocatable:: x_surtrans(:,:),x_gridsurf(:,:),x_surtransc(:,:),x_gridsurfc(:,:)
	! 1D GLL points and corresponding weight w
	real*8, allocatable::x_1D(:),w(:)
	! all required matrices for mapping
	real*8, allocatable:: A(:,:),Ainv(:,:),Amap(:,:),Apro(:,:)
	
	! the GLL sampled nontransformed surface
	allocate(x_surtrans(4*Ngrid,2))
	! the GLL sampled square surface
	allocate(x_gridsurf(4*Ngrid,2))
	
	! the GLL sampled nontransformed surface without the previous corner
	allocate(x_surtransc(4*Ngrid-4,2))
	! the GLL sampled square surface without the previous corner
	allocate(x_gridsurfc(4*Ngrid-4,2))
	! 1D GLL points and corresponding weight
	allocate(x_1D(Ngrid))
	allocate(w(Ngrid))
	
	call GLLPoints(x_1D,w,Ngrid)
	print*,"okey"
	! Right edge
	x_gridsurf(1:Ngrid,1)=1.
	x_gridsurf(1:Ngrid,2)=x_1D	
	! Top edge
	x_gridsurf((Ngrid+1):2*Ngrid,1)=-1.*x_1D
	x_gridsurf((Ngrid+1):2*Ngrid,2)=1.
	! Left edge
	x_gridsurf((2*Ngrid+1):3*Ngrid,1)=-1.
	x_gridsurf((2*Ngrid+1):3*Ngrid,2)=-1.*x_1D
	! Bottom edge
	x_gridsurf((3*Ngrid+1):4*Ngrid,1)=x_1D
	x_gridsurf((3*Ngrid+1):4*Ngrid,2)=-1.
	
	!call OneDTransform(x_surtrans(1:Ngrid,:),x_surf(1:Nsurf(1),:),x_gridsurf(1:Ngrid,:),Nsurf(1),Ngrid)	
	!do i=1,Ngrid
	!	print*,x_surtrans(i,:)
	!end do
	! let's distribute surfaces
	do i=1,4
		if(i==1)then
			kstart=1
			jstart=1
		else
			kstart=(i-1)*Ngrid+1
			jstart=sum(Nsurf(1:(i-1)))+1
		end if
		jend=sum(Nsurf(1:i))
		kend=Ngrid*i
		call OneDTransform(x_surtrans(kstart:kend,:),x_surf(jstart:jend,:),x_gridsurf(kstart:kend,:),Nsurf(i),Ngrid)	
		! let's fill it into one without the extra corner
		!x_surtransc((kstart-i+1):(kend-i-1),1)=x_surtrans(kstart:(kend-1),1)
		if (i==2) then 
			x_surtransc((kstart-i+1):(kend-i-1),1)=-1.*x_surtrans(kstart:(kend-1),1)
			x_surtransc((kstart-i+1):(kend-i-1),2)=x_surtrans(kstart:(kend-1),2)
		eLse if (i==3) then
				x_surtransc((kstart-i+1):(kend-i-1),1)=x_surtrans(kstart:(kend-1),1)
				x_surtransc((kstart-i+1):(kend-i-1),2)=-1.*x_surtrans(kstart:(kend-1),2)
			 else
				x_surtransc((kstart-i+1):(kend-i-1),1)=x_surtrans(kstart:(kend-1),1)
				x_surtransc((kstart-i+1):(kend-i-1),2)=x_surtrans(kstart:(kend-1),2)
		end if
		x_gridsurfc((kstart-i+1):(kend-i-1),1)=x_gridsurf(kstart:(kend-1),1)
		x_gridsurfc((kstart-i+1):(kend-i-1),2)=x_gridsurf(kstart:(kend-1),2)
	end do
	
	! let's generate the grid
	do i=1,Ngrid		
		if(i==1)then
			kstart=1
		else
			kstart=(i-1)*Ngrid+1
		end if
		kend=Ngrid*i
		! the x coordinate of the grid
		x_grid(kstart:kend,1)=x_1D
		! the y coordinate of the grid
		x_grid(kstart:kend,2)=-1.*x_1D(i)
	end do
	
	
	open(128,file='1DMap.dat',status='unknown')
	open(129,file='OMap.dat',status='unknown')
	open(130,file='OGrid.dat',status='unknown')
	open(131,file='2DMap.dat',status='unknown')
	do i=1,(4*Ngrid-4)
			write(128,*) x_gridsurfc(i,:)
			write(129,*) x_surtransc(i,:)
	!do i=1,4*Ngrid
	!		write(128,*) x_gridsurf(i,:)
	!		write(129,*) x_surtrans(i,:)
	end do 
	
	allocate(A((4*Ngrid-4),(4*Ngrid-4)))
	allocate(Ainv((4*Ngrid-4),(4*Ngrid-4)))
	allocate(Amap(Ngrid*Ngrid,(4*Ngrid-4)))
	allocate(Apro(Ngrid*Ngrid,(4*Ngrid-4)))
	
	!allocate(Amap((4*Ngrid-4),(4*Ngrid-4)))
	!allocate(Apro((4*Ngrid-4),(4*Ngrid-4)))
	
	
	! let's generate the coefficient matrix
	call matgen(A,x_gridsurfc(:,1),x_gridsurfc(:,2),(4*Ngrid-4))
	call invertSVD((4*Ngrid-4),A,Ainv)	
	call matgenxi(Amap,x_gridsurfc(:,1),x_grid(:,1),x_gridsurfc(:,2),x_grid(:,2),Ngrid*Ngrid,(4*Ngrid-4))
	call matmultnon(Amap,Ainv,Apro,Ngrid*Ngrid,(4*Ngrid-4),(4*Ngrid-4))

	!call matgenxi(Amap,x_gridsurfc(:,1),x_gridsurfc(:,1),x_gridsurfc(:,2),x_gridsurfc(:,2),(4*Ngrid-4),(4*Ngrid-4))
	!call matmultnon(Amap,Ainv,Apro,(4*Ngrid-4),(4*Ngrid-4),(4*Ngrid-4))
	
	! let's perform mapping
	! for x coordinate
	call matvectnon(Apro,x_surtransc(:,1),x_trans(:,1),Ngrid*Ngrid,(4*Ngrid-4))
	! for y coordinate
	call matvectnon(Apro,x_surtransc(:,2),x_trans(:,2),Ngrid*Ngrid,(4*Ngrid-4))
	
	open(128,file='1DMap.dat',status='unknown')
	open(129,file='OMap.dat',status='unknown')
	open(130,file='OGrid.dat',status='unknown')
	open(131,file='2DMap.dat',status='unknown')
	do i=1,(4*Ngrid-4)
			write(128,*) x_gridsurfc(i,:)
			write(129,*) x_surtransc(i,:)
	!do i=1,4*Ngrid
	!		write(128,*) x_gridsurf(i,:)
	!		write(129,*) x_surtrans(i,:)
	end do 
	
	do i=1,Ngrid*Ngrid
		write(130,*) x_grid(i,:)
		write(131,*) x_trans(i,:)
	end do
	
	! for x coordinate
	!call matvectnon(Apro,x_surtransc(:,1),x_gridsurfc(:,1),(4*Ngrid-4),(4*Ngrid-4))
	! for y coordinate
	!call matvectnon(Apro,x_surtransc(:,2),x_gridsurfc(:,2),(4*Ngrid-4),(4*Ngrid-4))
	
	
	!do i=1,(4*Ngrid-6)
	!	write(130,*) x_grid(i,:)
	!	write(131,*) x_gridsurfc(i,:)
	!end do
	
	
end subroutine TwoDMeshTransform