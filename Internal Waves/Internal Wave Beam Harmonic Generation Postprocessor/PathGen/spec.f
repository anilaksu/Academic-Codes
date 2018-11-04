subroutine spectrum1d(u,nx,ny,nz,ntot,ndims,sp)
c Calculate average 1-D spectrum in x-direction.
c ************************************************************
c
c  Matrix us assumed to be in spectral space
c      parameter ( nxh=nx/2, nxhp=nxh+1 )
      real sp(nx/2+1,nz)
      complex u(nx/2+1,ny,nz)

      nxh=nx/2
      nyh=ny/2
      nxhp=nxh+1
c
c  mean velocity,  in spectral it's the first point???
c
      write(66,26)
 26   format(//,25x,'plane mean',/)
      do k=1,nz
        kr=nz-k+1   !reverse index for printing from top to bottom
        write(66,11) kr,u(1,1,kr)
      end do
 11   format(1x,i5,2e15.4)
c
c-3-D simulation
      if (ndims.eq.3) then

c skip the zero mode in the x-direction
      do k=1,nz
        do i=2,nxhp
c       do i=1,nxhp
          sp(i,k)=0.
          do j=1,nyh
            sp(i,k)=sp(i,k)+cabs(u(i,j,k))**2
          end do
          do j=(nyh+2),ny
            sp(i,k)=sp(i,k)+cabs(u(i,j,k))**2
          end do
        end do
      end do
c
c  skip the zero mode in the y-direction
      do k=1,nz
        do i=1,1
          sp(i,k)=0.
c         do j=1,nyh
          do j=2,nyh
            sp(i,k)=sp(i,k)+cabs(u(i,j,k))**2
          end do
          do j=(nyh+2),ny
            sp(i,k)=sp(i,k)+cabs(u(i,j,k))**2
          end do
        end do
      end do
c
      do k=1,nz
        do i=1,nxhp
          sp(i,k)=sp(i,k)/float(ny)
        end do
      end do
c

c-2-D simulation (a lot simpler)
      else

c-Zero mode is skipped as in 3-D case
      j = 1
      do k=1,nz
        do i=2,nxhp
          sp(i,k) = sp(i,k)+cabs(u(i,j,k))**2
         enddo
      enddo

      endif
      
      return
      end

        subroutine spectrum2d(u,xsq,ysq,alpha,beta,
     +                        nx,ny,nz,ntot,sp)
c Calculate 2-D spectrum
c ************************************************************
c
c  Matrix us assumed to be in spectral space
c      parameter ( nxh=nx/2, nxhp=nxh+1 )
      parameter(np=20)
      real sp(np+1,nz),xsq(nx/2+1),ysq(ny)
      complex u(nx/2+1,ny,nz)

      nxh=nx/2
      nyh=ny/2
      nxhp=nxh+1
C-aklim is the maximum wavenumber magnitude
      aklim=sqrt((float(nxhp)*alpha)**2.+
     +          (float(nyh)*beta)**2.)    
c-Nulify spectra in all planes
      do k=1,nz
        do m=1,np
          sp(m,k) = 0.
        enddo
      enddo
c
c-ak = the wave# magnitude of (kx,ky)
c-Skip zero mode in x-direction

      do k=1,nz
c       do i=1,nxhp
        do i=2,nxhp
          do j=1,nyh
            ak = sqrt(xsq(i) + ysq(j))
            nbin = int(ak/aklim*float(np))+1
            sp(nbin,k)=sp(nbin,k)+cabs(u(i,j,k))**2
          end do
          do j=(nyh+2),ny
            nbin = int(ak/aklim*float(np))+1
            ak = sqrt(xsq(i) + ysq(j))
            sp(nbin,k)=sp(nbin,k)+cabs(u(i,j,k))**2
          end do
        end do
      end do
c
c-Skip zero mode in y-direction. Make up for all (1,j) points.
      do k=1,nz
        do i=1,1
c         do j=1,nyh
          do j=2,nyh
            ak = sqrt(xsq(i) + ysq(j))
            nbin = int(ak/aklim*float(np))+1
            sp(nbin,k)=sp(nbin,k)+cabs(u(i,j,k))**2
          end do
          do j=(nyh+2),ny
            nbin = int(ak/aklim*float(np))+1
            ak = sqrt(xsq(i) + ysq(j))
            sp(nbin,k)=sp(nbin,k)+cabs(u(i,j,k))**2
          end do
        end do
      end do
c
      return
      end

