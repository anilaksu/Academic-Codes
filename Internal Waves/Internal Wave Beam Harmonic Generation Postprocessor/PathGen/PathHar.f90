! the number of the reflections for fisrt harmonic

RefH1=12
 

 call PathPoints(Nh1,z0,z(N(1),1),nzr,2.*kx,r,dt,ds,N1,2.*w0)

print*,"NumPoint",Nh1

allocate(xh1(Nh1*RefH1,il))
allocate(zh1(Nh1*RefH1,il))


allocate(Cgxh1(Nh1*RefH1,il))
allocate(Cgzh1(Nh1*RefH1,il))

do i=1,il

  call RayPath(Nh1,z0,z(N(1),1),x(N(i),i),2.*kx,r,dt,ds,N1,2.*w0,xh1(:,i),zh1(:,i),RefH1,Cgxh1(:,i),Cgzh1(:,i))

 end do 

open(620,file='RayH11.dat',status='unknown')

do i=1,Nh1*RefH1

 write(620,*) xh1(i,1)/lambdax,zh1(i,1)/lambdax

end do


open(621,file='RayH15.dat',status='unknown')

do i=1,Nh1*RefH1

 write(621,*) xh1(i,5)/lambdax,zh1(i,5)/lambdax

end do


open(622,file='RayH19.dat',status='unknown')

do i=1,Nh1*RefH1

 write(622,*) xh1(i,9)/lambdax,zh1(i,9)/lambdax

end do


open(623,file='RayH113.dat',status='unknown')

do i=1,Nh1*RefH1

 write(623,*) xh1(i,13)/lambdax,zh1(i,13)/lambdax

end do


open(624,file='RayH117.dat',status='unknown')

do i=1,Nh1*RefH1

 write(624,*) xh1(i,17)/lambdax,zh1(i,17)/lambdax

end do


! path for second harmonic


! the number of the reflections for fisrt harmonic

RefH2=18
 

 call PathPoints(Nh2,z0,z(N(1),1),nzr,3.*kx,r,dt,ds,N1,3.*w0)

print*,"NumPoint",Nh2

allocate(xh2(Nh2*RefH2,il))
allocate(zh2(Nh2*RefH2,il))


allocate(Cgxh2(Nh2*RefH2,il))
allocate(Cgzh2(Nh2*RefH2,il))

do i=1,il

  call RayPath(Nh2,z0,z(N(1),1),x(N(i),i),3.*kx,r,dt,ds,N1,3.*w0,xh2(:,i),zh2(:,i),RefH2,Cgxh2(:,i),Cgzh2(:,i))

 end do 

open(632,file='RayH21.dat',status='unknown')

do i=1,Nh2*RefH2

 write(632,*) xh2(i,1)/lambdax,zh2(i,1)/lambdax

end do

open(633,file='RayH25.dat',status='unknown')

do i=1,Nh2*RefH2

 write(633,*) xh2(i,5)/lambdax,zh2(i,5)/lambdax

end do

open(634,file='RayH29.dat',status='unknown')

do i=1,Nh2*RefH2

 write(634,*) xh2(i,9)/lambdax,zh2(i,9)/lambdax

end do

open(635,file='RayH213.dat',status='unknown')

do i=1,Nh2*RefH2

 write(635,*) xh2(i,13)/lambdax,zh2(i,13)/lambdax

end do

open(636,file='RayH217.dat',status='unknown')

do i=1,Nh2*RefH2

 write(636,*) xh2(i,il)/lambdax,zh2(i,il)/lambdax

end do
