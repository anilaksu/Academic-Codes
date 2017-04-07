% this routine is generated to visualize 2D solution
 
SimData1=importdata('2DNumerical.dat');
x=SimData1(:,1);
y=SimData1(:,2);
AmpFharx=SimData1(:,3);



Nx=20;
Ny=20;

%% let's reshape vectors to plot them

AmpFharx_r = reshape(AmpFharx,Nx,Ny);
x_r = reshape(x,Nx,Ny);
y_r = reshape(y,Nx,Ny);


%% let's generate the countor plot
figure(1)

contourf(x_r,y_r,AmpFharx_r,'LineColor','none')
colorbar
title('A sample solution', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'F')
%pbaspect([3 1 1])
