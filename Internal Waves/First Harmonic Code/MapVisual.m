
% this routine is generated to visualize 2D solution
 
%SimData1=importdata('2DNumerical.dat');
%x=SimData1(:,1);
%y=SimData1(:,2);
%AmpFharx=SimData1(:,3);

SimData2=importdata('IncidentAmplitude.dat');
x=SimData2(:,1);
y=SimData2(:,2);
AmpFharx=SimData2(:,3);

SimData3=importdata('ReflectingAmplitude.dat');

AmpRefx=SimData3(:,3);

Nx=20;
Ny=20;

%% let's reshape vectors to plot them

AmpFharx_r = reshape(AmpFharx,Nx,Ny);
AmpRefx_r = reshape(AmpRefx,Nx,Ny);
x_r = reshape(x,Nx,Ny);
y_r = reshape(y,Nx,Ny);


%% let's generate the countor plot
figure(1)

contourf(x_r,y_r,AmpFharx_r,'LineColor','none')
colorbar
title('The amplitute of the first higher harmonic velocity  in x direction', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{2\omega_0}^{x}/A_0')

figure(2)

contourf(x_r,y_r,AmpRefx_r,'LineColor','none')
colorbar
title('The amplitute of the reflecting velocity  in x direction', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{ref}^{x}/A_0')