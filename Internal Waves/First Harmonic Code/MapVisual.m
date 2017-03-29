
% this routine is generated to visualize 2D solution
 
SimData1=importdata('E_har.dat');
%x=SimData1(:,1);
%y=SimData1(:,2);
AmpFharx=SimData1(:,3);

%let's clean the data
[ AmpFharx ] = cleanData(AmpFharx);

SimData2=importdata('IncidentAmplitude.dat');
x=SimData2(:,1);
y=SimData2(:,2);
AmpIncx=SimData2(:,3);

SimData3=importdata('ReflectingAmplitude.dat');

AmpRefx=SimData3(:,3);

SimData4=importdata('RHSForcing.dat');

RHS=SimData4(:,3);

[ RHS ] = cleanData(RHS);

Nx=40;
Ny=40;

%% let's reshape vectors to plot them

AmpFharx_r = reshape(AmpFharx,Nx,Ny);
AmpIncx_r = reshape(AmpIncx,Nx,Ny);
AmpRefx_r = reshape(AmpRefx,Nx,Ny);
RHS_r = reshape(RHS,Nx,Ny);
x_r = reshape(x,Nx,Ny);
y_r = reshape(y,Nx,Ny);


%% let's generate the countor plot
figure(1)

contourf(x_r,y_r,AmpFharx_r,'LineColor','none')
colorbar
title('The energy flux of the first higher harmonic wave  in \xi_{2\omega_0} direction', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'F_{2\omega_0}^{\xi_{2\omega_0}}/A_{0}^4')
pbaspect([3 1 1])
%% let's generate the countor plot
figure(2)

contourf(x_r,y_r,AmpIncx_r,'LineColor','none')
colorbar
title('The amplitute of the incident velocity  in x direction', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{2\omega_0}^{x}/A_0')
pbaspect([3 1 1])
figure(3)

contourf(x_r,y_r,AmpRefx_r,'LineColor','none')
colorbar
title('The amplitute of the reflecting velocity  in x direction', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{ref}^{x}/A_0')
pbaspect([3 1 1])
figure(4)

contourf(x_r,y_r,RHS_r,'LineColor','none')
colorbar
title('RHS Forcing', 'FontSize', 8)
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{ref}^{x}/A_0')
pbaspect([3 1 1])