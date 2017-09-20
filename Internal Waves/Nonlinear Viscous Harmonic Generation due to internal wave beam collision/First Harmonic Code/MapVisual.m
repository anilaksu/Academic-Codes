
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


SimData5=importdata('E_int.dat');

EnergyFlux=SimData5(:,3);
eta=SimData5(:,2)

SimData6=importdata('E_tot.dat');

EnergyFluxTot=SimData6(:,2);
xi=SimData6(:,1)
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
title(hcb,'F_{2\omega_0}/A_{0}^2')
pbaspect([3 1 1])

u_limit=40;
l_limit=31;
figure(5)

plot(eta(l_limit:u_limit),EnergyFlux(l_limit:u_limit),'LineWidth',2)
hold on
plot(eta(l_limit:u_limit),EnergyFlux(l_limit:u_limit),'r*')

title('The first higher harmonic wave at \xi_{2\omega_0}/\lambda_x =1')
xlabel('x/\lambda_x')
ylabel('A_{2\omega}^{x}/A_{0}^2')

figure(6)

plot(xi,EnergyFluxTot,'LineWidth',4)
%hold on
%plot(xi,EnergyFluxTot,'r*')

title('The first higher harmonic wave energy flux along the Ray path')
xlabel('\xi_{2\omega_0}/\lambda_x')
ylabel('F_{2\omega}^{tot}/A_{0}^4')