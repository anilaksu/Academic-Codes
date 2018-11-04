% created by Anil Aksu
% this routine is generated to compute the scattered electric field by
% circular PEC

clear all 
format long

phii = pi/4; % the angle of incidence
phio = pi/2; % the angle of observation
freq   = 2.*pi;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/(freq*1e6);     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber
a=100/k;   % the radius of PEC
r=a; % the distance
%n_up= 5; % the upper limit of the series


% the error convergence 
for n_up=1:400

a=0.1/k;   % the radius of PEC
% the incident field
E_inc=exp(-i*k*a*cos(phio-phii));
% the scattered electric field
[Escat_1] = fun_cylinder_PEC(freq, a, phii, phio, a,n_up);
% the percent error
err(n_up,1)=abs(real(E_inc+Escat_1)/real(E_inc));


a=10.*a; % the radius of PEC
% the incident field
E_inc=exp(-i*k*a*cos(phio-phii));
% the scattered electric field
[Escat_2] = fun_cylinder_PEC(freq, a, phii, phio, a,n_up);
% the percent error
err(n_up,2)=abs(real(E_inc+Escat_2)/real(E_inc));

a=10.*a; % the radius of PEC
% the incident field
E_inc=exp(-i*k*a*cos(phio-phii));
% the scattered electric field
[Escat_3] = fun_cylinder_PEC(freq, a, phii, phio, a,n_up);
% the percent error
err(n_up,3)=abs(real(E_inc+Escat_3)/real(E_inc));

a=10.*a; % the radius of PEC
% the incident field
E_inc=exp(-i*k*a*cos(phio-phii));
% the scattered electric field
[Escat_4] = fun_cylinder_PEC(freq, a, phii, phio, a,n_up);
% the percent error
err(n_up,4)=abs(real(E_inc+Escat_4)/real(E_inc));

end

figure(1)

semilogy(err(:,1),'LineWidth',2)
hold on
semilogy(err(:,2),'r','LineWidth',2)
hold on 
semilogy(err(:,3),'m','LineWidth',2)
hold on
semilogy(err(:,4),'c','LineWidth',2)

ylabel('the relative error')
xlabel('n number of terms')
title('The relative error for ka=0.1,1,10,100')
legend('ka=0.1','ka=1','ka=10','ka=100')