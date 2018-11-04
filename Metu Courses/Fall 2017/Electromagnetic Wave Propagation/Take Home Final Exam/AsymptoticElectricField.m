% created by Anil Aksu
% this routine is generated to compute the scattered electric field by
% circular PEC

clear all 
format long

phii = pi; % the angle of incidence
phio = linspace(-0.5*pi,.5*pi); % the angle of observation
freq   = 2.*pi;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/(freq*1e6);     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber
a=100/k;   % the radius of PEC
r=10.*a; % the distance
n_up=200;  % the upper limit of the series
%eta=

for i=1:length(phio)

% the scattered electric field given as series
[Escat_series(i)] = fun_cylinder_PEC(freq, a, phii, phio(i), r,n_up);
[ E_asym(i) ] = getAysmptoticElectric(freq, a, phii, phio(i), r);
end

figure(1)
polarplot(phio,abs(Escat_series),'r','LineWidth',2)
hold  on
polarplot(phio,abs(E_asym),'LineWidth',2)
title('Absolute Value of Scattered Electric Field for ka=100 at r=10a')
legend('Series Solution','Asymptotic Solution')

figure(2)
polarplot(phio,2.*pi*abs(Escat_series).^2,'r','LineWidth',2)
hold  on
polarplot(phio,2.*pi*abs(E_asym).^2,'LineWidth',2)
title('RCS for ka=100 at r=10a')
legend('Series Solution','Asymptotic Solution')