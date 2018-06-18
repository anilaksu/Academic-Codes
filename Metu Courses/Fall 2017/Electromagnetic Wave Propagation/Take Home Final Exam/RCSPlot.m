% created by Anil Aksu
% this routine is generated to RCS plot by circular PEC

clear all 
format long

phii = 0; % the angle of incidence
phio = linspace(0,2*pi);; % the angle of observation
freq   = 2.*pi;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/(freq*1e6);     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber
a=100/k;   % the radius of PEC
r=a; % the distance
n_up=400;  % the upper limit of the series

for i=1:length(phio)

% the RCS value 
RCS_1(i)= getRCS( freq, a, phii, phio(i),n_up )/lambda;

end

figure(1)
polarplot(phio,RCS_1)

title('RCS Plot \sigma/\lambda for ka=100')

