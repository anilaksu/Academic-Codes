function [ sigma ] = getRCS( freq, a, phii, phio,n_up )
% this function calculates RCS value at the far field, the formulation
% adopted here is based on the lecture notes of EE 523

% INPUTS:
% freq: frequency in MHz
% a   : radius of the cylinder (m)
% phii: Angle of incidence (radian)
% phio: Angle of observation point (radian)
% r   : Radius of observation point (distance) (m)

% OUTPUT:
% sigma: scattered field

freq   = freq*1e6;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/freq;     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber

phii = phii-pi;  

ka = k*a;

E0 = besselj(0,ka)/besselh(0,2,ka);
n = 1:n_up; % terms in the summation
E = abs(sum(2*cos(n*(phio-phii)).*besselj(n,ka)./besselh(n,2,ka)));

% RCS value 
sigma = (4/k)*E^2;

end

