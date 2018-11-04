function [Escat] = fun_cylinder_PEC(freq, a, phii, phio, r,n_up)
% ***********************************************************************
% Scattered field computation for a PEC circular cylinder (series solution)
% by Dr. Ozlem Ozgun & Dr. Mustafa Kuzuoglu
% ***********************************************************************

% INPUTS:
% freq: frequency in MHz
% a   : radius of the cylinder (m)
% phii: Angle of incidence (radian)
% phio: Angle of observation point (radian)
% r   : Radius of observation point (distance) (m)

% OUTPUT:
% Escat: scattered field

freq   = freq*1e6;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/freq;     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber

phii = phii-pi;  

ka = k*a;
kr = k*r;

E0 = besselj(0,ka)*besselh(0,2,kr)/besselh(0,2,ka);
n = 1:n_up; % terms in the summation
E = sum(2*((-1j).^n).*cos(n*(phio-phii)).*besselj(n,ka).*besselh(n,2,kr)./besselh(n,2,ka));

Escat = -(E0+E); % scattered field

end
