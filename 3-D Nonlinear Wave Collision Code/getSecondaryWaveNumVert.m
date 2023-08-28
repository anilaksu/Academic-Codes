function [ kz_su,kz_sd ] = getSecondaryWaveNumVert(N,f,omega_s, kx_s,ky_s )
% this function computes the wavenumber in vertical direction both for
% upward and downward propagating wave beam 

% the vertical wavenumber for downward propagating beam 
kz_sd = sqrt((N^2-omega_s^2)*(kx_s^2+ky_s^2)/(omega_s^2-f^2));
% the vertical wavenumber for downward propagating beam 
kz_su = -sqrt((N^2-omega_s^2)*(kx_s^2+ky_s^2)/(omega_s^2-f^2));

end

