function [ omega_1, Cg_1, zeta_1, J ] = getFrequency( N,f,kx,ky,kz )
% this function calculates the frequency, the group velocity and zeta direction 
% for given wave number and the stratification profile
    % the BV frequency N
    % the coriolis frequency f
    % the wavenumber vector (kx,ky,kz) 
% the frequency
omega_1 = sqrt(((N^2)*(kx^2+ky^2)+(f^2)*(kz^2))/(kx^2+ky^2+kz^2));

syms k_x k_y k_z 

% the frequency
omega = sqrt(((N^2)*(k_x^2+k_y^2)+(f^2)*(k_z^2))/(k_x^2+k_y^2+k_z^2));

omega= symfun(omega,[k_x k_y k_z]);

Cg_x = symfun(diff(omega,k_x),[k_x k_y k_z]);
Cg_y = symfun(diff(omega,k_y),[k_x k_y k_z]);
Cg_z = symfun(diff(omega,k_z),[k_x k_y k_z]);
% the wave number vector
k_1=[kx ky kz];
% the group velocity vector
Cg_1= double([Cg_x(kx,ky,kz) Cg_y(kx,ky,kz) Cg_z(kx,ky,kz)]);
% the zeta direction
zeta_1= cross(Cg_1, k_1);
% let's form the transformation matrix
J=[Cg_1/norm(Cg_1); k_1/norm(k_1); zeta_1/norm(zeta_1)];

end