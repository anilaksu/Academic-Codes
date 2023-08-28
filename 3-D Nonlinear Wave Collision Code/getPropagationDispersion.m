function [ P,D ] = getPropagationDispersion( N,f,kx,ky,sig )
% this function calculates the frequency, the group velocity and zeta direction 
% for given wave number and the stratification profile
    % the BV frequency N
    % the coriolis frequency f
    % the horizonntal wave number wavenumber vector (kx,ky) 
    % sig: the characteristic width of the envelope 
    % P : the measure of the propagation
    % D : the measure of the dispersion
% the epsilon term
eps=10^-6;
syms k_x k_y k_z omega w

% the frequency
omega = sqrt(((N^2)*(k_x^2+k_y^2)+(f^2)*(k_z^2))/(k_x^2+k_y^2+k_z^2));
omega = symfun(omega,[k_x k_y k_z]);
% the vertical wavenumber
kz = sqrt((N^2-w^2)*(k_x^2+k_y^2)/(w^2-f^2));
kz = symfun(kz,[k_x k_y w]);
% the  group velocity
Cg_x = symfun(diff(omega,k_x),[k_x k_y k_z]);
Cg_y = symfun(diff(omega,k_y),[k_x k_y k_z]);
Cg_z = symfun(diff(omega,k_z),[k_x k_y k_z]);
% the dispersion matrix
D_ij= [diff(Cg_x,k_x) diff(Cg_x,k_y) diff(Cg_x,k_z); diff(Cg_y,k_x) diff(Cg_y,k_y) diff(Cg_y,k_z); diff(Cg_z,k_x) diff(Cg_z,k_y) diff(Cg_z,k_z)];
D_ij = symfun(D_ij,[k_x k_y k_z]);
% the sample of frequencies
omega_1=linspace(f+eps,N-eps,100);

for i=1:length(omega_1)
    % the corresponding wavenumber in vertical 
    kz_1=kz(kx,ky,omega_1(i));
    P(i)= norm(double([Cg_x(kx,ky,kz_1) Cg_y(kx,ky,kz_1) Cg_z(kx,ky,kz_1)]))/(sig);
    D(i)= double(max(abs(eig(D_ij(kx,ky,kz_1))))/((sig)^2));

end

end

