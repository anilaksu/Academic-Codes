function [ gamma,phase] = getPhasesRatios( omega,f,kx,ky,kz,N,g,rho_0 )
% this function computes the phase differences for given carrier wavenumber
% and freuqency, it also computes the dimensiona ratio between the velocity
% vector and the perturbation density and pressure field.
% omega: the frequency of the given wave beam
% f: the coriolis parameter
% kx: the carrier wavenumber in x direction of the given wave beam
% ky: the carrier wavenumber in y direction of the given wave beam
% kz: the carrier wavenumber in z direction of the given wave beam

% the phase differences in radians
phase(1,1)= atan((f*ky)/(omega*kx)); % the phase difference of x velocity
phase(2,1)= atan(-(f*kx)/(omega*ky)); % the phase difference of y velocity
phase(3,1)= 0.; % the phase difference of z velocity

% the gamma coefficients dimensional
gamma(1,1)= -sqrt((kx/kz)^2. + ((f*ky)/(omega*kz))^2.)*(N^2-omega^2)/(omega^2-f^2); % ratio coefficient of x velocity
gamma(2,1)= -sqrt((ky/kz)^2. + ((f*kx)/(omega*kz))^2.)*(N^2-omega^2)/(omega^2-f^2); %  ratio coefficient of y velocity
gamma(3,1)= 1.; % ratio coefficient of z velocity
gamma(4,1)= (rho_0*(omega^2-N^2)/(omega*kz)); %  ratio coefficient of pertubation pressure 
gamma(5,1)= (rho_0*N^2)/(g*omega); %  ratio coefficient of perturbation density
gamma(6,1)= 0.5*rho_0*(gamma(1,1)^2+gamma(2,1)^2+1)*(omega^2)/(omega^2-f^2); %  ratio coefficient of mean energy

end

