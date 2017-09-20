function [ phi ] = getPhase( kx,kz,omega,x,z,t,phi_0 )
%this function calculates the phase of the complex exponential
%based on k_x x +k_z z-w t+ phi_0
phi=kx*x+kz*z-omega*t+phi_0;
end

