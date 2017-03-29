function [ kz ] = getBVProfile( N,omega,kx )
%this function calculates kx the wave number vector in z direction 
% Note that kz=kx(N^2/w^2-1)^0.5
% the dispersion relation is given as w^2=N^2kx^2/(kx^2+kz^2)
kz=kx*(N^2/omega^2-1)^0.5;

end

