function [ alpha ] = getFourierCoeff(knx,kny,ky,H,D,A,L)
%this function calculates fourier coefficients
% kny is the wave number of each mode
% ky is the wave number of excitation 
% D is the distance between each exciation source
% A is the excitation amplitude
% H is the domain height 
% L is the domain length 
% the wave length of excitation 
lambda=2*pi/ky;
% the fourier coefficients
alpha=cos((kny+ky)*(D-lambda/2))-cos((kny+ky)*(D+lambda/2))+cos((kny-ky)*(D+lambda/2));
alpha=alpha-cos((kny-ky)*(D-lambda/2));
alpha=(-2.*A/L)*alpha/cos(pi/2-knx*L);
end

