function [ E ] = getSecondHarmonic( omega,N,kx_0,sig )


% the angle of incidence 
theta= asin(omega/N);
% the angle of second harmonic wave
gamma=asin(2.*omega/N);
% the gaussian parameters
alpha=(sin(theta)/sig);
beta=(cos(theta)/sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%     The Spectral Amplitude Compuation   %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


syms x y a kx
F_2 = exp(-(x*alpha)^2)*exp(2*i*kx_0*x);
f_2four=fourier(F_2, x,kx);
fz=int(exp(-(beta*x)^2)*sin(cot(gamma)*kx*x), x, -Inf, 0);
% the fourier envelope in x direction
Envelope =symfun(f_2four,[kx]);
fz_hat=symfun(fz,[kx]);
Amp = symfun( 24*i*sin(theta)*Envelope*fz_hat/(sin(2*gamma)*kx),[kx]);
% the integral of envelope
Env_int=double(int(Envelope,kx, 0,100));
% the integral of energy 
%Energy_int=double(int(Energy,kx, 0,100));
% number of wave numbers
Nwave=301;
% the wave number spectrum 
kx_lin=linspace(1.5*kx_0,2.5*kx_0,Nwave);
% the 
Amp_samp=Amp(kx_lin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%         The Energy Integrataion         %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[ E ] = getSpectralEnegry( fz,kx,theta,Nwave );
% f is the spectral density 
E=0.;
for j=1:(Nwave-1)
   E=E+0.5*cos(gamma)*double((kx_lin(j+1)-kx_lin(j))*(kx_lin(j+1)*abs(Amp_samp(j+1))^2+kx_lin(j)*abs(Amp_samp(j))^2));
end

end

