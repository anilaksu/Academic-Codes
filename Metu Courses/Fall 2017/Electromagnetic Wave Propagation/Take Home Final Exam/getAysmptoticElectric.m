function [ E_asym ] = getAysmptoticElectric(freq, a, phii, phio, r)
%This function computes the asymptotic integral given in EE 523 final

phid = linspace(-0.5*pi,.5*pi,100); % the angle of integral
dphi=pi/100.; % the step size of the integral
freq   = freq*1e6;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/freq;     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber
phii = phii-pi;  

ka = k*a;
% the integral 
int=0.;

for i=1:(length(phid)-1)
    f(1) = getIntFunc( phid(i), k, a, phii, phio, r);
    f(2) = getIntFunc( phid(i+1), k, a, phii, phio, r);
   % the trapezoidal rule integral
   int=int+0.5*(f(1)+f(2))*dphi;
end
% the other terms multiplication
alpha=0.82046;
% the asymptotic electric field
E_asym=-0.5*alpha*a*(k^2)*int*sqrt(2*i/(pi*k*r));
end

