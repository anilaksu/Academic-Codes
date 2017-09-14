% created by Anil Aksu
% this routine is generated to solve and plot the solution of higher
% harmonics generation by T.R. Akylas

clear all 
format long

% the results by Anil Aksu
SecondAmp = xlsread('AnilResults_v1.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%       Dispersion Relation Parameters    %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the number of frequencies
Nfre=100;
% the BV frequency 
N=2.35;
% the freqency of the internal wave beam 
omega=linspace(0.02*N,0.499*N,Nfre);
% the wave length in x direction
lamx=0.875;
% the wave number in x direction
kx_0=2*pi/lamx;
% the half width 
sig=0.5;

for i=1:Nfre
    E(i,1) = getSecondHarmonic( omega(i),N,kx_0,sig );
    E_inc(i,1) = getIncidentWave( omega(i),N,kx_0,sig );
    E_sec(i,1)=E(i,1)/E_inc(i,1);
end

figure(5)
plot(omega/N,E)
title('Second harmonic energy')
xlabel('\omega/N')
ylabel('F_{2 \omega}/F_{inc}')

figure(6)
plot(omega/N,E_inc)
title('Incident internal wave energy')
xlabel('\omega/N')
ylabel('F_{2 \omega}/F_{inc}')

figure(7)
plot(omega/N,E_sec)
hold on
plot(SecondAmp(:,1),SecondAmp(:,2),'r*')
hold on
plot(SecondAmp(:,1),SecondAmp(:,3),'mo')
hold on
plot(SecondAmp(:,1),SecondAmp(:,4),'gs')
hold on
plot(SecondAmp(:,1),SecondAmp(:,5),'bp')
hold on
plot(SecondAmp(:,1),SecondAmp(:,6),'kh')

title('Normalized second harmonic energy flux')
xlabel('\omega/N')
ylabel('F_{2 \omega}/(A_{0}^2 F_{inc})')
legend('Inviscid Solution by TAL','Viscous Solution at \xi_{2 \omega} = 0.5 \lambda_x','Viscous Solution at \xi_{2 \omega} = \lambda_x','Viscous Solution at \xi_{2 \omega} = 1.5 \lambda_x','Viscous Solution at \xi_{2 \omega} = 2 \lambda_x','Viscous Solution at \xi_{2 \omega} = 2.5 \lambda_x')



 