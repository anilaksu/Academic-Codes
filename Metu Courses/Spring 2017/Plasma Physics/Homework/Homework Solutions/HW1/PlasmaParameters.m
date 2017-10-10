% created by Anil Aksu to calculate plasma parameters of some ionized gases

clear all 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%       Plasma related constants          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the speed of light 
c=3.*10^8;
% the charge of electron
e=1.602*10^-19;
% the mass of electron
m=9.109*10^-31;
% permitivity of vacuum 
eps=8.85*10^-12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%       Ionized Medium Parameters         %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  let's read them from excel document 
Plasma = xlsread('PlasmaParameters.xlsx');

% the number of ionized gas
Nion=length(Plasma(:,1));
% the debye length array
Debye=zeros(Nion,1);
% the number of particles in debye sphere
webge=zeros(Nion,1);
% the electron plasma frequency
omega_pe=zeros(Nion,1);
% the electron gyro frequency
omega_ce=zeros(Nion,1);

% let's calculate them 
for i=1:Nion
    % the debye length
    [ Debye(i,1) ] = getDebyeLength( eps,Plasma(i,3),Plasma(i,2),e );
    % the number of particles in debye sphere
    wedge(i,1)=Plasma(i,2)*Debye(i,1)^3.;
    % the plasma electron frequency 
    [ omega_pe(i,1) ] = getPlasmaFrequency( eps,Plasma(i,2),e,m );
    % the plasma gyro frequency 
    [ omega_ce(i,1) ] = getGyroFrequency(e,Plasma(i,4),m,c );
end