function [ E ] = getIncidentWave( omega,N,kx,sig )

% this function generates group velocity vector for given profile
% first let's calculate kz for given profile

[ kz ] = getBVProfile( N,omega,kx );
% the group velocity in x direction
Cg(1)=N/(kx^2.+kz^2.)^0.5-(N*kx^2.)/(kx^2.+kz^2.)^1.5;
% the group velocity in z direction
Cg(2)=-(N*kx*kz)/(kx^2.+kz^2.)^1.5;

% the magnitude of the group velocity
Cg_mag=(Cg(1)^2.+Cg(2)^2.)^0.5;

% the incident internal wave beam energy
E=Cg_mag*(pi^0.5)*sig*((N*kx)^2./(kz*omega^2.)^2.);

end

