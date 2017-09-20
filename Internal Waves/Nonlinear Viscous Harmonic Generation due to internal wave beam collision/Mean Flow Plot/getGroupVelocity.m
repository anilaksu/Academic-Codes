function [ Cg ] = getGroupVelocity( N,omega,kx )
% this function generates group velocity vector for given profile
% first let's calculate kz for given profile

[ kz ] = getBVProfile( N,omega,kx );
% the group velocity in x direction
Cg(1)=N/(kx^2.+kz^2.)^0.5-(N*kx^2.)/(kx^2.+kz^2.)^1.5;
% the group velocity in z direction
Cg(2)=-(N*kx*kz)/(kx^2.+kz^2.)^1.5;

end

