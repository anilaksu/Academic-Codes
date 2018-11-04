function [ u ] = getVelocityIncident( x_trans,A0,sig,alpha,kx,kz )
%this function calculates the velocity field
u(1)=A0*exp(0.5*(sig*alpha)^2.)*(0.5+0.5*erf((x_trans(1)-alpha*sig^2.)/sig*2.^0.5));
u(1)=u(1)*exp(-0.5*(x_trans(2)/sig)^2.)*exp(-alpha*x_trans(1));
% let's eleminate zeroth contour
if(u(1)<10^-2)
    u(1)=0.;
end
u(2)=u(1)*kx/kz;

end

