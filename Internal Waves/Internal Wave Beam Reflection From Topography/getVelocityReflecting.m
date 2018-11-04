function [u] = getVelocityReflecting( x_trans,A0,sig,alpha,kx,kz )
%this function calculates the reflecting velocity field
u(1)=A0*exp(-0.5*(x_trans(2)/sig)^2.)*exp(-alpha*x_trans(1));
% let's eleminate zeroth contour
if(u(1)<10^-2)
    u(1)=0.;
end
u(2)=u(1)*kx/kz;

end

