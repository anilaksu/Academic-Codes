function [ u ] = getVelocityFirstHarmonic( x_trans,x_int,gamma,alpha_har,G_har,alpha_xi,beta_xi,alpha_eta,beta_eta,ksi,kx_har,kz_har,phi_0 )
%this function calculates the velocity field for the first higher harmonic
% the coordinate transformation angle
theta_har=atan(-1.*kx_har/kz_har);
% the amplitude parameter
H_har=G_har*gamma*sqrt(pi/abs(alpha_xi))*exp((kz_har*sin(theta_har))^2./(4.*alpha_xi));
% the argumment in sinusoidal part
varphi=kz_har*cos(theta_har)*x_trans(2)+kz_har*x_int(2)-phi_0-kz_har*(beta_xi+alpha_har+ksi*x_trans(2))/(2.*alpha_xi);
% the argument in exponential part
phi=alpha_eta*x_trans(2)^2.+beta_eta*x_trans(2)-(beta_xi+alpha_har+ksi*x_trans(2))/(4.*alpha_xi)-(kz_har*sin(theta_har))^2./(4.*alpha_xi);

% let's calculate the amplitude finally
u(1)=phi;
%*sin(varphi);
% let's eleminate zeroth contour
% if(u(1)<10^-4)
%     u(1)=0.;
% end
u(2)=-1.*u(1)*kx_har/kz_har;

end

