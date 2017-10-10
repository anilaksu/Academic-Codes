% This script is developed to produce results for Analytical Shock given by
% B. M. Johnson on JFM Rapids and the code is developed by Anil Aksu 
format long;
% the inlet velocity
v_0=100.;
% the exit velocty 
v_1=30.;

% the thermal constants
Lx=1;
% the heat capacitance
cp=1003.;
% the initial temperature
T_0=200.;
% Gamma Constans
gamma=1.4;
% alpha constant
alpha=2.*Lx/(1.+gamma);
% the domain 
% the number of points
Nx=100;
% the error margin
eps=10^-1.;
% the velocity
x=linspace(-5,6,Nx);
%v=linspace(v_0,v_1+eps,Nx);
% the velocity and temperature profile
v=zeros(Nx,1);
T=zeros(Nx,1);
%x=zeros(Nx,1);
R=(gamma+1)/(gamma-1);
%for i=1:Nx
 % x(i,1) = getClosedFormShock( alpha,v(i),v_0,v_1 )
%end
for i=1:Nx
    % the velocity at each location 
    v(i,1)  = getShockVelocity( alpha,v_0,v_1,x(i) );
    % the temperature at each location
    T(i,1)=(R*v_0*v_1-v(i,1)^2.)/(2.*cp); 
end

figure(1)
plot(x/Lx,v/v_0,'LineWidth',2)
 grid on
 title('The velocity profile')
 ylabel('v/v_{0}')
 xlabel('x/L_{\kappa}')
 set(gca,'fontsize',16)
 
 figure(2)
 plot(x/Lx,T/T(1,1),'r','LineWidth',2)
 grid on
 title('The temperature profile')
 ylabel('T/T_{0}')
 xlabel('x/L_{\kappa}')
 set(gca,'fontsize',16)