% created by Anil Aksu
% this routine is generated to create plots for first higher harmonic

clear all 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%       Dispersion Relation Parameters    %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the time of the plot
time=0.4;
% the BV frequency 
N=2.35;
% the freqency of the internal wave beam 
%omega=N/(2^0.5);
omega=0.267*N;
% the wave length in x direction
lamx=0.875;
% the wave number in x direction
kx=2*pi/lamx;
% the wave number in z direction
kz = getBVProfile( N,omega,kx );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%             Grid Parameters             %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the number of points in x direction
Nx=200;
% the number of points in y direction
Ny=200;

% the length of the domain
L=25*lamx;
H=5*lamx;

% the mother interval 
x_mother=linspace(0,L,Nx);
y_mother=linspace(0,H,Ny);

% the center of excitation region
% x_ref(1)=2.*lamx;
% x_ref(2)=1.3*lamx;
x_ref(1)=0.;
x_ref(2)=0.;
%% the grid
x_grid=zeros(Nx*Ny,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%        Wave Energetics Parameters       %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the maximum amplitude of the x velocity
A0=3.;
% the half-width of the wave envelope
%sig=0.1014;
sig=0.4014;
% the kinematic viscosity
xnu=10^-7;

% let's generate the grid
for i=1:Nx
    for j=1:Ny
        % the counter
        count=(i-1)*Ny+j;
        x_grid(count,1)=x_mother(i);
        x_grid(count,2)=y_mother(j);
    end
end

% let's generate the group velocity
[ Cg ] = getGroupVelocity( N,omega,kx )

Cg=abs(Cg);

% the viscous decay coefficient
alpha=xnu*((N*kx)^2)/(2.*(Cg(1)^2.+Cg(2)^2.)*omega^2.);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%     The Incident Beam Calculations      %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% let's perform the coordinate transformation and velocity calculation
for i=1:Nx*Ny
    x_trans(i,:) = transformCoordinate( x_grid(i,:),x_ref,Cg );
    Amp(i,:) = getVelocityIncident( x_trans(i,:),A0,sig,alpha,kx,kz );
    %u(i,:) = getVelocity( x_trans(i,:),A0,sig,alpha,kx,kz );
    % the phase of complex exponential 
    phi = getPhase( kx,-1.*kz,omega,x_grid(i,1),x_grid(i,2),time,0.);
    % let's multiply amplitude with cos(phi)
    u(i,:)=Amp(i,:)*cos(phi);
end

%% let's reshape vectors to plot them
U_r = reshape(u(:,1),Ny,Nx);
W_r = reshape(u(:,2),Ny,Nx);
Ampx_r = reshape(Amp(:,1),Ny,Nx);
Ampz_r = reshape(Amp(:,2),Ny,Nx);
x_r = reshape(x_grid(:,1),Ny,Nx);
y_r = reshape(x_grid(:,2),Ny,Nx);

%figure(1)

%plot(x_grid(:,1)/lamx,x_grid(:,2)/lamx,'*')


%% let's generate the countor plot
figure(1)

contourf(x_r,y_r,Ampx_r/A0,'LineColor','none')
colorbar
title('The amplitute of the velocity in x direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{inc}^{x}/A_0')

%% let's generate the countor plot
figure(2)

contourf(x_r,y_r,Ampz_r/A0,'LineColor','none')
colorbar
title('The amplitute of the velocity in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{inc}^{z}/A_0')

figure(3)

contourf(x_r,y_r,U_r/A0,'LineColor','none')
colorbar
title('The velocity in x direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'u_{inc}^{x}/A_0')
pbaspect([2 1 1])
%% let's generate the countor plot
figure(4)

contourf(x_r,y_r,W_r/A0,'LineColor','none')
colorbar
title('The velocity in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'w_{inc}^{z}/A_0')
pbaspect([2 1 1])
%% let's obtain the vicous dissipation along the ray path
% first let's get the point where the centerline of the internal wave beam
% intersects the top boundary 
[ x_int ] = getReflectionPoint( x_ref,kx,kz,H )
Nray=70;

%x_r = getRayPath( Cg,x_ref,Nray,0.5*lamx );
x_r = getRayPath(x_ref,x_int,Nray );
% let's perform the coordinate transformation and velocity calculation
for i=1:Nray
    xr_trans(i,:) = transformCoordinate( x_r(i,:),x_ref,Cg );
    u_r(i,:) = getVelocityIncident( xr_trans(i,:),A0,sig,alpha,kx,kz );
end

figure(5)

plot(xr_trans(:,1)/lamx,u_r(:,1)/max(u_r(:,1)),'LineWidth',2)

title('The energy flux along the Ray path till the level of free slip surface')
xlabel('\xi /\lambda_x')
ylabel('F(\xi)/F_{max}')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%  The First Higher Harmomic Parameters   %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency 
omega_har=2.*omega;
% the wave number in x direction
kx_har=2.*kx;
% the wave number in z direction
[ kz_har ] = getBVProfile( N,omega_har,kx_har );
% the group velocities
[ Cg_har ] = getGroupVelocity( N,omega_har,kx_har );

%% the coefficients defined in Appendix C of Nonlinear viscous higher generation 
%% by Anil Aksu 2017

% the gamma coeffiecient 
gamma=-2.*(kz_har*omega)^2./(sqrt(Cg_har(1)^2.+Cg_har(2)^2.)*kx*N^2.);
% the viscous decay coefficient for the first higher harmonic
alpha_har=xnu*((N*kx_har)^2)/(2.*(Cg_har(1)^2.+Cg_har(2)^2.)*omega_har^2.);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%  The cofficient from first coordinate   %
%            transfromation               %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the rotation angle for the incident internal wave beam
theta=abs(atan(Cg(2)/Cg(1)));

% the amplitude parameter
    D_har=u_r(Nray,1)*A0*exp(alpha*(x_int(2)*sin(theta)+x_int(1)*cos(theta)));
% the argument in exponential 
arg1=(-x_int(1)*x_int(2)*sin(2.*theta)-(x_int(1)*sin(theta))^2.-(x_int(2)*cos(theta))^2.)/(2.*sig^2.);
% the coordinate mapping parameters
alpha_x=-1.*(sin(theta)/sig)^2.;
beta_x=x_int(1)*(sin(theta)/sig)^2.+0.5*x_int(2)*sin(2.*theta)/(sig^2.);
%-2.*alpha*cos(theta);
alpha_z=-1.*(cos(theta)/sig)^2.;
beta_z=x_int(2)*(cos(theta)/sig)^2.+0.5*x_int(1)*sin(2.*theta)/(sig^2.);

% the right handside forcing of first higher harmonic
F_har=zeros(Nx*Ny,1)
for i=1:Nx*Ny
    
    % the function inside the exponential part of the forcing
    arg2=(alpha_x*x_grid(i,1)+beta_x)*x_grid(i,1)+alpha_z*x_grid(i,2)^2+beta_z*x_grid(i,2);
    % the exponential part of first higher harmonic
    F_har(i,1)=D_har*exp(arg2+arg1)*sin(kz_har*x_grid(i,2)+2*kz*x_int(2));
    
    % let's get rid of the zeroth contour
    if( abs(F_har(i,1))<10^-2)
        F_har(i,1)=0.;
    end 
end

%% let's reshape vectors to plot them
Fhar_r = reshape(F_har(:,1),Nx,Ny);
x_r = reshape(x_grid(:,1),Nx,Ny);
y_r = reshape(x_grid(:,2),Nx,Ny);

%% let's generate the countor plot for the first harmonic
figure(6)

contourf(x_r,y_r,Fhar_r/A0^2,'LineColor','none')
colorbar
pbaspect([2 1 1])
title('The RHS forcing field of the first higher harmonic')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'F_{2\omega}/A_{0}^{2}')

% let's generate the cross-sectional profile of the first higher harmonic 
% velocity field
H_har=D_har*exp(alpha_x*x_int(1)^2.+beta_x*x_int(1)+alpha_z*x_int(2)^2.+beta_z*x_int(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%  The cofficient from second coordinate  %
%            transfromation               %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % the rotation angle for the incident internal wave beam
 theta_har=atan(Cg_har(2)/Cg_har(1));
% 
% % the updated argument parameter
% arg=arg+alpha_x*x_int(1)^2.+beta_x*x_int(1)+alpha_z*x_int(2)^2.+beta_z*x_int(1)
% % the amplitude parameter
% G_har=D_har*exp(arg);
% 
% the coordinate mapping parameters
alpha_xi=alpha_x*cos(theta_har)^2.+alpha_z*sin(theta_har)^2.;
beta_xi=2.*alpha_x*cos(theta_har)*x_int(1)+beta_x*cos(theta_har)+2.*alpha_z*sin(theta_har)*x_int(2)+beta_z*sin(theta_har);
alpha_eta=alpha_x*sin(theta_har)^2.+alpha_z*cos(theta_har)^2.;
beta_eta=-2.*alpha_x*sin(theta_har)*x_int(1)-beta_x*sin(theta_har)+2.*alpha_z*cos(theta_har)*x_int(2)+beta_z*cos(theta_har);
ksi=2.*(alpha_z-alpha_x)*sin(theta_har)*cos(theta_har);

% let's generate the cross-sectional profile of the first higher harmonic 
% velocity field
G_har=D_har*exp(alpha_x*x_int(1)^2.+beta_x*x_int(1)+alpha_z*x_int(2)^2.+beta_z*x_int(2));
arg3=(kz_har*sin(theta_har))^2./(4.*alpha_xi);
H_har=D_har*gamma*sqrt(-1.*pi/alpha_xi);
%*exp(arg1+arg3);
% the number of points along the wave number direction
N_xi=100;
% the cooridantes in the first harmonic rotation
xi_har=linspace(0.5*lamx,1.5*lamx,3);
eta_har=linspace(-2.*lamx,2.*lamx,N_xi);
% the velocity profile
u_har=zeros(N_xi,3);
v_har=zeros(N_xi,3);

for i=1:3
    for j=1:N_xi
        % the argument inside the exponential
        arg_exp=alpha_eta*eta_har(j)^2.+beta_eta*eta_har(j);
        arg_exp=arg_exp-(beta_xi+alpha_har+ksi*eta_har(j))^2./(4.*alpha_xi);
        % arg_exp=arg_exp+(kz_har*sin(theta_har))^2./(4.*alpha_xi);
        % the argument inside the sinusoidal function
        arg_sin=kz_har*cos(theta_har)*eta_har(j)+kz_har*x_int(2)-2.*kz*x_int(2);
        arg_sin=arg_sin-kz_har*sin(theta_har)*(beta_xi+alpha_har+ksi*eta_har(j))/(2.*alpha_xi);
        u_har(j,i)=H_har*exp(arg_exp)*sin(arg_sin);
    end
end

figure(7)

plot(eta_har,u_har(:,1)/A0^2,'LineWidth',4)

title('The first higher harmonic wave at \xi_{2\omega_0}/\lambda_x =1')
xlabel('x/\lambda_x')
ylabel('A_{2\omega}^{x}/A_{0}^2')


% %let's perform the coordinate transformation and velocity calculation
% for i=1:Nx*Ny
%     x_trans(i,:) = transformCoordinate( x_grid(i,:),x_int,Cg_har );
%     AmpFHar(i,:) = getVelocityFirstHarmonic( x_trans(i,:),x_int,gamma,alpha_har,G_har,alpha_xi,beta_xi,alpha_eta,beta_eta,ksi,kx_har,kz_har,-2.*kz*H );
%     % the phase of complex exponential 
%     phi_har = getPhase( kx_har,kz_har,omega_har,x_grid(i,1),x_grid(i,2),time,0.);
%     % let's multiply amplitude with cos(phi)
%     u_har(i,:)=AmpFHar(i,:);
%     %*cos(phi_har);
% end

% %% let's reshape vectors to plot them
% Uhar_r = reshape(u_har(:,1),Ny,Nx);
% Whar_r = reshape(u_har(:,2),Ny,Nx);
% AmpFHarx_r = reshape(AmpFHar(:,1),Ny,Nx);
% AmpFHarz_r = reshape(AmpFHar(:,2),Ny,Nx);
% 
% x_r = reshape(x_grid(:,1),Ny,Nx);
% y_r = reshape(x_grid(:,2),Ny,Nx);
% 
% %% let's generate the countor plot for the first harmonic
% figure(14)
% 
% contourf(x_r,y_r,Whar_r/A0,'LineColor','none')
% colorbar
% 
% title('The velocity of the first higher harmonic in z direction')
% xlabel('x/\lambda_x')
% ylabel('z/\lambda_x')
% hcb=colorbar
% title(hcb,'w_{2\omega}^{z}/A_0')