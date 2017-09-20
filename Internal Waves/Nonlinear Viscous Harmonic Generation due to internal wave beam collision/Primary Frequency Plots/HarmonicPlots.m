% created by Anil Aksu
% this routine is generated to create plots for internal wave beam problems

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
omega=N/(2^0.5);

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
L=12*lamx;
H=5*lamx;

% the mother interval 
x_mother=linspace(0,L,Nx);
y_mother=linspace(0,H,Ny);

% the center of excitation region
x_ref(1)=2.*lamx;
x_ref(2)=1.3*lamx;

%% the grid
x_grid=zeros(Nx*Ny,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%        Wave Energetics Parameters       %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the maximum amplitude of the x velocity
A0=0.3;
% the half-width of the wave envelope
%sig=0.1014;
sig=0.4014;
% the kinematic viscosity
xnu=10^-5;

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

%% let's generate the countor plot
figure(4)

contourf(x_r,y_r,W_r/A0,'LineColor','none')
colorbar
title('The velocity in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'w_{inc}^{z}/A_0')
%% let's obtain the vicous dissipation along the ray path
% first let's get the point where the centerline of the internal wave beam
% intersects the top boundary 
[ x_int ] = getReflectionPoint( x_ref,kx,kz,H )
Nray=70;

%x_r = getRayPath( Cg,x_ref,Nray,0.5*lamx );
x_r = getRayPath(x_ref,x_int,Nray );
% let's perform the coordinate transformation and velocity calculation
for i=1:Nray
    xr_trans(i,:) = transformCoordinate( x_r(i,:),x_r,Cg );
    u_r(i,:) = getVelocityIncident( xr_trans(i,:),A0,sig,alpha,kx,kz );
end

figure(5)

plot(xr_trans(:,1),u_r(:,1)/max(u_r(:,1)),'LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%    The Reflecting Beam Calculations     %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% let's generate the group velocity
[ Cg ] = getGroupVelocity( N,omega,kx )
% let's perform the coordinate transformation and velocity calculation
for i=1:Nx*Ny
    x_trans(i,:) = transformCoordinate( x_grid(i,:),x_int,Cg );
    AmpRef(i,:) = getVelocityReflecting( x_trans(i,:),u_r(Nray,1),sig,alpha,kx,kz );
    %u(i,:) = getVelocity( x_trans(i,:),A0,sig,alpha,kx,kz );
    % the phase of complex exponential 
    % note that the phase difference phi_0 = -2.*kz*H added
    phiRef = getPhase( kx,kz,omega,x_grid(i,1),x_grid(i,2),time,-2.*kz*H);
    % let's multiply amplitude with cos(phi)
    uRef(i,:)=AmpRef(i,:)*cos(phiRef);
end

%% let's reshape vectors to plot them
URef_r = reshape(uRef(:,1),Ny,Nx);
WRef_r = reshape(uRef(:,2),Ny,Nx);
AmpRefx_r = reshape(AmpRef(:,1),Ny,Nx);
AmpRefz_r = reshape(AmpRef(:,2),Ny,Nx);

x_r = reshape(x_grid(:,1),Ny,Nx);
y_r = reshape(x_grid(:,2),Ny,Nx);


%% let's generate the countor plot
figure(6)

contourf(x_r,y_r,AmpRefx_r/A0,'LineColor','none')
colorbar
title('The amplitute of the velocity of the reflecting wave beam in x direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{ref}^{x}/A_0')

%% let's generate the countor plot
figure(7)

contourf(x_r,y_r,AmpRefz_r/A0,'LineColor','none')
colorbar
title('The amplitute of the velocity of the reflecting wave beam in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{ref}^{z}/A_0')

figure(8)

contourf(x_r,y_r,URef_r/A0,'LineColor','none')
colorbar
title('The velocity of the reflecting wave beam in x direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'u_{ref}^{x}/A_0')

%% let's generate the countor plot
figure(9)

contourf(x_r,y_r,WRef_r/A0,'LineColor','none')
colorbar

title('The velocity of the reflecting wave beam in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'w_{ref}^{z}/A_0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%  The incident and reflecting velocity   %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% the amplitude 
Amp_tot=Amp+AmpRef;
u_tot=u+uRef;
% let's reshape it 
Utot_r = reshape(u_tot(:,1),Ny,Nx);
Wtot_r = reshape(u_tot(:,2),Ny,Nx);
Amptotx_r = reshape(Amp_tot(:,1),Ny,Nx);
Amptotz_r = reshape(Amp_tot(:,2),Ny,Nx);

%% let's generate the countor plot
figure(10)

contourf(x_r,y_r,Amptotx_r/A0,'LineColor','none')
colorbar
title('The amplitute of the velocity of the total wave beam in x direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{tot}^{x}/A_0')

%% let's generate the countor plot
figure(11)

contourf(x_r,y_r,Amptotz_r/A0,'LineColor','none')
colorbar
title('The amplitute of the velocity of the total wave fieldd in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{tot}^{z}/A_0')

figure(12)

contourf(x_r,y_r,Utot_r/A0,'LineColor','none')
colorbar
title('The velocity of the total wave beam in x direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'u_{tot}^{x}/A_0')

%% let's generate the countor plot
figure(13)

contourf(x_r,y_r,Wtot_r/A0,'LineColor','none')
colorbar

title('The velocity of the total wave beam in z direction')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'w_{tot}^{z}/A_0')
