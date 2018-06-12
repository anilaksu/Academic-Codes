
% Copyright (c) 2018 Anil A. Aksu 2018                                        
% Developed by Anil A. Aksu, Tubitak Space Technologies Research Institute
% All rights are reserved
% Contributors:                                                             
% Ozgur Gundoğan  
% this routine is generated to create internal wave beam field interacting with 
% the topography, the results and the analysis is submitted to Journal of Computational Physics 

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
omega=0.39*N;
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
L=20*lamx;
H=4.5*lamx;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %
%         Agnesi Properties            %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the amplitude of agnesi
A_agnesi=0.03;
% the location of agnesi
x_0=L/2;
% the deviation of agnesi
sigma=L/16;

% % the plot of agnesi
% figure(9)
% plot(x_mother,f)
% title('The derivative of Witch of Agnesi')
% xlabel('df/dx')
% ylabel('z/\lambda_x')
% the dispacement due to disturbace
% Ampdistx=zeros(Nx,Ny) ;
% Ampdistz_r=zeros(Nx,Ny) ;
% % the deviation amplitude due to distortion 
% for i=1:Ny
%     for j=1:Nx
%         Ampdistz_r(i,j)=f(j,1)*Amptotx_r(i,j);
%         Ampdistx_r(i,j)= Ampdistz_r(i,j)*kz/kx;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%     The Incident Beam Calculations      %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% let's perform the coordinate transformation and velocity calculation
for i=1:Nx*Ny
    x_trans(i,:) = transformCoordinate( x_grid(i,:),x_ref,Cg );
    Amp(i,:) = getVelocityIncident( x_trans(i,:),A0,sig,alpha,kx,kz );
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
    % the incident disturbance
     % the derivative of agnesi for the incident internal wave beam
    f = getTumorAgnesi( A_agnesi,sigma,x_0,x_int,x_trans(i,:),kx,kz);
    % the disturbance due to incident internal wave beam in z direction
    Ampdist_inc(i,2)=f*Amp(i,1);
    % the disturbance due to incident internal wave beam in x direction
    Ampdist_inc(i,1)= Ampdist_inc(i,2)*kz/kx;
    % the reflecting disturbance
     % the derivative of agnesi for the incident internal wave beam
    f = getTumorAgnesi( A_agnesi,sig,x_0,x_int,x_trans(i,:),-kx,kz);
    % the disturbance due to incident internal wave beam in z direction
    Ampdist_ref(i,2)=f*AmpRef(i,1);
    % the disturbance due to incident internal wave beam in x direction
    Ampdist_ref(i,1)= Ampdist_ref(i,2)*kz/kx;
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
AmpIdistx_r = reshape(Ampdist_inc(:,1),Ny,Nx);
AmpIdistz_r = reshape(Ampdist_inc(:,2),Ny,Nx);
AmpRdistx_r = reshape(Ampdist_ref(:,1),Ny,Nx);
AmpRdistz_r = reshape(Ampdist_ref(:,2),Ny,Nx);

x_r = reshape(x_grid(:,1),Ny,Nx);
y_r = reshape(0-x_grid(:,2),Ny,Nx);

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
figure(1)


subplot(2,2,1)   
contourf(x_r,y_r,Amptotx_r/A0,'LineColor','none')
colorbar
%title('The amplitute of the velocity of the total wave beam in x direction')
title('(a)')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{tot}^{x}/A_0')

%% let's generate the countor plot
subplot(2,2,2) 
contourf(x_r,y_r,Amptotz_r/A0,'LineColor','none')
colorbar
%title('The amplitute of the velocity of the total wave fieldd in z direction')
title('(b)')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{tot}^{z}/A_0')

subplot(2,2,3) 
contourf(x_r,y_r,AmpIdistx_r/A0,'LineColor','none')
colorbar
%title('The amplitute of the disturbance velocity in x direction')
title('(c)')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{tot}^{x}/A_0')

%% let's generate the countor plot
subplot(2,2,4) 
contourf(x_r,y_r,AmpIdistz_r/A0,'LineColor','none')
colorbar
%title('The amplitute of the disturbance velocity in z direction')
title('(d)')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'A_{tot}^{z}/A_0')


