% created by Anil Aksu
% this routine is generated to plot the analytical solution for the small
% disturbance on the top of rectangular domain which is meant to model a
% brain tumor 

clear all 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Excitation Parameters          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the freqency of the excitation
omega=0.176;
% the wave length in z direction
lamz=0.875/2.;
% the wave number in z direction
kz=2*pi/lamz;
% the tumor amplitude
A=0.01;
% the distance between the excitation regions
D=8.*lamz;


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
L=20*lamz;
H=20*lamz;
% the location of tumor
x_0=L/4.;

% the mother interval 
x_mother=linspace(0,L,Nx);
y_mother=linspace(-H/2,H/2,Ny);


%% the grid
x_grid=zeros(Nx*Ny,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%        Wave Energetics Parameters       %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the maximum amplitude of the excitation
A0=0.3;
% the half-width of the wave envelope
%sig=4.014;
sig=5.;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                         %
% %          Elasticity Parameters          %
% %                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % lame constant 
% lambda=40.38*10^-5;
% % the displacement field in
% u_x=zeros(Nx*Ny,2);
% 
% 
% % let's generate the grid
% for i=1:Nx
%     for j=1:Ny
%         % the counter
%         count=(i-1)*Ny+j;
%         x_grid(count,1)=x_mother(i);
%         x_grid(count,2)=y_mother(j);
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                         %
% %     Generated Displacement Field        %
% %                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [ u_x, k_ny, k_nx, alpha_n ] = PrimaryDispField( omega,lamz,D,A0,L,H,Nx,Ny );
% % tumor itself
% [ f ] = getTumor( A,sig,x_0,L,Nx );
  [ f ] = getTumorAgnesi( A,sig,x_0,L,Nx );
% % disturbance due to tumor
% [ u_xt ] = getTumorDisturbace( A,f,u_x,L,H,omega,Nx,Ny);

% the combined function 
[u_x,res,x_grid] = mainTumor( A,A0,sig,lamz,D,L,H,omega,Nx,Ny,x_0);
%% let's reshape vectors to plot them
U_r = reshape(u_x(:,1),Ny,Nx);
V_r = reshape(u_x(:,2),Ny,Nx);
% Ut_r = reshape(u_xt(:,1),Ny,Nx);
% Vt_r = reshape(u_xt(:,2),Ny,Nx);
Ut_r = reshape(res(:,1),Ny,Nx);
Vt_r = reshape(res(:,2),Ny,Nx);
x_r = reshape(x_grid(:,1),Ny,Nx);
y_r = reshape(x_grid(:,2),Ny,Nx);

%% let's generate the countor plot
figure(1)

contourf(x_r/lamz,y_r/lamz,U_r/A0,'LineColor','none')
colorbar
title('The displacement in x direction')
xlabel('x/\lambda_z')
ylabel('y/\lambda_z')
hcb=colorbar
title(hcb,'u_{x}/A_0')

%% let's generate the countor plot
figure(2)

contourf(x_r/lamz,y_r/lamz,V_r/A0,'LineColor','none')
colorbar
title('The displacement in y direction')
xlabel('x/\lambda_z')
ylabel('y/\lambda_z')
hcb=colorbar
title(hcb,'u_{y}/A_0')

% figure(3)
% plot(alpha_n)
% title('Fourier Coefficients')
% xlabel('n')
% ylabel('\lambda_n')

figure(4)

contourf(x_r/lamz,y_r/lamz,Ut_r/A0,'LineColor','none')
colorbar
title('The displacement due to tumor in x direction')
xlabel('x/\lambda_z')
ylabel('y/\lambda_z')
hcb=colorbar
title(hcb,'u_{x}/A_0')

%% let's generate the countor plot
figure(5)

contourf(x_r/lamz,y_r/lamz,Vt_r/A0,'LineColor','none')
colorbar
title('The displacement due to tumor in y direction')
xlabel('x/\lambda_z')
ylabel('y/\lambda_z')
hcb=colorbar
title(hcb,'u_{y}/A_0')

% tumor disturbance 
figure(6)
plot(x_mother,f)