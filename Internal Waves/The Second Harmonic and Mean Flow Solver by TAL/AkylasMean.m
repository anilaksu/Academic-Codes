clear all 
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%       Dispersion Relation Parameters    %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the BV frequency 
N=2.35;
% the freqency of the internal wave beam 
omega=0.5*N;
% the wave length in x direction
lamx=0.875;
% the wave number in x direction
kx_0=2*pi/lamx;
% the half width 
sig=0.5;

% the angle 
theta=asin(omega/N);
% alpha 
alpha=sin(theta)/sig;
% alpha 
beta=cos(theta)/sig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Mean Flow Compuations          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% syms x z a kx kpx x_2
% F_x = exp(-(x*alpha)^2);
% F_z = exp(-(z*beta)^2);
%  % the stream function
% % Q_0=4*sin(2*cot(theta)*kx_0*z)*F_x*F_z/cos(theta);
% % Q0=symfun(Q_0,[x,z]);
% % U=real(diff(Q0,z));
% % W=real(diff(Q0,x));
% w_0=2*kx_0*cos(2*cot(theta)*kx_0*z)*F_x*F_z/cos(theta);
% W=symfun(w_0,[x,z]);
% the number of points in x direction
Nx=200;
% the number of points in y direction
Ny=200;
% the length of the domain
L=20*lamx;
H=5*lamx;

% the mother interval 
x_mother=linspace(-L/2,L/2,Nx);
y_mother=linspace(0,-H,Ny);

% let's generate the grid
for i=1:Nx
    for j=1:Ny
        % the counter
        count=(i-1)*Ny+j;
        x_grid(count,1)=x_mother(i);
        x_grid(count,2)=y_mother(j);
    end
end
u_mean=zeros(Nx*Ny,1);
% let's perform the coordinate transformation and velocity calculation
for i=1:Nx*Ny
    % to filter noise from plot
    %if(abs(U( x_grid(i,1),x_grid(i,2)))/N>0.01)
       % u_mean(i,1) = real( U( x_grid(i,1),x_grid(i,2)))/N;
       % w_mean(i,1) = real( W( x_grid(i,1),x_grid(i,2)))/N;
    %end
    w_mean(i,1)=2*kx_0*cos(2*cot(theta)*kx_0*x_grid(i,2))*exp(-(x_grid(i,1)*alpha)^2)*exp(-(x_grid(i,2)*beta)^2)/(N*cos(theta));
end

%% let's reshape vectors to plot them
%U_r = reshape(u_mean(:,1),Ny,Nx);
W_r = reshape(w_mean(:,1),Ny,Nx);

x_r = reshape(x_grid(:,1),Ny,Nx);
y_r = reshape(x_grid(:,2),Ny,Nx);



% %% let's generate the countor plot
% figure(1)
% 
% contourf((x_r+L/2)/lamx,(y_r+H)/lamx,U_r,'LineColor','none')
% colorbar
% colorbar
% title('The mean velocity in x direction')
% xlabel('x/\lambda_x')
% ylabel('z/\lambda_x')
% hcb=colorbar
% title(hcb,'u_{mean}/A_{0}^2')
% pbaspect([2 1 1])

%% let's generate the countor plot
figure(2)

contourf((x_r+L/2)/lamx,(y_r+H)/lamx,W_r,'LineColor','none')
colorbar
colorbar
title('The mean velocity in z direction by TAL')
xlabel('x/\lambda_x')
ylabel('z/\lambda_x')
hcb=colorbar
title(hcb,'w_{mean}/A_{0}^2')
pbaspect([2 1 1])



