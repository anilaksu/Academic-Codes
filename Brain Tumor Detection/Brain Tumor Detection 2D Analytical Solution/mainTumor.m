function [u_x,res,x_grid] = mainTumor( A,A0,sig,lamz,D,L,H,omega,Nx,Ny,x_0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%        Wave Energetics Parameters       %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the maximum amplitude of the excitation
%A0
% the half-width of the wave envelope
%sig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Excitation Parameters          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the freqency of the excitation
%omega
% the wave length in z direction
%lamz
% the wave number in z direction
kz=2*pi/lamz;
% the tumor amplitude
%A
% the distance between the excitation regions
%D

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
%A0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Elasticity Parameters          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lame constant 
lambda=40.38*10^-5;

% let's generate the grid
for i=1:Nx
    for j=1:Ny
        % the counter
        count=(i-1)*Ny+j;
        x_grid(count,1)=x_mother(i);
        x_grid(count,2)=y_mother(j);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%     Generated Displacement Field        %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ u_x, k_ny, k_nx, alpha_n ] = PrimaryDispField( omega,lamz,D,A0,L,H,Nx,Ny );
% tumor itself
[ f ] = getTumorAgnesi( A,sig,x_0,L,Nx );
% [ f ] = getTumor( A,sig,x_0,L,Nx );
% disturbance due to tumor
[ u_xt ] = getTumorDisturbace( A,f,u_x,L,H,omega,Nx,Ny);
%% let's reshape vectors to plot them
% U_r = reshape(u_x(:,1),Ny,Nx);
% V_r = reshape(u_x(:,2),Ny,Nx);
% Ut_r = reshape(u_xt(:,1),Ny,Nx);
% Vt_r = reshape(u_xt(:,2),Ny,Nx);
% x_r = reshape(x_grid(:,1),Ny,Nx);
% y_r = reshape(x_grid(:,2),Ny,Nx);


res = (u_x+u_xt)/A0;
end

