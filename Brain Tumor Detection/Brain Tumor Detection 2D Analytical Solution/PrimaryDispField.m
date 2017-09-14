function [ u_x, k_ny, k_nx, alpha_n ] = PrimaryDispField( omega,lamz,D,A,L,H,Nx,Ny )
%this function computes the primary displacement field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Excitation Parameters          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time: the time of the plot
% omega  the freqency of the excitation
% lamz: the wave length in z direction
% the excitaion amplitude
% D:the distance between the excitation regions
% the wave number in z direction
kz=2*pi/lamz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%             Grid Parameters             %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the number of points in x direction
% Nx
% the number of points in y direction
% Ny
% the length of the domain
%L
%H

% the mother interval 
x_mother=linspace(0,L,Nx);
y_mother=linspace(-H/2,H/2,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Elasticity Parameters          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lame constant 
lambda=40.38*10^-5;
% the displacement field in
u_x=zeros(Nx*Ny,2);


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

% modes
ny=12;
% the wave number in y direction 
for i=1:ny
    k_ny(i,1)=getArgY(i,H );
    % the wave number in x direction 
    k_nx(i,1)=(omega^2/lambda-k_ny(i,1)^2)^0.5;
    % fourier coefficients
    alpha_n(i,1)= getFourierCoeff( k_nx(i,1),k_ny(i,1),kz,H,D,A,L);
end

% let's perform the displacement field computation
for j=1:ny
    for i=1:Nx*Ny    
        % the argument of cosine in y direction
        phi_y=k_ny(j,1)*x_grid(i,2);
        % the argument in x direction 
        phi_x=k_nx(j,1)*x_grid(i,1)-k_nx(j,1)*L+pi/2;   
        % let's multiply amplitude with cos(phi_x)*sin(phi_y)
        u_x(i,1)=u_x(i,1)+alpha_n(j,1)*cos(phi_x)*sin(phi_y);
        u_x(i,2)=u_x(i,2)+(k_ny(j,1)/k_nx(j,1))*alpha_n(j,1)*sin(phi_x)*sin(phi_y);
    end
end 

end

