function [ u_xt ] = getTumorDisturbace(A,f,u_x,L,H,omega,Nx,Ny)
% this function calculates the disturbance due to tumor
%f: tumor disturbance
% u_x: displacement field

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
% modes
ny=12;
% alpha coefficients
alpha_n=zeros(ny,1);
for i=1:ny
    k_nx(i,1)=getArgY(i,H );
    % the wave number in x direction 
    k_ny(i,1)=(omega^2/lambda-k_nx(i,1)^2)^0.5;
    % fourier coefficients
    for j=1:Nx
        alpha_n(i,1)=  alpha_n(i,1)-2.*sin(k_nx(i,1)*x_mother(1,j))*f(j,1)*u_x(j,1)/L;
    end
end

% the displacement field in
u_xt=zeros(Nx*Ny,2);

% let's generate the grid
for i=1:Nx
    for j=1:Ny
        % the counter
        count=(i-1)*Ny+j;
        x_grid(count,1)=x_mother(i);
        x_grid(count,2)=y_mother(j);
    end
end

for j=1:ny
    for i=1:Nx*Ny    
        % the argument of cosine in y direction
        phi_y=k_ny(j,1)*x_grid(i,2);
        % the argument in x direction 
        phi_x=k_nx(j,1)*x_grid(i,1)-k_nx(j,1)*L+pi/2;   
        % let's multiply amplitude with cos(phi_x)*sin(phi_y)
        u_xt(i,1)=u_xt(i,1)+alpha_n(j,1)*cos(phi_x)*sin(phi_y);
        u_xt(i,2)=u_xt(i,2)+(k_ny(j,1)/k_nx(j,1))*alpha_n(j,1)*sin(phi_x)*sin(phi_y);
    end
end 
% let's multiply it with amplitude
u_xt=A*u_xt;

end

