%% this script is developed to solve 1-D wind driven flow with only
%% bottom friction developed by Bulut Çaðdaþ

% the domain length 
Lx=1000*10^3;
Ly=1000*10^3;

% other supplementary parameters
r=10^-6;
D=1000.;
beta=2.*10^-11;
tau=0.1;
rho=1000.;
% the step size in x direction
dx=5*10^3;
% the time step size
dt=0.001;
% the number of intervals in x direction
Nx=round(Lx/dx);
% the grid
x=linspace(0,Lx,Nx+1);
% the number of time steps
Nt=100;
% let's generate the solver matrix for phi
[ Mp ] = getPhiMatrix( Nx,dx,Ly);
% the first derivative matrix
[ Diffx ] = getFirstDervMatrix( Nx,dx );

% Z variable
Z=zeros(Nx-1,Nt);
% phi variable
phi=zeros(Nx-1,Nt);
% the time array
time=zeros(Nt,1);
% time integration
for i=2:Nt
   % the time array
   time(i,1)=i*dt;
   % the time integration of Z
   Z(:,i)=dt*((1-r)*Z(:,i-1)-beta*Diffx*phi(:,i-1)-ones(Nx-1,1)*(tau*pi)/(rho*D*Ly));
   %-ones(Nx-1,1)*(tau*pi)/(rho*D*Ly);
   % the time update of phi  
   phi(:,i)=Mp*Z(:,i);
end

figure(1)

plot(x(2:Nx),Z(:,Nt),'r*','LineWidth',2)
hold on
plot(x(2:Nx),Z(:,Nt/2),'b','LineWidth',2)
hold on
plot(x(2:Nx),Z(:,5),'mo','LineWidth',2)

title('Distribution of Z')
xlabel('x (km)')
ylabel('Z')
legend('t=1','t=5','t=10')

figure(2)

plot(x(2:Nx),phi(:,Nt),'r*','LineWidth',2)
hold on
plot(x(2:Nx),phi(:,Nt/2),'b','LineWidth',2)
hold on
plot(x(2:Nx),phi(:,5),'mo','LineWidth',2)

title('Distribution of Z')
xlabel('x (km)')
ylabel('Z')
legend('t=1','t=5','t=10')