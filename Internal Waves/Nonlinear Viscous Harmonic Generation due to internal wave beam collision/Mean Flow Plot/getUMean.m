function [ u_mean ] = getUMean( dv_meandz,Nx,L )
%this function generates mean velocity profile in x direction
%by satisfying the continuity equation
u_mean=zeros(Nx,1);
% the step size
dx=L/Nx
% let's generate it 
for i=2,Nx
    u_mean(i,1)= u_mean(i-1,1)+0.5*(dv_meandz(i,1)+ dv_meandz(i-1,1))*dx;
end

end

