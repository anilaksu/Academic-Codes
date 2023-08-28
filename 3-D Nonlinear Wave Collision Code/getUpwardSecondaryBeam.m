function [ u, p, rho ] = getUpwardSecondaryBeam( F_coll, x_coll, lamx, Nx, Ny )
% this function computes upward propagating secondary wave beam 
% Nx, Ny and Nz are number of grid points in x, y and z directions
% respectively
for i=1:Nx
    for j=1:Ny
     %x coordinate   
     x(i,j)=2.*(i-1)*lamx/Nx-lamx;     
     %y coordinate  
     y(i,j)=2.*(j-1)*lamx/Ny-lamx; 
     % the velocity field
     u(i,j)= double(F_coll(x(i,j)+x_coll(1),y(i,j)+x_coll(2), x_coll(3)));
      % the pressure field
     p(i,j)= double(F_coll(x(i,j)+x_coll(1),y(i,j)+x_coll(2), x_coll(3)));
      % the density field
     rho(i,j)= double(F_coll(x(i,j)+x_coll(1),y(i,j)+x_coll(2), x_coll(3)));
    end
end

end

