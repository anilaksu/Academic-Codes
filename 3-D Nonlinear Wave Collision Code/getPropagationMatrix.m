function [ D_xy ] = getPropagationMatrix( Cg, Lx, Ly, Lz, Nx, Ny, Nz )
% this function generates propagation matrix for the secondary
% inertia-gravity wave beams
% Cg: Group velocity vector
% Lx, Ly, Lz: Domain length in x, y, z directions respectively
% Nx, Ny, Nz: Number of grid points in x, y, z directions respectively

Dx = zeros(Nx,Nx);         % differentiation matrix in x direction
Dy = zeros(Ny,Ny);         % differentiation matrix in x direction
Dz = zeros(Nz,Nz);         % differentiation matrix in x direction

for i=1:(Nx-1)
   if(sign(Cg(1)) == 1) % if beam is leftward propagating
       Dx(1,1) = 1;        % No flux boundary condition on the right boundary
       Dx(i+1,i) = -Cg(1)*(Nx-1)/Lx; % First derivative with first order scheme in x direction
       Dx(i+1,i+1) = Cg(1)*(Nx-1)/Lx;  % First derivative with first order scheme in x direction
   else
       Dx(Nx,Nx) = 1;        % No flux boundary condition on the right boundary
       Dx(i,i) = -Cg(1)*(Nx-1)/Lx; % First derivative with first order scheme in x direction
       Dx(i,i+1) = Cg(1)*(Nx-1)/Lx;% First derivative with first order scheme in x direction
   end
end

for i=1:(Ny-1)
   if(sign(Cg(2)) == 1) % if beam is onward propagating
       Dy(1,1) = 1;        % No flux boundary condition on the front boundary
       
       if (i == (Ny-1))                    % Last node backward difference
           Dy(i+1,i) = -Cg(2)*(Ny-1)/Ly;   % First derivative with first order scheme in x direction
           Dy(i+1,i+1) = Cg(2)*(Ny-1)/Ly;  % First derivative with first order scheme in x direction
       else
           Dy(i+1,i) = -Cg(2)*(Ny-1)/Ly; % First derivative with first order scheme in x direction
           Dy(i+1,i+1) = Cg(2)*(Ny-1)/Ly;  % First derivative with first order scheme in x direction
       end
   else
       Dy(Ny,Ny) = 1;        % No flux boundary condition on the right boundary
       Dy(i,i) = -Cg(2)*(Ny-1)/Ly; % First derivative with first order scheme in x direction
       Dy(i,i+1) = Cg(2)*(Ny-1)/Ly;% First derivative with first order scheme in x direction
   end
end



D_xy = zeros(Nx*Ny,Nx*Ny);   % System Matrix 


for j = 1:Ny
    i_lstart = (j-2)*Nx+1; % Starting index on left block of matrix
    i_lend = (j-1)*Nx;     % end index on left block of matrix
    
    i_start = (j-1)*Nx+1; % Starting index
    i_end = j*Nx;         % end index
    
    i_rstart = j*Nx+1;    % Starting index on right block of matrix
    i_rend = (j+1)*Nx;    % end index on right block of matrix
    
    if (j == 1)
        D_xy(i_start:i_end,i_start:i_end) = Dx + eye(Nx)*Dy(j,j);        % Let's fill with system matrix
        D_xy(i_start:i_end,i_rstart:i_rend) =  eye(Nx)*Dy(j,j+1);      % Let's fill with system matrix
    elseif( j == Ny)
        D_xy(i_start:i_end,i_start:i_end) = Dx + eye(Nx)*Dy(j,j);        % Let's fill with system matrix
        D_xy(i_start:i_end,i_lstart:i_lend) =  eye(Nx)*Dy(j,j-1);      % Let's fill with system matrix
        j
    else
        D_xy(i_start:i_end,i_start:i_end) = Dx + eye(Nx)*Dy(j,j);        % Let's fill with system matrix
        D_xy(i_start:i_end,i_rstart:i_rend) = eye(Nx)*Dy(j,j+1);       % Let's fill with system matrix
        D_xy(i_start:i_end,i_lstart:i_lend) = eye(Nx)*Dy(j,j-1);       % Let's fill with system matrix
        j
    end
            
end


end

