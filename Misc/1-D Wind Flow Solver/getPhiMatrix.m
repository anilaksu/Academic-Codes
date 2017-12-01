function [ Mp ] = getPhiMatrix( Nx,dx,Ly)
% this function generates phi matrix that solves the second equation
Mp=toeplitz([2 -1 zeros(1, Nx-3)])/(dx^2)-eye(Nx-1)*(pi/Ly)^2;
% it has to be inverted to reduce the computational effort for the rest of
% the computation
Mp=inv(Mp);

end

