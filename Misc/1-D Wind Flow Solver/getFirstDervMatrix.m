function [ Diffx ] = getFirstDervMatrix( Nx,dx )
% this function generates first derivative matrix

Diffx=-1.*eye(Nx-1)/dx;
for i=1:Nx-2
    Diffx(i,i+1)=1/dx;
end

end

