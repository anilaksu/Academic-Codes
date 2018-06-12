function [ f ] = getTumorAgnesi( A,sig,x_0,x_int,x_trans,kx,kz)
% this function returns the shape of the tumor
theta=atan(kx/kz);
% back rotated x coordinate
x=x_trans(1)*cos(theta)-x_trans(2)*sin(theta);

% function itself
f=-2.*A*(x-(x_0-x_int(1)))/((1+(sig*(x-(x_0-x_int(1))))^2)^2);

end


