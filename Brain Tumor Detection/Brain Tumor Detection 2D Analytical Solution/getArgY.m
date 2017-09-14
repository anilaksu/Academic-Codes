function [ k_y ] = getArgY( n,H )
%this function returns nth argument in y direction
% n is the order of the argument in y direction 
% y is the coordinate
% H  is the domain length in y direction 
k_y=(pi/H +2*(n-1)*pi/H);

end

