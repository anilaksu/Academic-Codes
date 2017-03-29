function [ x_r ] = getRayPath(x_start,x_end,N )
%this function generates the ray path for given group velocity profile
% Cg is the group velocity vector 
% x_start is the start point obviously
% N is the number of the points
% dx is the step size

% it is distributed linearly
x_r(:,1)=linspace(x_start(1),x_end(1),N);
x_r(:,2)=linspace(x_start(2),x_end(2),N);



end

