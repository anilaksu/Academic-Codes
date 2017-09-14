function [ f ] = getTumorAgnesi( A,sig,x_0,L,Nx )
% this function returns the shape of the tumor

% the mother interval 
x=linspace(0,L,Nx);
    % function itself
    for i=1:Nx
        f(i,1)=-2.*A*(x(i)-x_0)/((1+(sig*(x(i)-x_0))^2)^2);
    end

end


