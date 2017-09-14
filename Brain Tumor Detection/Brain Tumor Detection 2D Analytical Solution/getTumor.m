function [ f ] = getTumor( A,sig,x_0,L,Nx )
% this function returns the shape of the tumor

% the mother interval 
x=linspace(0,L,Nx);

% function itself
for i=1:Nx
    f(i,1)=A*(x(i)-x_0)*exp(0.5*((x(i)-x_0)/(sig))^2)/sig^2;
   %  f(i,1)=A*exp(0.5*((x(i)-x_0)/(sig))^2);
end

end

