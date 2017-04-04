function [ x_ref ] = getReflectionPoint( x_start,kx,kz,H )
%finds the reflection point for given start point x_start
%where H is the domain height 
theta=atan(kz/kx);
% the height is alrady known 
x_ref(2)=H;
%but delta x is not it is calculate as 
x_ref(1)=x_start(1)+(H-x_start(2))*tan(theta);

end

