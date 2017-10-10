function [ v ] = getShockVelocity( alpha,v_0,v_1,x )
%this function calculates the velocity at a given location

v_left=v_0;
v_right=v_1;
v_mid=0.5*(v_0+v_1);

% the first error
x_it = getClosedFormShock( alpha,v_mid,v_0,v_1 );
err=abs(x-x_it);
while(err>10^-2.)
    % the error on right point
    [ x_right ] = getClosedFormShock( alpha,v_right,v_0,v_1 );
    f_right=x-x_right;
    % the error on mid point
    [ x_mid ] = getClosedFormShock( alpha,v_mid,v_0,v_1 );
    f_mid=x-x_mid;
    if(f_mid*f_right<0)
        v_left=v_mid;
    else 
        v_right=v_mid;
    end 
    % let's update mid velocity
    v_mid=0.5*(v_left+v_right);
    % let's calculate the error
    x_it= getClosedFormShock( alpha,v_mid,v_0,v_1 );
    err=abs(x-x_it);
end
% the final update
v=v_mid;

end

