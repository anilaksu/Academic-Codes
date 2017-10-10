function [ x ] = getClosedFormShock( alpha,v,v_0,v_1 )
%this function computes the analytical solution for 1-D shock relation

% the first argument inside the logarithm
arg1=(v_0-v)^(v_0/(v_0-v_1));
% the first argument inside the logarithm
arg2=(v-v_1)^(-1.*v_1/(v_0-v_1));
x= alpha*log(arg1*arg2);

end

