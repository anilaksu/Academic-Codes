function [ l_d ] = getDebyeLength( eps,T,n,e )
%this function calculates the debye length under given parameters
l_d=sqrt(eps*T/(n*(e^2)));

end

