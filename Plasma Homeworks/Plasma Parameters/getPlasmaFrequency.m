function  [ omega_pe ] = getPlasmaFrequency( eps,n,e,m )
%this function calculates the plasma frequency under given parameters
omega_pe=sqrt(4.*pi*n*(e^2)/m);

end

