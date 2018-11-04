function [ f ] = getIntFunc( phid, k, a, phii, phio, r)
% this function calculates the function inside the asymptotic integral

f=cos(phid-phii)*exp(-i*k*a*cos(phid-phii))*exp(-i*k*sqrt((r*cos(phio)-a*cos(phid))^2+(r*sin(phio)-a*sin(phid))^2));

end

