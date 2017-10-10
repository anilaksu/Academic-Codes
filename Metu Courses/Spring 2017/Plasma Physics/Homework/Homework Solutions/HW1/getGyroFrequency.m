function  [ omega_ce ] = getGyroFrequency(e,B,m,c )
%this function calculates the gyro frequency under given parameters
omega_ce=e*B/(m*c);

end

