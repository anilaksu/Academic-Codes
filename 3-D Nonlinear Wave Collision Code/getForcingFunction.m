function [ F_coll ] = getForcingFunction( J_1,J_2,sig_1,sig_2,x_cen_1,x_cen_2 )
% this function generates the forcing function due to collision
% J_1: the transformation matrix for the first wave
% J_2: the transformation matrix for the second wave
% sig_1: the half width of the first wave beam
% sig_2: the half width of the second wave beam
% x_cen_1: the center of forcing of the first beam
% x_cen_2: the center of forcing of the second beam

syms x y z 

% transformed center of forcings 
x_trc_1= J_1*x_cen_1;
x_trc_2= J_2*x_cen_2;
% let's transform the coordinate
x_tr_1=J_1*[x-x_cen_1(1); y-x_cen_1(2); z-x_cen_1(3)];
eta_1=symfun(x_tr_1(2),[x y z]);
zeta_1=symfun(x_tr_1(3),[x y z]);

x_tr_2=J_2*[x-x_cen_2(1); y-x_cen_2(2); z-x_cen_2(3)];
eta_2=symfun(x_tr_2(2),[x y z])
zeta_2=symfun(x_tr_1(3),[x y z]);

% arg_1=((J_1(2,1)^2+J_1(3,1)^2)/(2.*sig_1^2)+(J_2(2,1)^2+J_2(3,1)^2)/(2.*sig_2^2))*x^2;
% arg_2=((J_1(2,2)^2+J_1(3,2)^2)/(2.*sig_1^2)+(J_2(2,2)^2+J_2(3,2)^2)/(2.*sig_2^2))*y^2;
% arg_3=((J_1(2,3)^2+J_1(3,3)^2)/(2.*sig_1^2)+(J_2(2,3)^2+J_2(3,3)^2)/(2.*sig_2^2))*z^2;
% arg_4=((J_1(2,1)*J_1(2,2)+J_1(3,1)*J_1(3,2))/(sig_1^2)+(J_2(2,1)*J_2(2,2)+J_2(3,1)*J_2(3,2))/(sig_2^2))*x*y;
% arg_5=((J_1(2,1)*J_1(2,3)+J_1(3,1)*J_1(3,3))/(sig_1^2)+(J_2(2,1)*J_2(2,3)+J_2(3,1)*J_2(3,3))/(sig_2^2))*x*z;
% arg_6=((J_1(2,2)*J_1(2,3)+J_1(3,2)*J_1(3,3))/(sig_1^2)+(J_2(2,2)*J_2(2,3)+J_2(3,2)*J_2(3,3))/(sig_2^2))*y*z;
% arg_7=(-(J_1(2,1)*x_trc_1(2)+J_1(3,1)*x_trc_1(3))/(sig_1^2)-(J_2(2,1)*x_trc_2(2)+J_2(3,1)*x_trc_2(3))/(sig_2^2))*x;
% arg_8=(-(x_trc_1(2)*J_1(2,3)+x_trc_1(3)*J_1(3,3))/(sig_1^2)+(x_trc_2(2)*J_2(2,3)+x_trc_2(3)*J_2(3,3))/(sig_2^2))*z;
% arg_9=(-(J_1(2,2)*x_trc_1(2)+J_1(3,2)*x_trc_1(3))/(sig_1^2)-(J_2(2,2)*x_trc_2(2)+J_2(3,2)*x_trc_2(3))/(sig_2^2))*y;
% arg_10=(x_trc_1(2)^2+x_trc_1(3)^2)/(2*sig_1^2)+(x_trc_2(2)^2+x_trc_2(3)^2)/(2.*sig_2^2);
% the exponential terms of the collision
F_coll=exp(-(eta_1^2+zeta_1^2)/(2.*sig_1^2)-(eta_2^2+zeta_2^2)/(2.*sig_2^2));
%F_coll=exp(-arg_1-arg_2-arg_3-arg_4-arg_5-arg_6-arg_7-arg_8-arg_9-arg_10);
F_coll=symfun(F_coll,[x y z]);
F_coll = matlabFunction(F_coll);

end

