function [ theta_1, theta_2 ] = inverseKinematices( d1,d2,x_1,y_1 )
% this function calculates the angle values using inverse kinematics
% method
%x_1 = d1 * cos(THETA1) + d2 * cos(THETA1 + THETA2); 
%y_1 = d1 * sin(THETA1) + d2 * sin(THETA1 + THETA2); 


% the second arms angle
theta_2=acos((x_1^2+y_1^2-d1^2-d2^2)/(2.*d1*d2))
% the first arms angle
theta_1=acot((-d2*sin(theta_2)*x_1+(d1+d2*cos(theta_2))*y_1)/(d2*sin(theta_2)*y_1+(d1+d2*cos(theta_2))*x_1))
end

