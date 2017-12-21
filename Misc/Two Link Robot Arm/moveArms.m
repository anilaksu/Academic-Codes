function [ x_1,x_2,y_1,y_2 ] = moveArms( d1,d2,theta_0,alpha_0,theta_f,alpha_f,N )
% this function calculates the movement of the arm
theta=linspace(theta_0,theta_f,N);
alpha=linspace(alpha_0,alpha_f,N);

% the end location of the first arm
x_1=d1*cosd(theta);
y_1=d1*sind(theta);
% the end location of the second arm
x_2 = d1 * cos(theta) + d2 * cos(theta + alpha); 
y_2 = d1 * sin(theta) + d2 * sin(theta+alpha); 

end

