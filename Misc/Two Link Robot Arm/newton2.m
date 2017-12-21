function [ thetaF,alphaF,Angles ] = newton2( d1,d2,theta0,alpha0,xfinal,yfinal,N  )
%this function performs Newton Raphson method
Angles=zeros(2,N);
% Vector
G=zeros(2,N);
% Jacobian matrice
J=zeros(2,2);
% the inital angles
Angles(1,1)=theta0;
Angles(2,1)=alpha0;

for i=1:(N-1)
   % the end location of the second arm
    G(1,i)=xfinal - d1 * cosd(Angles(1,i)) - d2 * cosd(Angles(1,i) + Angles(2,i)); 
    G(2,i)=yfinal - d1 * sind(Angles(1,i)) - d2 * sind(Angles(1,i)+Angles(2,i)); 
    % the jacobian matrice
    J(1,1)=d1 * sind(Angles(1,i)) + d2 * sind(Angles(1,i) + Angles(2,i)); 
    J(1,2)= d2 * sind(Angles(1,i) + Angles(2,i)); 
    J(2,1)=- d1 * cosd(Angles(1,i)) - d2 * cosd(Angles(1,i)+Angles(2,i));
    J(2,2)= - d2 * cosd(Angles(1,i)+Angles(2,i));
    % the angles at next iteration
    Angles(:,i+1)=Angles(:,i)-inv(J)*G(:,i);
end

thetaF=Angles(1,N);
alphaF=Angles(2,N);
end

