function [ x_trans ] = transformCoordinate( x,x_ref,Cg )
%this function performs an orthogonal coordinate transformation
% x_ref is the reference point where the coordinate system is rotated
% the group velocity vector is used as reference

% the rotation angle 
theta=atan(Cg(2)/Cg(1));

% the rotation matrix
A=zeros(2,2);

%% the first row of the rotation matrix
A(1,1)=cos(theta);
A(1,2)=sin(theta);
%% the second row of the rotation matrix
A(2,1)=-sin(theta);
A(2,2)=cos(theta);

%% the position vector with respect to x_ref
delx=x-x_ref;

%% since the matrix is small let's perfom the matrix vector multiplication explicitly
x_trans(1)=A(1,1)*delx(1)+A(1,2)*delx(2);
x_trans(2)=A(2,1)*delx(1)+A(2,2)*delx(2);

end

