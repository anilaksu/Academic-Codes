clear all 
format long


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%          Excitation Parameters          %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the freqency of the excitation
omega=0.176;
% the wave length in z direction
lamz=0.875/2.;
% the wave number in z direction
kz=2*pi/lamz;
% the distance between the excitation regions
D=8.*lamz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%             Grid Parameters             %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the number of points in x direction
Nx=200;
% the number of points in y direction
Ny=200;

% the length of the domain
L=20*lamz;
H=20*lamz;

% the mother interval 
x_mother=linspace(0,L,Nx);
y_mother=linspace(-H/2,H/2,Ny);


%% the grid
x_grid=zeros(Nx*Ny,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%        Wave Energetics Parameters       %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the maximum amplitude of the excitation
A0=0.3;

for A=0.1:0.01:0.3
    for sig=2:0.1:20        
        for x_0=L/4:L/8:3*L/4
            [u_x,res,x_grid] = mainTumor( A,A0,sig,lamz,D,L,H,omega,Nx,Ny,x_0);
            % convert to a single vector by concat
            res2 = reshape(res,[Nx*Ny*2,1]);
%             % normalize matrix within itself before writing file
            res2Norm = (res2 - min(res2))/(max(res2) - min(res2));
            fileName= strcat('D:\braintumor_data_agnesi\',num2str(A),'_',num2str(sig),'_',num2str(x_0),'.txt')
            if(max(res2)>100)
                error('cortladý')
            end
            csvwrite(fileName,res2Norm);
        end
    end
end

% kanka banu taþ gibi 



