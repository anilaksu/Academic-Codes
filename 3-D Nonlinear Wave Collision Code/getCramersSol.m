function [u_xp,u_yp,u_zp,p_p,rho_p, u_xm,u_ym,u_zm,p_m,rho_m] = getCramersSol( F_coll,beta,phi_f,N_1,g_0,f_0,omega_0, rho_0, kx_0,ky_0,kz_0 )
% this function computes the solution by Cramer's rule
% here rho: rho_0
% here g : g_0

syms x y z kx ky kz N f omega Nx Ny Nz Nr rho g 

% % the system matrix 
% A= [i*kx i*ky i*kz 0 0; -i*omega -f 0 i*kx/rho 0 ; f -i*omega 0  i*ky/rho 0; 0 0 -i*omega i*kz/rho g/rho; 0 0 -rho*N^2/g 0 -i*omega];
% % the nonlinear forcing in u 
% F_x= [0 i*ky i*kz 0 0; Nx -f 0 i*kx/rho 0 ; Ny -i*omega 0  i*ky/rho 0; Nz 0 -i*omega i*kz/rho g/rho; Nr 0 -rho*N^2/g 0 -i*omega];
% % the nonlinear forcing in v
% F_y= [i*kx 0 i*kz 0 0; -i*omega Nx 0 i*kx/rho 0 ; f Ny 0  i*ky/rho 0; 0 Nz -i*omega i*kz/rho g/rho; 0 Nr -rho*N^2/g 0 -i*omega];
% % the nonlinear forcing in w
% F_z= [i*kx i*ky 0 0 0; -i*omega -f Nx i*kx/rho 0 ; f -i*omega Ny  i*ky/rho 0; 0 0 Nz i*kz/rho g/rho; 0 0 Nr 0 -i*omega];
% % the nonlinear forcing in p
% F_p = [i*kx i*ky i*kz 0 0; -i*omega -f 0 Nx 0 ; f -i*omega 0  Ny 0; 0 0 -i*omega Nz g/rho; 0 0 -rho*N^2/g Nr -i*omega];
% %F_rho= [i*kx i*ky i*kz 0 0; -i*omega -f 0 i*kx/rho_0 0 ; f -i*omega 0  i*ky/rho_0 0; 0 0 -i*omega i*kz/rho_0 g/rho_0; 0 0 -rho_0*N^2/g 0 -i*omega];
% F_r= [i*kx i*ky i*kz 0 0; -i*omega -f 0 i*kx/rho Nx ; f -i*omega 0  i*ky/rho Ny; 0 0 -i*omega i*kz/rho Nz; 0 0 -rho*N^2/g 0 Nr];
% 
% M= symfun(A,[kx ky kz omega]);
% 
% F_r= symfun(det(F_r),[kx ky kz omega])
% 
% u= det(M);
% u_f =symfun(u,[kx ky kz omega]);

F_coll=symfun(F_coll*exp(i*(kx_0*x+ky_0*y+kz_0*z)),[x y z]);
F_2 = symfun(F_coll,[x y z]);
f_1four=symfun(fourier(F_2, z,kz),[x y kz]);
f_2four=symfun(fourier(f_1four, y,ky),[x ky kz]);
f_3four=symfun(fourier(f_2four, x,kx),[kx ky kz]);

[ omega, sig_u, sig_v, sig_w, sig_p, sig_r, kz_p ] = getNonlinearForcingFunctions(beta,phi_f, N_1, f_0, omega_0, kx_0, ky_0 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%    Upward Propagating Secondary Beam    %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the velocity in x direction
Envelope =symfun(sig_u*f_3four*dirac(kz+kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
u_xp= symfun(Env_xyz,[x y z]);
u_xp = matlabFunction(u_xp);

% the velocity in y direction
Envelope =symfun(sig_v*f_3four*dirac(kz+kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
u_yp= symfun(Env_xyz,[x y z]);
u_yp = matlabFunction(u_yp);

% the velocity in z direction
Envelope =symfun(sig_w*f_3four*dirac(kz+kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
u_zp= symfun(Env_xyz,[x y z]);
u_zp = matlabFunction(u_zp);

% the pressure
Envelope =symfun(sig_p*f_3four*dirac(kz+kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
p_p= symfun(Env_xyz,[x y z]);
p_p = matlabFunction(p_p);

% the density
Envelope =symfun(sig_r*f_3four*dirac(kz+kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
rho_p= symfun(Env_xyz,[x y z]);
rho_p = matlabFunction(rho_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%  Downward Propagating Secondary Beam    %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the velocity in x direction
Envelope =symfun(sig_u*f_3four*dirac(kz-kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
u_xm= symfun(Env_xyz,[x y z]);
u_xm = matlabFunction(u_xm);

% the velocity in y direction
Envelope =symfun(sig_v*f_3four*dirac(kz-kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
u_ym= symfun(Env_xyz,[x y z]);
u_ym = matlabFunction(u_ym);

% the velocity in z direction
Envelope =symfun(sig_w*f_3four*dirac(kz-kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
u_zm= symfun(Env_xyz,[x y z]);
u_zm = matlabFunction(u_zm);

% the pressure
Envelope =symfun(sig_p*f_3four*dirac(kz-kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
p_m= symfun(Env_xyz,[x y z]);
p_m = matlabFunction(p_m);

% the density
Envelope =symfun(sig_r*f_3four*dirac(kz-kz_p),[kx ky kz]);
% % % the integral of envelope
Env_z=symfun(ifourier(Envelope,kz,z),[kx ky z]);
Env_yz=symfun(ifourier(Env_z,ky,y),[kx y z]);
Env_xyz=symfun(ifourier(Env_yz,kx,x),[x y z]);
% % the integral of energy 
rho_m= symfun(Env_xyz,[x y z]);
rho_m = matlabFunction(rho_m);
end

