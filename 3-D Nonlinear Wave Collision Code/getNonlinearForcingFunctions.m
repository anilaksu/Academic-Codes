function [ omega_app, sig_up, sig_vp, sig_wp, sig_pp, sig_rp, kz_p ] = getNonlinearForcingFunctions(beta,phi_f, N, f, omega_0, kx_0, ky_0 )
% this function computes the sigmas of the nonlinear forcing functions

% the gravity 
g=9.81;
% the density 
rho_0=1000.;

syms kx ky kz  

% the wavenumber in z direction according to dispersion relation !!!!!!
kz_0 = sqrt((N^2-omega_0^2)*(kx_0^2+ky_0^2)/(omega_0^2-f^2))
% the wavenumbers in z direction
kz_p = sqrt((N^2-omega_0^2)*(kx^2+ky^2)/(omega_0^2-f^2));
kz_p = symfun(kz_p,[kx ky]);
kz_m = -1.*sqrt((N^2-omega_0^2)*(kx^2+ky^2)/(omega_0^2-f^2));
kz_m = symfun(kz_m,[kx ky]);
% the derivatives of kz wrt to kx and ky
kz_x = symfun(diff(kz_p,kx),[kx ky]);
kz_y = symfun(diff(kz_p,ky),[kx ky]);
% the first order taylor approximation of kz
kz_p = kz_0+ kz_x(kx_0, ky_0)*(kx-kx_0)+kz_y(kx_0, ky_0)*(ky-ky_0);
kz_p = symfun(kz_p,[kx ky kz]);
%kz_p = matlabFunction(kz_p);

% the frequency
omega = sqrt(((N^2)*(kx^2+ky^2)+(f^2)*(kz^2))/(kx^2+ky^2+kz^2));
omega = symfun(omega,[kx ky kz]);

% the  group velocity
Cg_x = symfun(diff(omega,kx),[kx ky kz]);
Cg_y = symfun(diff(omega,ky),[kx ky kz]);
Cg_z = symfun(diff(omega,kz),[kx ky kz]);

% the expansion of frequency around given wavenumber vector
omega_app=Cg_x(kx_0, ky_0,kz_0)*(kx-kx_0)+Cg_y(kx_0, ky_0,kz_0)*(ky-ky_0)+Cg_z(kx_0, ky_0,kz_0)*(kz-kz_0);
omega_app= symfun(omega_app,[kx ky kz]);

% the sigma coefficients for forcing in u velocity
sig_uu=symfun(-i*beta(1)*exp(-i*phi_f(1))*((N*ky)^2-(ky^2+kz_p^2)*omega_0^2)/omega_0,[kx ky kz]);
sig_uv=symfun(i*beta(2)*exp(-i*phi_f(2))*(kx*ky*(N^2-omega_0^2)+i*kx^2*omega_0*f)/omega_0,[kx ky kz]);
sig_uw=symfun(i*beta(3)*exp(-i*phi_f(3))*kz*(kx*omega_0+i*ky*f),[kx ky kz]);
sig_ur=symfun(beta(4)*exp(-i*phi_f(4))*g*kz*(ky*f-i*kx*omega_0)/(rho_0*omega_0),[kx ky kz]);
% the total forcing coefficient
sig_up=symfun((sig_uu+sig_uv+sig_uw+sig_ur)/(2.*(kx_0^2+ky_0^2+kz_0^2)*abs(Cg_z(kx_0, ky_0,kz_0))),[kx ky kz]);
sig_up = matlabFunction(sig_up);

% the sigma coefficients for forcing in v velocity
sig_vu=symfun(-i*beta(1)*exp(-i*phi_f(1))*(kx*ky*(N^2-omega_0^2)-i*kx^2*omega_0*f)/omega_0,[kx ky kz]);
sig_vv=symfun(-i*beta(2)*exp(-i*phi_f(2))*(kx^2*(N^2-omega_0^2)-i*kz^2*omega_0^2)/omega_0,[kx ky kz]);
sig_vw=symfun(-i*beta(3)*exp(-i*phi_f(3))*kx*kz*(omega_0-i*f),[kx ky kz]);
sig_vr=symfun(-beta(4)*exp(-i*phi_f(4))*g*kz*(kx*f+i*ky*omega_0)/(rho_0*omega_0),[kx ky kz]);
% the total forcing coefficient
sig_vp=symfun((sig_vu+sig_vv+sig_vw+sig_vr)/(2.*(kx_0^2+ky_0^2+kz_0^2)*abs(Cg_z(kx_0, ky_0,kz_0))),[kx ky kz]);
sig_vp = matlabFunction(sig_vp);

% the sigma coefficients for forcing in w velocity
sig_wu=symfun(-beta(1)*exp(-i*phi_f(1))*kz*(i*kx*omega_0+ky*f),[kx ky kz]);
sig_wv=symfun(beta(2)*exp(-i*phi_f(2))*kz*(-i*kx*omega_0+ky*f),[kx ky kz]);
sig_ww=symfun(i*beta(3)*exp(-i*phi_f(3))*omega_0*(kx^2+ky^2),[kx ky kz]);
sig_wr=symfun(i*beta(4)*exp(-i*phi_f(4))*g*(kx^2+ky^2)/(rho_0),[kx ky kz]);
% the total forcing coefficient
sig_wp=symfun((sig_wu+sig_wv+sig_ww+sig_wr)/(2.*(kx_0^2+ky_0^2+kz_0^2)*abs(Cg_z(kx_0, ky_0,kz_0))),[kx ky kz]);
sig_wp = matlabFunction(sig_wp);

% the sigma coefficients for forcing in pressure
sig_pu=symfun(beta(1)*exp(-i*phi_f(1))*(N^2-omega_0^2)*(i*kx*omega_0+ky*f)/omega_0,[kx ky kz]);
sig_pv=symfun(beta(2)*exp(-i*phi_f(2))*(N^2-omega_0^2)*(-i*ky*omega_0+kx*f)/omega_0,[kx ky kz]);
sig_pw=symfun(i*beta(3)*exp(-i*phi_f(3))*kz*(f^2-omega_0^2),[kx ky kz]);
sig_pr=symfun(i*beta(4)*exp(-i*phi_f(4))*g*kz*(f^2-omega_0^2)/(rho_0*omega_0),[kx ky kz]);
% the total forcing coefficient
sig_pp=symfun(rho_0*(sig_pu+sig_pv+sig_pw+sig_pr)/(2.*(kx_0^2+ky_0^2+kz_0^2)*abs(Cg_z(kx_0, ky_0,kz_0))),[kx ky kz]);
sig_pp = matlabFunction(sig_pp);

% the sigma coefficients for forcing in density
sig_ru=symfun(-i*beta(1)*exp(-i*phi_f(1))*N^2*kz*(i*kx*omega_0+ky*f)/(g*omega_0),[kx ky kz]);
sig_rv=symfun(i*beta(2)*exp(-i*phi_f(2))*i*N^2*kz*(-i*ky*omega_0+kx*f)/(g*omega_0),[kx ky kz]);
sig_rw=symfun(-beta(3)*exp(-i*phi_f(3))*N^2*(kx^2+ky^2)/g,[kx ky kz]);
sig_rr=symfun(-beta(4)*exp(-i*phi_f(4))*((kx^2+ky^2+kz^2)*omega_0^2-kz^2*f^2)/(rho_0*omega_0),[kx ky kz]);
% the total forcing coefficient
sig_rp=symfun(rho_0*(sig_ru+sig_rv+sig_rw+sig_rr)/(2.*(kx_0^2+ky_0^2+kz_0^2)*abs(Cg_z(kx_0, ky_0,kz_0))),[kx ky kz]);
sig_rp = matlabFunction(sig_rp);
end

