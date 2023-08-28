function [ beta,phi_f ] = getNonlinearForcingBetas( gamma_1,gamma_2,phase_1,phase_2,k_1,k_2 )
% this function computes beta coefficients and the phase shifts of the
% nonlinear forcing term 

% the sub beta terms in the nonlinear forcing term in momentum equations
for j=1:3
    % the first three terms
    for i=1:3
        % the beta terms in each nonlinear forcing 
        beta_vel(i,j)=gamma_2(i)*gamma_1(j)*k_1(i);
    end
     for i=1:3
        % the beta terms in each nonlinear forcing 
        beta_vel(i+3,j)=gamma_1(i)*gamma_2(j)*k_2(i);
     end 
     % let's compute the corresponding lambda and chi coefficients 
    lambda_s(j)=0.;
    chi_s(j)=0.;
    for i=1:3
        lambda_s(j)=lambda_s(j)+beta_vel(i,j)*sin(phase_1(j)+phase_2(i))+beta_vel(i+3,j)*sin(phase_2(j)+phase_1(i));
        chi_s(j)=chi_s(j)+beta_vel(i,j)*cos(phase_1(j)+phase_2(i))+beta_vel(i+3,j)*cos(phase_2(j)+phase_1(i));
    end
    % the beta coefficient and the phase difference of the nonlinear forcing
    beta(j)=sqrt(lambda_s(j)^2+chi_s(j)^2);
    phi_f(j)=atan(chi_s(j)/lambda_s(j));
end

% the sub beta terms in the nonlinear forcing term in energy equation
% the first three terms
for i=1:3
    % the beta terms in each nonlinear forcing 
    beta_vel(i,4)=gamma_2(i)*gamma_1(5)*k_1(i);
end
for i=1:3
    % the beta terms in each nonlinear forcing 
    beta_vel(i+3,4)=gamma_1(i)*gamma_2(5)*k_2(i);
end 
% let's compute the corresponding lambda and chi coefficients 
lambda_s(4)=beta_vel(1,4)*sin(phase_2(1))+beta_vel(2,4)*sin(phase_2(2))+beta_vel(4,4)*sin(phase_1(1))+beta_vel(5,4)*sin(phase_1(2));
chi_s(4)=0.;
for i=1:3
    chi_s(4)=chi_s(4)+beta_vel(i,4)*cos(phase_2(i))+beta_vel(i+3,4)*cos(phase_1(i));
end
% the beta coefficient and the phase difference of the nonlinear forcing
beta(4)=sqrt(lambda_s(4)^2+chi_s(4)^2);
phi_f(4)=atan(chi_s(4)/lambda_s(4));

end

