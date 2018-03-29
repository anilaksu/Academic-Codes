% this routine is generated to visualize 2D solution
 clear all
 format long
SimData1=importdata('E_tot.dat');
x=SimData1(:,1);
E_field=SimData1(:,2);

figure(1)

plot(x,E_field)
title('The electric field in z direction')
xlabel('z/\lambda_z')
ylabel('E_{field}/A_0')


