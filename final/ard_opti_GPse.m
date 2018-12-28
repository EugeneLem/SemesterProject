clear all;
%load('data_for_opti.mat') %Data from a simulation with no GP
load('Data_Fault_GP.mat')
Z=Z(:,[true,true,true,true,false,true]);


noisevar =    1.0e-07 *[0.499394606930038   0.000250000000000   0.331496659462403   0.004129682494018];
sigma_f = 1.0e+03 *[4.832024100700585   0.010013856207360   0.009999147113626   0.002790605058578];   %size: 1x4
sigma_l = 1.0e+04 *[1.199673163015586   0.001068857665558   0.000999885788767   0.001447417602870;...
       0.114051481875275   0.001000002438980   0.000999981881745   0.000863749522226;...
       0.114051481875275   0.001000002438980   0.000999981881745   0.000863749522226;...
       0.114051481875277   0.001000002438980   0.000999981881745   0.000863749522226;...
       0.114051481875277   0.001000002438980   0.000999981881745   0.000863749522226]; %size: 5x4


lb = [0, 0, 0, 0, 0, 0, 0];
ub = [inf, inf, inf, inf, inf, inf, inf];
z_0_1 = [1.0e-07, 9, 9, 9, 9, 9, 9];
z_0_2 = [1.0e-07, 1, 2, 2, 2, 2, 2 ];
z_0_3 = [1.0e-07, 1, 2, 2, 2, 2, 2 ];
z_0_4 = [1.0e-07, 1, 1, 1, 1, 1, 1];

param_opti_2 = fmincon(@(z) ard_opti_GPse_cost(z, Z, Y_2), z_0_2 , [], [], [], [], lb, ub);

param_opti_1 = fmincon(@(z) ard_opti_GPse_cost(z, Z, Y_1), z_0_1 , [], [], [], [], lb, ub); 
%param_opti_2 = fmincon(@(z) ard_opti_GPse_cost(z, Z, Y_2), z_0_2 , [], [], [], [], lb, ub);
param_opti_3 = fmincon(@(z) ard_opti_GPse_cost(z, Z, Y_3), z_0_3 , [], [], [], [], lb, ub);
param_opti_4 = fmincon(@(z) ard_opti_GPse_cost(z, Z, Y_4), z_0_4 , [], [], [], [], lb, ub);

noivar = [param_opti_1(1), param_opti_2(1), param_opti_3(1), param_opti_4(1)];
sig_f = [param_opti_1(2), param_opti_2(2), param_opti_3(2), param_opti_4(2)];
sig_l = [param_opti_1(3:end)', param_opti_2(3:end)', param_opti_3(3:end)', param_opti_4(3:end)'];





















