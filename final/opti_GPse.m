clear all;
%load('data_for_opti.mat') %Data from a simulation with no GP
load('Data_Fault_GP.mat')
Z=Z(:,[true,true,true,true,false,true]);

lb = [0, 0, 0];
ub = [inf, inf, inf];
z_0_1 = [0.000001, 10, 10];

param_opti_1 = fmincon(@(z) opti_GPse_cost(z, Z, Y_1), z_0_1 , [], [], [], [], lb, ub); 
param_opti_2 = fmincon(@(z) opti_GPse_cost(z, Z, Y_2), z_0_1 , [], [], [], [], lb, ub);
param_opti_3 = fmincon(@(z) opti_GPse_cost(z, Z, Y_3), z_0_1 , [], [], [], [], lb, ub);
param_opti_4 = fmincon(@(z) opti_GPse_cost(z, Z, Y_4), z_0_1 , [], [], [], [], lb, ub);

noivar = [param_opti_1(1), param_opti_2(1), param_opti_3(1), param_opti_4(1)];
sig_f = [param_opti_1(2), param_opti_2(2), param_opti_3(2), param_opti_4(2)];
sig_l = [param_opti_1(3), param_opti_2(3), param_opti_3(3), param_opti_4(3)];





















