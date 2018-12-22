load('Data_noFault.mat') 
Z= Z(:,[true, true, true, true, false, true]);

L2NW_case = 2;
lb = [0, 0, 0];
ub = [inf , inf, inf];
z_0 =       1.0e+03 *[ 0.000033035090513   1.582306505888283   0.000422549219831];

param_opti_1 = fmincon(@(z) NW_cost(z, Z, Y_1, L2NW_case), z_0 , [], [], [], [], lb, ub); 

param_opti_2 = fmincon(@(z) NW_cost(z, Z, Y_2, L2NW_case), z_0 , [], [], [], [], lb, ub);

param_opti_3 = fmincon(@(z) NW_cost(z, Z, Y_3, L2NW_case), z_0 , [], [], [], [], lb, ub);

param_opti_4 = fmincon(@(z) NW_cost(z, Z, Y_4, L2NW_case), z_0 , [], [], [], [], lb, ub);

lambda = [param_opti_1(1), param_opti_2(1), param_opti_3(1), param_opti_4(1)]
sigma_f = [param_opti_1(2), param_opti_2(2), param_opti_3(2), param_opti_4(2)]
sigma_l = [param_opti_1(3), param_opti_2(3), param_opti_3(3), param_opti_4(3)]