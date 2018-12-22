clear all;
format long
import casadi.*
opti = casadi.Opti();
%load('data_for_opti.mat') %Data from a simulation with no GP
load('Data_noFault.mat')    %size of Z to check
%load('Data_Faulty_L2NW.mat')




%---------------------------------------------------------------------------
% Z = Z(:,[true,true,true,true,false,true]);
% L2NW_deriv(z_star', Z', Y_1', [10.695781461698118, 1.677582380126192 ,10^(-10)], 2)

%kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));


sig_f =[  20.575591520093276,  17.878167302000445,  16.199466399148097,  17.938357317647856];
sig_l =[  99.966954902978543,  34.241347755939984,  78.546996930435938,  88.147249188528974];
%sig_f =   1.0e+02 *[ 2.639696625733967   0.170294093647814   0.085880808896654   0.101553751261396];
%sig_l =[ 99.991949107543803  25.640539532091616  34.284915666906144  25.862479178807501];

%% ----------------Y_1-------------------------------------------------------
import casadi.*

opti = casadi.Opti();

lambda = opti.variable(1, 1);
%sig_l = opti.variable(1, 1);
%sig_f = opti.variable(1, 1);


opti.subject_to( 0 <= lambda );
obj = MX.zeros(1,1);

for i = 2:size(Z,1)
    obj = obj + (L2NW_mean(Z(i,:), Z(1:i-1,:), Y_1(1:i-1), [sig_l ,sig_f, lambda], 3)-Y_1(i))^2;   
end
opti.minimize(obj);

ops = struct;
ops.ipopt.print_level = 0;
ops.ipopt.tol = 1e-3;
opti.solver('ipopt', ops);

sol = opti.solve();
lambda_1 = sol.value(lambda);
sig_l_1 = sol.value(sig_l);
sig_f_1 = sol.value(sig_f);
sol.value(obj)

%% ----------------Y_2-------------------------------------------------------
opti = casadi.Opti();

lambda = opti.variable(1, 1);
% sig_l = opti.variable(1, 1);
% sig_f = opti.variable(1, 1);

opti.subject_to( 0 <= lambda );
obj = MX.zeros(1,1);

for i = 2:size(Z,1)
    obj = obj + (L2NW_mean(Z(i,:), Z(1:i-1,:), Y_2(1:i-1), [sig_l(2),sig_f(2), lambda], 3)-Y_2(i))^2;   
end
opti.minimize(obj);

ops = struct;
ops.ipopt.print_level = 0;
ops.ipopt.tol = 1e-3;
opti.solver('ipopt', ops);

sol = opti.solve();
lambda_2 = sol.value(lambda);
sig_l_2 = sol.value(sig_l);
sig_f_2 = sol.value(sig_f);
sol.value(obj)
%----------------Y_3-------------------------------------------------------
opti = casadi.Opti();

lambda = opti.variable(1, 1);
% sig_l = opti.variable(1, 1);
% sig_f = opti.variable(1, 1);

opti.subject_to( 0 <= lambda );
obj = MX.zeros(1,1);

for i = 2:size(Z,1)
    obj = obj + (L2NW_mean(Z(i,:), Z(1:i-1,:), Y_3(1:i-1), [sig_l(3),sig_f(3), lambda], 3)-Y_3(i))^2;   
end
opti.minimize(obj);

ops = struct;
ops.ipopt.print_level = 0;
ops.ipopt.tol = 1e-3;
opti.solver('ipopt', ops);

sol = opti.solve();
lambda_3 = sol.value(lambda);
sig_l_3 = sol.value(sig_l);
sig_f_3= sol.value(sig_f);
sol.value(obj)

%----------------Y_4-------------------------------------------------------
opti = casadi.Opti();

lambda = opti.variable(1, 1);
% sig_l = opti.variable(1, 1);
% sig_f = opti.variable(1, 1);

opti.subject_to( 0 <= lambda );
obj = MX.zeros(1,1);

for i = 2:size(Z,1)
    obj = obj + (L2NW_mean(Z(i,:), Z(1:i-1,:), Y_4(1:i-1), [sig_l(4),sig_f(4), lambda], 3)-Y_1(i))^2;   
end
opti.minimize(obj);

ops = struct;
ops.ipopt.print_level = 0;
ops.ipopt.tol = 1e-3;
opti.solver('ipopt', ops);

sol = opti.solve();
lambda_4 = sol.value(lambda);
sig_l_4 = sol.value(sig_l);
sig_f_4 = sol.value(sig_f);
sol.value(obj)

%--------------------
lambda_opt = [lambda_1, lambda_2, lambda_3, lambda_4]
% sigma_l_opt = [sig_l_1,sig_l_2,sig_l_3,sig_l_4]
% sigma_f_opt = [sig_f_1,sig_f_2,sig_f_3,sig_f_4]
