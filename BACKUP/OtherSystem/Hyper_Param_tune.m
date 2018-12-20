clear all;
load('data_for_opti.mat') %Data from a simulation with no GP

sigmaL_bound = [9,100]; %Found with some iteration

params_1 = hyperparameters('fitrgp',Z,Y_1);
params_1(4).Range = sigmaL_bound;
params_1(4).Optimize = true;        %Tell him to optimize KernelScale as well
gprMdl_1 = fitrgp(Z, Y_1, 'OptimizeHyperparameters', params_1);

params_2 = hyperparameters('fitrgp',Z,Y_2);
params_2(4).Range = sigmaL_bound;
params_2(4).Optimize = true;
gprMdl_2 = fitrgp(Z, Y_2, 'OptimizeHyperparameters', params_2);

params_3 = hyperparameters('fitrgp',Z,Y_3);
params_3(4).Range = sigmaL_bound;
params_3(4).Optimize = true;
gprMdl_3 = fitrgp(Z, Y_3, 'OptimizeHyperparameters', params_3);

params_4 = hyperparameters('fitrgp',Z,Y_4);
params_4(4).Range = sigmaL_bound;
params_4(4).Optimize = true;
gprMdl_4 = fitrgp(Z, Y_1, 'OptimizeHyperparameters', params_4);


%Get optimal value
sigma_l = [ gprMdl_1.KernelInformation.KernelParameters(1),  gprMdl_2.KernelInformation.KernelParameters(1), gprMdl_3.KernelInformation.KernelParameters(1), gprMdl_4.KernelInformation.KernelParameters(1)]
sigma_f = [ gprMdl_1.KernelInformation.KernelParameters(2),  gprMdl_2.KernelInformation.KernelParameters(2), gprMdl_3.KernelInformation.KernelParameters(2), gprMdl_4.KernelInformation.KernelParameters(2)]


