function [GP_mu,GP_var] = GP_SE_mean_var(x, x_train, y_train, trueKern, noiseVar)
% Compute mean and variance of GP with squared exponential kernel, 
% 
% WARNING!: Only works for  x  is only one testing point meaning for many testing point
% this function needs to be called in the loop
% WARNING!: Current implementation is only for scalar, needs to  be
% extended for ND dimension of input case
%
% x : testing point
% x_train: training inputs
% y_train: training outputs
% trueKern : Kernel according to GPmat package
% noiseVar : Noise variance such that pdinv(K + noiseVar*eye(n))


% calculate cost for GP
K_cost = kernCompute(trueKern, x_train) + eye(size(x_train, 1))*noiseVar;
k_star_cost = kernCompute(trueKern, x(:,:), x_train);
% The means of the prediction
pdinv_K_cost = pdinv(K_cost);
alpha = pdinv_K_cost*y_train;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% mean of the GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GP_mu = k_star_cost*alpha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% variance of the GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variance of the prediction
GP_var = kernDiagCompute(trueKern, x(:)) - sum(k_star_cost*pdinv_K_cost.*k_star_cost, 2);


end