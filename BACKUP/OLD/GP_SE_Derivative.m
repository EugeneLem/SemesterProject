function grad_mu = GP_SE_mean_var(x, x_train, y_train, trueKern, noiseVar)
% Estimate the matrix A_k and B_k by deriving the norm
% 
% x : testing point
% x_train: training inputs
% y_train: training outputs
% trueKern : Kernel according to GPmat package
% noiseVar : Noise variance such that pdinv(K + noiseVar*eye(n))



K_cost = kernCompute(trueKern, x_train) + eye(size(x_train, 1))*noiseVar;
k_star_cost = kernCompute(trueKern, x(:,:), x_train);
% The means of the prediction
pdinv_K_cost = pdinv(K_cost);

alpha = pdinv_K_cost*y_train;


% yPred = k_star_cost*alpha;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Calculat gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% mean of the GP's gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = size(x_train,2);
n = size(x_train,1);

grad_mu_tmp = zeros(d, n);
del_k_del_x = zeros(d, n);


x_train_cons = x_train';

for i = 1:n
grad_mu_tmp(:,i) = -diag(trueKern.inverseWidth)*(x(:)-x_train_cons(:,i)); % since hyp2.cov =  [ln l1, ln l2, ..., ln \sigma] format
del_k_del_x(:,i) = grad_mu_tmp(:,i) * k_star_cost(i);
end

grad_mu = grad_mu_tmp*(k_star_cost'.*alpha);
end