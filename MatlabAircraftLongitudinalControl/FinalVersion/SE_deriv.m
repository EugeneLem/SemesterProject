function dm = SE_deriv(z_star, Z_train, Y_train, noisevar, sigma_f, sigma_l)

kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));   

alpha= (kfcn(Z_train, Z_train, sigma_f, sigma_l) + noisevar* eye(size(Y_train,1))) \ Y_train ;
K_star = kfcn(z_star, Z_train, sigma_f, sigma_l);
X_tilde = z_star-Z_train;

lambda = 1/(sigma_l^2)*eye(length(z_star));

dm= -(lambda* X_tilde')*(K_star' .* alpha);


end