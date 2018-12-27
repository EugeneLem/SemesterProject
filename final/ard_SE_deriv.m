function dm = ard_SE_deriv(z_star, Z_train, Y_train, noisevar, sigma_f, sigma_l)

if length(sigma_l) == 1 %SE
    kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));   

    alpha= (kfcn(Z_train, Z_train, sigma_f, sigma_l) + noisevar* eye(size(Y_train,1))) \ Y_train ;
    K_star = kfcn(z_star, Z_train, sigma_f, sigma_l);
    X_tilde = z_star-Z_train;

    lambda = 1/(sigma_l^2)*eye(length(z_star));

    dm= -(lambda* X_tilde')*(K_star' .* alpha);

else  %SEARD
    sigma_l = sigma_l'; %transform to row vector
    kfcn_ard = @(XN,XM,sigmaF,sig_l ) (sigmaF^2)*exp(-(pdist2(XN./sig_l, XM./sig_l).^2)/2);

    alpha= (kfcn_ard(Z_train, Z_train, sigma_f, sigma_l) + noisevar* eye(size(Y_train,1))) \ Y_train ;
    K_star = kfcn_ard(z_star, Z_train, sigma_f, sigma_l);
    X_tilde = z_star-Z_train;

    lambda = diag(1./(sigma_l.^2));

    dm= -(lambda* X_tilde')*(K_star' .* alpha);   
    
    
end

end