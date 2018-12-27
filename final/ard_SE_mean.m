function m = ard_SE_mean(z_star, Z_train, Y_train, noisevar, sigma_f, sigma_l)
    
if length(sigma_l) == 1 %SE
    kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));
    
    K_star = kfcn(z_star, Z_train, sigma_f, sigma_l);
    K = kfcn(Z_train, Z_train, sigma_f, sigma_l) + noisevar* eye(size(Y_train,1)) ;
    
    m=K_star*(K\Y_train);
    
    
else  %SEARD
    sigma_l = sigma_l'; %transform to row vector
    kfcn_ard = @(XN,XM,sigmaF,sig_l ) (sigmaF^2)*exp(-(pdist2(XN./sig_l, XM./sig_l).^2)/2);
    
    K_star = kfcn_ard(z_star, Z_train, sigma_f, sigma_l);
    K = kfcn_ard(Z_train, Z_train, sigma_f, sigma_l) + noisevar* eye(size(Y_train,1)) ;
    
    m=K_star*(K\Y_train);  
    
    
end

end