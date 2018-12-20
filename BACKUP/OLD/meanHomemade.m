function m = meanHomemade(z_star, Z_train, Y_train, noisevar, sigma_f, sigma_l)
    
    kfcn = @(XN,XM,theta) (theta(2)^2)*exp(-(pdist2(XN,XM).^2)/(2*theta(1)^2));
    
    K_star = kfcn(z_star, Z_train, [sigma_f, sigma_l]);
    K = kfcn(Z_train, Z_train, [sigma_f, sigma_l]) + noisevar* eye(size(Y_train,1)) ;
    
    m=K_star*(K\Y_train);


end