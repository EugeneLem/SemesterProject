function dm = L2NW_deriv(z_star, Z_train, Y_train, param, type)
%type is for the type of kernel function used

switch type
    case 1
        %In case 1: k(v) = 1-|v|, param = [h,lambda]
        E = (pdist2(z_star,Z_train)/param(1)^2)'; 
        i = E<1;
        dm =  ( sum( -(2/param(1)^2)* Y_train(i).*(z_star - Z_train(i)) ,1)/( param(2) + sum(1-E(i),1) ) ...
            +(2/param(1)^2)*sum(z_star - Z_train(i),1) * sum( Y_train(i).*(1-E(i)) ,1)/( param(2)+sum(1-E(i),1) )^2 )' ; 
    case 2
    %in this case, k(v) from GP, with param = [sigmaL sigmaF lambda]
    kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));
    k = kfcn(z_star,Z_train,param(2), param(1)); % k(i) = k(x_star - Z(i)) 
    dk = -1/(param(1)^2)*(z_star-Z_train).*k';   % dk(i) = d( k(x_star - Z(i))) / d(x_star)
    
    dm = sum(Y_train'.*dk',2)/(param(3)+sum(k))  -  (sum(dk,1)*sum(Y_train'.*k)/(param(3) + sum(k))^2)';
        
end