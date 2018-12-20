function m = L2NW_mean(z_star, Z_train, Y_train, param, type)
%type is for the type of kernel function used

switch type
    case 1
        %In case 1: k(v) = 1-|v|, param = [h,lambda]
        E = (pdist2(z_star,Z_train)/param(1)^2)'; 
        i = E<1;
        m = sum( Y_train(i).*(1-E(i)) ,1)/( param(2)+sum(1-E(i),1) );
        
        case 2
        %in this case, k(v) from GP, with param = [sigmaL sigmaF lambda]
            kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));
            m = sum(kfcn(z_star,Z_train,param(2), param(1)).*Y_train')/ (param(3)+sum(kfcn(z_star,Z_train,param(2), param(1))));
        case 3
        %Casadi Adaptation of case 2
            import casadi.*
            dist = [];
            for i=1:size(Z_train,1) %Casadi doesn't like vecnorm and pdist2
                dist = [dist; sqrt(sum( (z_star-Z_train(i,:)).^2 )) ]; %replaced the norm for casadi, maybe to reverse
            end
            
            kfcn_cas = @(dist,sigmaF, L) (sigmaF^2)*exp(-(dist.^2)/(2*L^2));
            m = sum(kfcn_cas(dist,param(2), param(1)).*Y_train)/ (param(3)+sum(kfcn_cas(dist,param(2), param(1))));


end