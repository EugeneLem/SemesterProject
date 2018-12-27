function m = ard_L2NW_mean(z_star, Z_train, Y_train, param)%, type)
%type is for the type of kernel function used
%param = [lambda sigmaF sigmaL]
    if length(param) ==3 %only 1 sig_l
        sig_L = param(3);
        sig_F = param(2);
        lambda = param(1);
        kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));
        m = sum(kfcn(z_star,Z_train,sig_F, sig_L).*Y_train')/ (lambda+sum(kfcn(z_star,Z_train,sig_F, sig_L)));

    else    %seard
        sig_L = param(3:end);
        sig_F = param(2);
        lambda = param(1);
        kfcn_ard = @(XN,XM,sigmaF,sig_l ) (sigmaF^2)*exp(-(pdist2(XN./sig_l, XM./sig_l).^2)/2);
        m = sum(kfcn_ard(z_star,Z_train,sig_F, sig_L).*Y_train')/ (lambda+sum(kfcn_ard(z_star,Z_train,sig_F, sig_L)));                
    end

% switch type
%     case 1
%         %In case 1: k(v) = 1-|v|, param = [h,lambda]
%         E = (pdist2(z_star,Z_train)/param(1)^2)'; 
%         i = E<1;
%         m = sum( Y_train(i).*(1-E(i)) ,1)/( param(2)+sum(1-E(i),1) );
%         
%         case 2  %in this case, k(v) from GP, with param = [lambda sigmaF sigmaL]
%             if length(param) ==3 %only 1 sig_l
%                 sig_L = param(3);
%                 sig_F = param(2);
%                 lambda = param(1);
%                 kfcn = @(XN,XM,sigmaF, L) (sigmaF^2)*exp(-(pdist2(XN,XM).^2)/(2*L^2));
%                 m = sum(kfcn(z_star,Z_train,sig_F, sig_L).*Y_train')/ (lambda+sum(kfcn(z_star,Z_train,sig_F, sig_L)));
% 
%             else    %seard
%                 sig_L = param(3:end);
%                 sig_F = param(2);
%                 lambda = param(1);
%                 kfcn_ard = @(XN,XM,sigmaF,sig_l ) (sigmaF^2)*exp(-(pdist2(XN./sig_l, XM./sig_l).^2)/2);
%                 m = sum(kfcn_ard(z_star,Z_train,sig_F, sig_L).*Y_train')/ (lambda+sum(kfcn_ard(z_star,Z_train,sig_F, sig_L)));                
%             end
%         case 3
%             
%         %Casadi Adaptation of case 2
%                 sig_L = param(3);
%                 sig_F = param(2);
%                 lambda = param(1);
%             
%             import casadi.*
%             dist2 = [];
%             for i=1:size(Z_train,1) %Casadi doesn't like vecnorm and pdist2
%                 dist2 = [dist2; sum( (z_star-Z_train(i,:)).^2 ) ]; %replaced the norm for casadi, no square roo as we will square afterward
%             end
%             
%             kfcn_cas = @(dist,sigmaF, L) (sigmaF^2)*exp(-(dist)/(2*L^2));
%             m = sum(kfcn_cas(dist2,sig_F, sig_L).*Y_train)/ (lambda+sum(kfcn_cas(dist2, sig_F, sig_L)));


end