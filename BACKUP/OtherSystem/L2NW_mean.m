function m = L2NW_mean(z_star, Z_train, Y_train, param, type)
%type is for the type of kernel function used

switch type
    case 1
        %In case 1: k(v) = 1-|v|, param = [h,lambda]
        E = (pdist2(z_star,Z_train)/param(1)^2)'; 
        i = E<1;
        m = sum( Y_train(i).*(1-E(i)) ,1)/( param(2)+sum(1-E(i),1) );

end