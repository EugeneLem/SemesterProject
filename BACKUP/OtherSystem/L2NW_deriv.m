function dm = L2NW_deriv(z_star, Z_train, Y_train, param, type)
%type is for the type of kernel function used

switch type
    case 1
        %In case 1: k(v) = 1-|v|, param = [h,lambda]
        E = (pdist2(z_star,Z_train)/param(1)^2)'; 
        i = E<1;
        dm =  ( sum( -(2/param(1)^2)* Y_train(i).*(z_star - Z_train(i)) ,1)/( param(2) + sum(1-E(i),1) ) ...
            +(2/param(1)^2)*sum(z_star - Z_train(i),1) * sum( Y_train(i).*(1-E(i)) ,1)/( param(2)+sum(1-E(i),1) )^2 )' ; 

end