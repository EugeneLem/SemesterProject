function cost = opti_GPse_cost(param, Z, Y)

noisevar = param(1);
sig_f = param(2);
sig_l = param(3);

cost = 0;
for i=2:size(Y,1)
    cost = cost + (Y(i) - ...
        ard_SE_mean(Z(i,:), Z(1:i-1,:), Y(1:i-1), noisevar, sig_f, sig_l) )^2;  
end

end


               