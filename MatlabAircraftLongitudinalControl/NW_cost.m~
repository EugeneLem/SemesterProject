function obj = NW_cost(var,Z, Y1, L2NW_case)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
lambda = var(1);
sig_f = var(2);
sig_l = var(3:end);

obj = 0;
for i = 2:size(Z,1)
    obj = obj + (L2NW_mean(Z(i,:), Z(1:i-1,:), Y1(1:i-1), [sig_l, sig_f, lambda], L2NW_case)-Y1(i))^2;   
end
end

