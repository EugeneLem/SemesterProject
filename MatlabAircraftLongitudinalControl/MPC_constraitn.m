function [c,ceq] = MPC_constraitn(var, horizon, X1, U0, Z, Y_1, Y_2, Y_3, Y_4, lambda, sigma_f, sigma_l)
%var = [x(1,k=1), ... , x(1, k=h), u2(0), ...,  u2(16)]
%      [x(2,k=1), ... , x(2, k=h), u2(1), ...,  u2(17)]
%      [x(3,k=1), ... , x(3, k=h), u2(2), ...,  u2(18)]
%      [x(4,k=1), ... , x(4, k=h), u2(3), ...,  u2(19)]

x = var(1:4, 1:horizon);
u = var(1:4, horizon+1:(horizon+ceil(horizon/4)) );
u = u(1:horizon)';

param = [sigma_l(1), sigma_f(1), lambda(1);...
    sigma_l(2), sigma_f(2), lambda(2);...
    sigma_l(3), sigma_f(3), lambda(3);...
    sigma_l(4), sigma_f(4), lambda(4)];

ceq = [];%[x(:,1)-X1];%;...
    %U0-u(1)];           %/!\ u(1) is u0
% for i=2:horizon
%     ceq = [ceq;  x(:,i)- [L2NW_mean([x(:,i-1)',u(i)] , Z, Y_1, param(1,:), 3);...   %/!\ u(1) is u0
%                         L2NW_mean([x(:,i-1)',u(i)], Z, Y_2, param(2,:), 3);...
%                         L2NW_mean([x(:,i-1)',u(i)], Z, Y_3, param(3,:), 3);...
%                         L2NW_mean([x(:,i-1)',u(i)], Z, Y_4, param(4,:), 3)]];
% 
%     
% end

c = [];

end

