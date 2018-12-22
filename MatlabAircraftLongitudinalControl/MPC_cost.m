function cost = MPC_cost(var, horizon, Q, R, R_delta)
%var = [x(1,k=1), ... , x(1, k=h), u2(1), ...,  u2(17)]
%      [x(2,k=1), ... , x(2, k=h), u2(2), ...,  u2(18)]
%      [x(3,k=1), ... , x(3, k=h), u2(3), ...,  u2(19)]
%      [x(4,k=1), ... , x(4, k=h), u2(4), ...,  u2(20)]

x = var(1:4, 1:horizon);
u = var(1:4, horizon+1:(horizon+ceil(horizon/4)) );
u = u(1:horizon)';

cost = x(:)'*Q*x(:) + u(:)'*R*u(:) + (u(2:end)-u(1:end-1))'*R_delta*(u(2:end)-u(1:end-1)); %

end

