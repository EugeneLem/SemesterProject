format long

%%System Constant==========================================================
%System
A=[0.86719   ,6.6936*10^(-5)     ,-0.19095  ,0; ...
    -0.027773, 0.99895           , 0.89264, -1.9609; ...
      0.20146, -2.1676*10^(-4)     , 0.88379,  0; ...
      0.2    ,  0                  ,  0      , 1];
  
A_Fault=[0.86719   ,6.6936*10^(-5)     ,-0.19095  ,0; ...
    -0.027773, 0.99895           , 0.89264, -1.9609; ...
      0.20146, -2.1676*10^(-4)     , 0.88379,  0; ...
      0.2    ,  0                  ,  0      , 1]; 
  
B=[-3.7758 *10^(-3)  , -9.0408*10^(-3);...
    0                ,0;...
      -1.2629*10^(-4), -3.2794*10^(-4);...
     0              ,0];
 
B_Fault=[3.7758 *10^(-3)  , -9.0408*10^(-3);...
    0                ,0;...
      1.2629*10^(-4), -3.2794*10^(-4);...
     0              ,0]; 
 
 
%Cost function 
Q=diag([34.38, 0, 0, 103.13]);
R=zeros(2,2);
R_delta=diag([0.02,10]);
P_k=zeros(6,6);
%horizon
horizon=20;
Ts=0.2;
%constraint
C_u = [1,0 ; 0,1 ; -1,0 ; 0,-1];
c_u = [15;4.6;25;10.4];%Need to check is elevator and THS are not inverted
C_u_delta = [-1 0   1 0;...
            0 -1  0 1;...
            1  0 -1 0;...
            0  1  0 -1];
c_u_delta = [37;0.236;37;0.236]*Ts;%Need to check is elevator and THS are not inverted

% Nominal MPC =============================================================

x=sdpvar(4,horizon);
u=sdpvar(2,horizon);
u_old=sdpvar(2,1);%For the input change constraint, we need the 0 input
%ops = sdpsettings( 'verbose', 0);
ops = sdpsettings('solver', 'mosek','verbose', 0);
con = C_u_delta*[u_old;u(:,1)] <= c_u_delta;
obj = 0;   
for i = 1:horizon-1 
    con = [con , x(:,i+1) == A*x(:,i) + B*u(:,i)]; % System dynamics
    con = [con , C_u*u(:,i) <= c_u];                 % Input constraints
    con = [con , C_u_delta*[u(:,i);u(:,i+1)] <= c_u_delta];  % Input change constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i) + (u(:,i+1)-u(:,i))'*R_delta*(u(:,i+1)-u(:,i));  % Cost function

end
obj = obj+[u(:,horizon);x(:,horizon)]'*P_k*[u(:,horizon);x(:,horizon)];
ctrl_nominal = optimizer(con, obj,ops, {x(:,1),u_old}, u(:,1));    


% Broken MPC =============================================================
x_hat=sdpvar(4,horizon);
u_hat=sdpvar(2,horizon);
A_k = sdpvar(4,4,'full');
B_k = sdpvar(4,2);
d_k = sdpvar(4,1);
u_old_hat=sdpvar(2,1);%For the input change constraint, we need the 0 input
%ops = sdpsettings( 'verbose', 0);
ops = sdpsettings('solver', 'mosek','verbose', 0);
con_hat = C_u_delta*[u_old_hat;u_hat(:,1)] <= c_u_delta;
obj_hat = 0;   
for i = 1:horizon-1 
    con_hat = [con_hat , x_hat(:,i+1) == A_k*x_hat(:,i) + B_k*u_hat(:,i) + d_k]; % System dynamics
    con_hat = [con_hat , C_u*u_hat(:,i) <= c_u];                 % Input constraints
    con_hat = [con_hat , C_u_delta*[u_hat(:,i);u_hat(:,i+1)] <= c_u_delta];  % Input change constraints
    obj_hat = obj_hat + x_hat(:,i)'*Q*x_hat(:,i) + u_hat(:,i)'*R*u_hat(:,i) + (u_hat(:,i+1)-u_hat(:,i))'*R_delta*(u_hat(:,i+1)-u_hat(:,i));  % Cost function

end
obj_hat = obj_hat+[u_hat(:,horizon);x_hat(:,horizon)]'*P_k*[u_hat(:,horizon);x_hat(:,horizon)];
ctrl_broken = optimizer(con_hat, obj_hat,ops, {x_hat(:,1),u_old_hat, A_k, B_k, d_k}, u_hat(:,1));





% Simulation constant =====================================================
stepNumber=300;
k_fault=2;
k_switch = k_fault+2;
X0 = [1;0.1;1;1];
U0=[0;0];

X_noFault = zeros(4,stepNumber); 
U_noFault=zeros(2,stepNumber-1); 
X_noFault(:,1) = X0;
u_noFault(:,1) = U0;

X_Fault = zeros(4,stepNumber); 
U_Fault = zeros(2,stepNumber-1); 
X_Fault(:,1) = X0;
u_Fault(:,1) = U0;

X_Fault_Perfect = zeros(4,stepNumber); 
U_Fault_Perfect = zeros(2,stepNumber-1); 
X_Fault_Perfect(:,1) = X0;
u_Fault_Perfect(:,1) = U0;

X_Fault_GP = zeros(4,stepNumber); 
U_Fault_GP = zeros(2,stepNumber-1); 
X_Fault_GP(:,1) = X0;
u_Fault_GP(:,1) = U0;
error_GP = zeros(4,stepNumber-1);

X_Fault_L2NW = zeros(4,stepNumber); 
U_Fault_L2NW = zeros(2,stepNumber-1); 
X_Fault_L2NW(:,1) = X0;
u_Fault_L2NW(:,1) = U0;
error_L2NW = zeros(4,stepNumber-1);


% System simu =============================================================

for k=1:stepNumber-1                    
    %----------no fault---------------
    if k==1
        [U_noFault(:,k),infeasible]=ctrl_nominal{X_noFault(:,k),U0};
    else
        [U_noFault(:,k),infeasible]=ctrl_nominal{X_noFault(:,k),U_noFault(:,k-1)};
    end
    X_noFault(:,k+1)=A*X_noFault(:,k)+B*U_noFault(:,k);
    
    
    %----------fault with old MPC---------------
    if k==1
        [U_Fault(:,k),infeasible] = ctrl_nominal{X_Fault(:,k),U0};
        X_Fault(:,k+1)=A*X_Fault(:,k)+B*U_Fault(:,k);
    elseif k < k_fault
        [U_Fault(:,k),infeasible]=ctrl_nominal{X_Fault(:,k),U_Fault(:,k-1)};
        X_Fault(:,k+1)=A*X_Fault(:,k)+B*U_Fault(:,k);
    else
        [U_Fault(:,k),infeasible] = ctrl_nominal{X_Fault(:,k),U_Fault(:,k-1)};
        X_Fault(:,k+1)=A_Fault*X_Fault(:,k)+B_Fault*U_Fault(:,k);
    end
    
    %----------fault with perfect MPC-------------
    if k==1
        [U_Fault_Perfect(:,k),infeasible] = ctrl_nominal{X_Fault_Perfect(:,k),U0};
        X_Fault_Perfect(:,k+1)=A*X_Fault_Perfect(:,k)+B*U_Fault_Perfect(:,k);
    elseif k < k_fault
        [U_Fault_Perfect(:,k),infeasible]=ctrl_nominal{X_Fault_Perfect(:,k),U_Fault_Perfect(:,k-1)};
        X_Fault_Perfect(:,k+1)=A*X_Fault_Perfect(:,k)+B_Fault*U_Fault_Perfect(:,k);
    else
        [U_Fault_Perfect(:,k),infeasible] = ctrl_broken{X_Fault_Perfect(:,k),U_Fault_Perfect(:,k-1), A_Fault, B_Fault, zeros(4,1)};
        X_Fault_Perfect(:,k+1)=A_Fault*X_Fault_Perfect(:,k)+B_Fault*U_Fault_Perfect(:,k);
    end
    
    
    
        %----------fault with GP-------------
    if k==1
        [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),U0};
        X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*U_Fault_GP(:,k);
    elseif k < k_fault
        [U_Fault_GP(:,k),infeasible]=ctrl_nominal{X_Fault_GP(:,k),U_Fault_GP(:,k-1)};
        X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*U_Fault_GP(:,k);
    elseif k < k_switch
        [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),U_Fault_GP(:,k-1)};
        X_Fault_GP(:,k+1)=A_Fault*X_Fault_GP(:,k)+B_Fault*U_Fault_GP(:,k);
    else
        Y_GP_1 = X_Fault_GP(1,k_fault+1:k-1)';
        Y_GP_2 = X_Fault_GP(2,k_fault+1:k-1)';
        Y_GP_3 = X_Fault_GP(3,k_fault+1:k-1)';
        Y_GP_4 = X_Fault_GP(4,k_fault+1:k-1)';
        Z_GP = [X_Fault_GP(:,k_fault:k-2)', U_Fault_GP(:,k_fault:k-2)'];
        
        noisevar = 0.00001;
        sigma_f = 2.17*10^(3)*ones(1,4);
        sigma_l = 1.14*10^(4)*ones(1,4);
%         sigma_f = [1.677582380126192, 80.087113811947248, 1.371605612665715, 6.573873522600369];
%         sigma_l = [10.695781461698118, 99.935986868827754,  9.046354763354064, 32.505687565655741];
%         noisevar = 0.00000001;
%         
        X_estimation = [SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1));...
            SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2));...
            SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3));...
            SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))];
        
        
        error_GP(:,k) = abs(X_estimation-X_Fault_GP(:,k)) ;
        
        if  max(error_GP(:,k)) < 0.1
%            disp("Validated")
            
            AB_k = [SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1))';...
                SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2))';...
                SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3))';...
                SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(:,k-1)'], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))'
                ];
            A_k = AB_k(:,1:4);
            B_k = AB_k(:,5:6);
            d_k = X_Fault_GP(:,k) - (A_k*X_Fault_GP(:,k-1) + B_k*U_Fault_GP(:,k-1));
            
            
            [U_Fault_GP(:,k),infeasible] = ctrl_broken{X_Fault_GP(:,k),U_Fault_GP(:,k-1), A_k, B_k, d_k};
            X_Fault_GP(:,k+1)=A_Fault*X_Fault_GP(:,k)+B_Fault*U_Fault_GP(:,k);
        else
            disp("NOT Validated")
            %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
            [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),U_Fault_GP(:,k-1)};
            X_Fault_GP(:,k+1)=A_Fault*X_Fault_GP(:,k)+B_Fault*U_Fault_GP(:,k);
        end
    end
        
        
        
        
        %----------fault with L2NW-------------
    if k==1
        [U_Fault_L2NW(:,k),infeasible] = ctrl_nominal{X_Fault_L2NW(:,k),U0};
        X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*U_Fault_L2NW(:,k);
    elseif k < k_fault
        [U_Fault_L2NW(:,k),infeasible]=ctrl_nominal{X_Fault_L2NW(:,k),U_Fault_L2NW(:,k-1)};
        X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*U_Fault_L2NW(:,k);
    elseif k < k_switch
        [U_Fault_L2NW(:,k),infeasible] = ctrl_nominal{X_Fault_L2NW(:,k),U_Fault_L2NW(:,k-1)};
        X_Fault_L2NW(:,k+1)=A_Fault*X_Fault_L2NW(:,k)+B_Fault*U_Fault_L2NW(:,k);
    else
        param = [5, 10; 5, 10; 5, 10; 5, 10]; %[h, lambda]
        
        Z_L2NW = [X_Fault_L2NW(:,k_fault:k-2)', U_Fault_L2NW(:,k_fault:k-2)'];
        Y_L2NW_1 = X_Fault_L2NW(1,k_fault+1:k-1)' - ((A(1,:)*Z_L2NW(:,1:4)' + B(1,:)*Z_L2NW(:,5:6)'))';
        Y_L2NW_2 = X_Fault_L2NW(2,k_fault+1:k-1)' - ((A(2,:)*Z_L2NW(:,1:4)' + B(2,:)*Z_L2NW(:,5:6)'))';
        Y_L2NW_3 = X_Fault_L2NW(3,k_fault+1:k-1)' - ((A(3,:)*Z_L2NW(:,1:4)' + B(3,:)*Z_L2NW(:,5:6)'))';
        Y_L2NW_4 = X_Fault_L2NW(4,k_fault+1:k-1)' - ((A(4,:)*Z_L2NW(:,1:4)' + B(4,:)*Z_L2NW(:,5:6)'))';
        
        
        X_estimation_L2NW = [L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_1,param(1,:),1);...
            L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_2,param(2,:),1);...
            L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_3,param(3,:),1);...
            L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_4,param(4,:),1)...
            ];
        
        
        error_L2NW(:,k) = abs([X_Fault_L2NW(:,k) -  X_estimation_L2NW]);
        
        if  max(error_L2NW(:,k)) < 300

            O_k = [L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_1,param(1,:),1)';...
                L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_2,param(2,:),1)';...
                L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_3,param(3,:),1)';...
                L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(:,k-1)'],Z_L2NW,Y_L2NW_4,param(4,:),1)';...
                ];
            O_A_k = O_k(:,1:4);
            O_B_k = O_k(:,5:6);
            O_d_k = X_Fault_L2NW(:,k) - (A*X_Fault_L2NW(:,k-1) + B*U_Fault_L2NW(:,k-1) + O_A_k*X_Fault_L2NW(:,k-1) + O_B_k*U_Fault_L2NW(:,k-1));


            [U_Fault_L2NW(:,k), infeasible] = ctrl_broken(X_Fault_L2NW(:,k),U_Fault_L2NW(:,k-1), A+O_A_k, B+O_B_k, O_d_k);
            if infeasible ~= 0
                disp('problem')
            end

            X_Fault_L2NW(:,k+1)=A_Fault*X_Fault_L2NW(:,k)+B_Fault*U_Fault_L2NW(:,k);
        else
            disp('L2NW  NOT validated')
            [U_Fault_L2NW(:,k),infeasible] = ctrl_nominal{X_Fault_L2NW(:,k),U_Fault_L2NW(:,k-1)};
            X_Fault_L2NW(:,k+1)=A_Fault*X_Fault_L2NW(:,k)+B_Fault*U_Fault_L2NW(:,k);
            
            
        end
    end
      
        
    %---------------ALL----------------------------------
    if mod(k,10) == 0
        disp(k)
    end
    
end




%%
% Plotting ================================================================
t=Ts*(0:stepNumber-1);

figure
subplot(2,2,1)
plot(t,X_noFault(1,:),t,X_Fault(1,:),t,X_Fault_Perfect(1,:),'o',t,X_Fault_GP(1,:),t,X_Fault_L2NW(1,:))
title('Pitch Angle rad/s')
xlabel('time s')
subplot(2,2,2)
plot(t,X_noFault(2,:),t,X_Fault(2,:),t,X_Fault_Perfect(2,:),'o',t,X_Fault_GP(2,:),t,X_Fault_L2NW(2,:))
title('Air speed m/s')
xlabel('time s')
subplot(2,2,3)
plot(t,X_noFault(3,:),t,X_Fault(3,:),t,X_Fault_Perfect(3,:),'o',t,X_Fault_GP(3,:),t,X_Fault_L2NW(3,:))
title('Pitch Rate rad')
xlabel('time s')
subplot(2,2,4)
plot(t,X_noFault(4,:),t,X_Fault(4,:),t,X_Fault_Perfect(4,:),'o',t,X_Fault_GP(4,:),t,X_Fault_L2NW(4,:))
title('AoA rad')
xlabel('time s')
legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP', 'L2NW_1')

figure
subplot(2,1,1)
plot(t(1:end-1),U_noFault(1,:),t(1:end-1),U_Fault(1,:),...
    t(1:end-1),U_Fault_Perfect(1,:),'o',t(1:end-1),U_Fault_GP(1,:),t(1:end-1),U_Fault_L2NW(1,:))
title('Elevator deg')
xlabel('time s')
subplot(2,1,2)
plot(t(1:end-1),U_noFault(2,:),t(1:end-1),U_Fault(2,:)...
    ,t(1:end-1),U_Fault_Perfect(2,:),'o',t(1:end-1),U_Fault_GP(2,:),t(1:end-1),U_Fault_L2NW(2,:))
title('THS deg ')
xlabel('time s')
legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP')

% figure
% subplot(2,1,1)
% plot(t(1:end-1),matrix_ratio(1,:))
% title('A_k')
% subplot(2,1,2)
% plot(t(1:end-1),matrix_ratio(2,:))
% title('B_k')


dU_noFault = (U_noFault(:,2:end)-U_noFault(:,1:end-1))/Ts;
dU_Fault = (U_Fault(:,2:end)-U_Fault(:,1:end-1))/Ts;
dU_Fault_Perfect = (U_Fault_Perfect(:,2:end)-U_Fault_Perfect(:,1:end-1))/Ts;
dU_Fault_GP = (U_Fault_GP(:,2:end)-U_Fault_GP(:,1:end-1))/Ts;
dU_Fault_L2NW = (U_Fault_L2NW(:,2:end)-U_Fault_L2NW(:,1:end-1))/Ts;

figure
subplot(2,1,1)
plot(t(1:end-2),dU_noFault(1,:),t(1:end-2),dU_Fault(1,:),...
    t(1:end-2),dU_Fault_Perfect(1,:),'o',t(1:end-2),dU_Fault_GP(1,:),t(1:end-2),dU_Fault_L2NW(1,:))
title('Elevator slew rate def/s')
xlabel('time s')
subplot(2,1,2)
plot(t(1:end-2),dU_noFault(2,:),t(1:end-2),dU_Fault(2,:),...
    t(1:end-2),dU_Fault_GP(2,:),'o',t(1:end-2),dU_Fault_GP(2,:),t(1:end-2),dU_Fault_L2NW(2,:))
title('THS slew rate deg/s')
xlabel('time s')
legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP')

figure

subplot(2,2,1)
plot(t(1:end-1),error_GP(1,:),t(1:end-1),error_L2NW(1,:))
title('Error on Pitch Angle rad/s')
xlabel('time s')
legend('GP', 'L2NW 1')
subplot(2,2,2)
plot(t(1:end-1),error_GP(2,:),t(1:end-1),error_L2NW(2,:))
title('Error on Air speed m/s')
xlabel('time s')
legend('GP', 'L2NW 1')
subplot(2,2,3)
plot(t(1:end-1),error_GP(3,:),t(1:end-1),error_L2NW(3,:))
title('Error on Pitch Rate rad')
xlabel('time s')
legend('GP', 'L2NW 1')
subplot(2,2,4)
plot(t(1:end-1),error_GP(4,:),t(1:end-1),error_L2NW(4,:))
title('Error on AoA rad')
xlabel('time s')
legend('GP', 'L2NW 1')



