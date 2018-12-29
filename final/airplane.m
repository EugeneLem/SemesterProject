format long

load('Data_Fault_GP.mat') %load data from when there are no test
Z = Z(:, [1,2,3,4,6]);   %For faulty actuator, we do not want the first input   


do_noFault              = false;   %list of calculation we want to do:
do_Fault                = false;
do_Fault_perfect        = false;
do_Fault_GP             = true;
do_Fault_GP_moreData    = true;
do_L2NW_GP              = false;
do_Fault_GPard          = true;
do_Fault_GPard_moreData = false;
do_L2NW_GPard           = false;

%% System Constant==========================================================
%System
A=[0.86719   ,6.6936*10^(-5)     ,-0.19095  ,0; ...
    -0.027773, 0.99895           , 0.89264, -1.9609; ...
      0.20146, -2.1676*10^(-4)     , 0.88379,  0; ...
      0.2    ,  0                  ,  0      , 1];
  
B=[-3.7758 *10^(-3)  , -9.0408*10^(-3);...
    0                ,0;...
      -1.2629*10^(-4), -3.2794*10^(-4);...
     0              ,0];
 
%Cost function 
Q=diag([34.38, 0, 0, 103.13]);
R=zeros(2,2);
R_delta=diag([0.02,10]);
P_k=zeros(6,6);

%horizon
horizon=20;
Ts=0.2;


Q_huge = blkdiag(kron(eye(horizon-1), Q),P_k(1:4,1:4));
R_huge = blkdiag(kron(eye(horizon-1), R(2,2)),P_k(6,6)); %Only P_k(6,6) because we only are interrested in u(2), and other at 0
R_delta_huge = blkdiag(kron(eye(horizon-1), R_delta(2,2)));

%constraint
C_u = [1,0 ; 0,1 ; -1,0 ; 0,-1];        %Constraint on input
c_u = [15;4.6;25;10.4]; 
C_u_delta = [-1 0   1 0;...
            0 -1  0 1;...
            1  0 -1 0;...
            0  1  0 -1];
c_u_delta = [37;0.236;37;0.236]*Ts;% constraint on input change rate


C_u_redu = [ 1 ; -1];       %redu are when we only regulate one input: no need to constraint the stuck one
c_u_redu = [4.6;10.4];
C_u_delta_redu = [-1  1; 1  -1];
c_u_delta_redu = [0.236; 0.236]*Ts;

%Inequality constraint for fmincon
C_ineq_2_2 = zeros(2*horizon-2,horizon); %the end of C_ineq matrix
for i=1:horizon-1
    C_ineq_2_2(2*i-1:2*i,i:i+1)  = C_u_delta_redu;

end

C_ineq = [zeros(4*horizon-2, 4*horizon),...
    [kron(eye(horizon),C_u_redu);C_ineq_2_2]];
c_ineq = [kron(ones(horizon,1), c_u_redu); kron(ones(horizon-1,1), c_u_delta_redu)];

%% Nominal MPC =============================================================

x=sdpvar(4,horizon);
u=sdpvar(2,horizon);
u_old=sdpvar(2,1);%For the input change constraint, we need the 0 input
%ops = sdpsettings( 'verbose', 0);
ops = sdpsettings('solver', 'mosek','verbose', 0);
con = C_u_delta*[u_old;u(:,1)] <= c_u_delta;    %First input change rate constraint
obj = 0;   
for i = 1:horizon-1 
    con = [con , x(:,i+1) == A*x(:,i) + B*u(:,i)];          % System dynamics
    con = [con , C_u*u(:,i) <= c_u];                        % Input constraints
    con = [con , C_u_delta*[u(:,i);u(:,i+1)] <= c_u_delta]; % Input change constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i) + (u(:,i+1)-u(:,i))'*R_delta*(u(:,i+1)-u(:,i));  % Cost function

end
obj = obj+[u(:,horizon);x(:,horizon)]'*P_k*[u(:,horizon);x(:,horizon)];
ctrl_nominal = optimizer(con, obj,ops, {x(:,1),u_old}, u(:,1));    

%% Broken MPC =============================================================
%for when we linearize the system
x_hat=sdpvar(4,horizon);
u_hat=sdpvar(2,horizon);
A_k = sdpvar(4,4,'full');       %Now the identified matrix are input of the problem
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

%% Simulation constant =====================================================
stepNumber = 1000;       %HERE
k_fault=2;               %We initialize everything to zero
k_switch = k_fault+2;
X0 = [1;0.1;1;1];
U0=[0;0];

% X_noFault = zeros(4,stepNumber); 
% U_noFault=zeros(2,stepNumber-1); 
% X_noFault(:,1) = X0;
% U_noFault(:,1) = U0;
% 
% X_Fault = zeros(4,stepNumber); 
% U_Fault = zeros(2,stepNumber-1); 
% X_Fault(:,1) = X0;
% U_Fault(:,1) = U0;
% 
% X_Fault_Perfect = zeros(4,stepNumber); 
% U_Fault_Perfect = zeros(2,stepNumber-1); 
% X_Fault_Perfect(:,1) = X0;
% U_Fault_Perfect(:,1) = U0;

X_Fault_GP = zeros(4,stepNumber); 
U_Fault_GP = zeros(2,stepNumber-1); 
X_Fault_GP(:,1) = X0;
U_Fault_GP(:,1) = U0;
error_GP = zeros(4,stepNumber-1);
error_mat_GP = zeros(1,stepNumber-1);

X_Fault_GP_moreData = zeros(4,stepNumber); 
U_Fault_GP_moreData = zeros(2,stepNumber-1); 
X_Fault_GP_moreData(:,1) = X0;
U_Fault_GP_moreData(:,1) = U0;
error_GP_moreData = zeros(4,stepNumber-1);
error_mat_GP_moreData = zeros(1,stepNumber-1);

X_Fault_GPard = zeros(4,stepNumber); 
U_Fault_GPard = zeros(2,stepNumber-1); 
X_Fault_GPard(:,1) = X0;
U_Fault_GPard(:,1) = U0;
error_GPard = zeros(4,stepNumber-1);
error_mat_GPard = zeros(1,stepNumber-1);

X_Fault_GPard_moreData = zeros(4,stepNumber); 
U_Fault_GPard_moreData = zeros(2,stepNumber-1); 
X_Fault_GPard_moreData(:,1) = X0;
U_Fault_GPard_moreData(:,1) = U0;
error_GPard_moreData = zeros(4,stepNumber-1);
error_mat_GPard_moreData = zeros(1,stepNumber-1);

X_Fault_L2NWGP = zeros(4,stepNumber); 
U_Fault_L2NWGP = zeros(2,stepNumber-1); 
X_Fault_L2NWGP(:,1) = X0;
U_Fault_L2NWGP(:,1) = U0;
error_L2NWGP = zeros(4,stepNumber-1);
z_new = zeros(4,horizon+5);

X_Fault_L2NWGP_ard = zeros(4,stepNumber); 
U_Fault_L2NWGP_ard = zeros(2,stepNumber-1); 
X_Fault_L2NWGP_ard(:,1) = X0;
U_Fault_L2NWGP_ard(:,1) = U0;
error_L2NWGP_ard = zeros(4,stepNumber-1);
z_new_ard = zeros(4,horizon+5);


%% ==System simu =============================================================

for k=1:stepNumber-1                    
   %% ----------no fault--------------- 
    %This is a reference case where there is no problem on the airplane
    if do_noFault == true
        if k==1
            [U_noFault(:,k),infeasible]=ctrl_nominal{X_noFault(:,k),U0};
        else
            [U_noFault(:,k),infeasible]=ctrl_nominal{X_noFault(:,k),U_noFault(:,k-1)};
        end
        X_noFault(:,k+1)=A*X_noFault(:,k)+B*U_noFault(:,k);
	end
    
   %% ----------fault with old MPC---------------
        %This is a reference case, where the control continue the same dispite
    %the fault    
    if do_Fault == true
        if k==1
            [U_Fault(:,k),infeasible] = ctrl_nominal{X_Fault(:,k),U0};
            X_Fault(:,k+1)=A*X_Fault(:,k)+B*U_Fault(:,k);
        elseif k < k_fault
            [U_Fault(:,k),infeasible]=ctrl_nominal{X_Fault(:,k),U_Fault(:,k-1)};
            X_Fault(:,k+1)=A*X_Fault(:,k)+B*U_Fault(:,k);
        else
            [U_Fault(:,k),infeasible] = ctrl_nominal{X_Fault(:,k),[U_Fault(1,k_fault-1);U_Fault(2,k-1)]};
            X_Fault(:,k+1)=A*X_Fault(:,k)+B*[U_Fault(1,k_fault-1);U_Fault(2,k)];
        end
    end
    
   %% ----------fault with perfect MPC-------------
    %This is a reference case, where the controller know excatly what the
    %fault is (A stay the same, B is [zeros(4,1), B(:,2)], and d_k is
    %B(:,1)*u(1,k_fault-1) )
    if do_Fault_perfect == true
        if k==1
            [U_Fault_Perfect(:,k),infeasible] = ctrl_nominal{X_Fault_Perfect(:,k),U0};
            X_Fault_Perfect(:,k+1)=A*X_Fault_Perfect(:,k)+B*U_Fault_Perfect(:,k);
        elseif k < k_fault
            [U_Fault_Perfect(:,k),infeasible]=ctrl_nominal{X_Fault_Perfect(:,k),U_Fault_Perfect(:,k-1)};
            X_Fault_Perfect(:,k+1)=A*X_Fault_Perfect(:,k)+B*U_Fault_Perfect(:,k);
        else
            [U_Fault_Perfect(:,k),infeasible] = ctrl_broken{X_Fault_Perfect(:,k),[U_Fault_Perfect(1,k_fault-1);U_Fault_Perfect(2,k-1)], A, [zeros(4,1), B(:,2)], B(:,1)*U_Fault_Perfect(1,k_fault-1)};
            X_Fault_Perfect(:,k+1)=A*X_Fault_Perfect(:,k)+B*[U_Fault_Perfect(1,k_fault-1);U_Fault_Perfect(2,k)];
        end
    end   
    
   %% ----------fault with GPse-------------
        %At the moment of the fault, we use a GPse linearisation of the
        %identified systme
    if do_Fault_GP ==true
        if k==1
            [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),U0};
            X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*U_Fault_GP(:,k);
        elseif k < k_fault
            [U_Fault_GP(:,k),infeasible]=ctrl_nominal{X_Fault_GP(:,k),U_Fault_GP(:,k-1)};
            X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*U_Fault_GP(:,k);
        elseif k < k_switch
            [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k-1)]};
            X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k)];
        else
            Y_GP_1 = X_Fault_GP(1,k_fault+1:k-1)';
            Y_GP_2 = X_Fault_GP(2,k_fault+1:k-1)';
            Y_GP_3 = X_Fault_GP(3,k_fault+1:k-1)';
            Y_GP_4 = X_Fault_GP(4,k_fault+1:k-1)';
            Z_GP = [X_Fault_GP(:,k_fault:k-2)', U_Fault_GP(2,k_fault:k-2)'];

            %hyperParameters
            noisevar = 1.0e-07 * [ 0.819920219120628   0.000049368579611   0.222888489932055   0.013943756158247];
            sigma_f = 1.0e+03 *[1.083289684800573   0.009040478282701   0.009999700200313   0.009997077666096];
            sigma_l = 1.0e+03 *[2.595789568810230   3.881067542862866   0.009999004618254   0.010000727078402];
            
            X_estimation = [ard_SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_1, noisevar(1), sigma_f(1), sigma_l(1));...
                ard_SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_2, noisevar(2), sigma_f(2), sigma_l(2));...
                ard_SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_3, noisevar(3), sigma_f(3), sigma_l(3));...
                ard_SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_4, noisevar(4), sigma_f(4), sigma_l(4))];

            error_GP(:,k) = abs(X_estimation-X_Fault_GP(:,k)) ;

%             if  max(error_GP(:,k)) < 0.1 %For testin purpose, the model
%             is always validated
%     %            disp("Validated")

            AB_k = [ard_SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_1, noisevar(1), sigma_f(1), sigma_l(1))';...
                ard_SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_2, noisevar(2), sigma_f(2), sigma_l(2))';...
                ard_SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_3, noisevar(4), sigma_f(3), sigma_l(3))';...
                ard_SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_4, noisevar(4), sigma_f(4), sigma_l(4))'
                ];
            A_k = AB_k(:,1:4);
            B_k = [zeros(4,1), AB_k(:,5)];
            d_k = X_Fault_GP(:,k) - (A_k*X_Fault_GP(:,k-1) + B_k*U_Fault_GP(:,k-1));

            error_mat_GP(k) = norm(A_k-A);


            [U_Fault_GP(:,k),infeasible] = ctrl_broken{X_Fault_GP(:,k),[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k-1)], A_k, B_k, d_k};
            X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k)];
%             else
%                 disp("NOT Validated")
%                 %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
%                 [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k-1)]};
%                 X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k)];
%             end
        end
    end
    
   %% ----------GP with more data---------------------------------
    %At the moment of the fault, we use a GPse linearisation of the
    %identified system, but we use data from when the system was doing ok
    %at the begining
    if do_Fault_GP_moreData ==true
        if k==1
            [U_Fault_GP_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_moreData(:,k),U0};
            X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*U_Fault_GP_moreData(:,k);
        elseif k < k_fault
            [U_Fault_GP_moreData(:,k),infeasible]=ctrl_nominal{X_Fault_GP_moreData(:,k),U_Fault_GP_moreData(:,k-1)};
            X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*U_Fault_GP_moreData(:,k);
        elseif k < k_switch
            [U_Fault_GP_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_moreData(:,k),[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k-1)]};
            X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k)];
        else
            data_number = 100;
            Y_GP_1 = [X_Fault_GP_moreData(1,k_fault+1:k-1)'; Y_1(end-data_number-1+k:end)];
            Y_GP_2 = [X_Fault_GP_moreData(2,k_fault+1:k-1)'; Y_2(end-data_number-1+k:end)];
            Y_GP_3 = [X_Fault_GP_moreData(3,k_fault+1:k-1)'; Y_3(end-data_number-1+k:end)];
            Y_GP_4 = [X_Fault_GP_moreData(4,k_fault+1:k-1)'; Y_4(end-data_number-1+k:end)];
            Z_GP = [X_Fault_GP_moreData(:,k_fault:k-2)', U_Fault_GP_moreData(2,k_fault:k-2)'; Z(end-data_number-1+k:end,:) ];
%             Y_GP_1 = [X_Fault_GP_moreData(1,k_fault+1:k-1)'; Y_1(end-data_number-1+k:end)];
%             Y_GP_2 = [X_Fault_GP_moreData(2,k_fault+1:k-1)'; Y_2(end-data_number-1+k:end)];
%             Y_GP_3 = [X_Fault_GP_moreData(3,k_fault+1:k-1)'; Y_3(end-data_number-1+k:end)];
%             Y_GP_4 = [X_Fault_GP_moreData(4,k_fault+1:k-1)'; Y_4(end-data_number-1+k:end)];
%             Z_GP = [X_Fault_GP_moreData(:,k_fault:k-2)', U_Fault_GP_moreData(2,k_fault:k-2)'; Z( end-data_number-1+k:end,[true,true,true,true,false,true]) ];

            
            
            %hyperParameters
            noisevar = 1.0e-07 * [ 0.819920219120628   0.000049368579611   0.222888489932055   0.013943756158247];
            sigma_f = 1.0e+03 *[1.083289684800573   0.009040478282701   0.009999700200313   0.009997077666096];
            sigma_l = 1.0e+03 *[2.595789568810230   3.881067542862866   0.009999004618254   0.010000727078402];
            error_GP(:,k) = abs(X_estimation-X_Fault_GP(:,k)) ; 
            
            
            
            X_estimation = [ard_SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_1, noisevar(1), sigma_f(1), sigma_l(1));...
                ard_SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_2, noisevar(2), sigma_f(2), sigma_l(2));...
                ard_SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_3, noisevar(3), sigma_f(3), sigma_l(3));...
                ard_SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_4, noisevar(4), sigma_f(4), sigma_l(4))];


            error_GP_moreData(:,k) = abs(X_estimation-X_Fault_GP_moreData(:,k)) ;

%             if  max(error_GP_moreData(:,k)) < 0.1 %Let's assume always
%             validated
%     %            disp("Validated")

            AB_k_moreData = [ard_SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_1, noisevar(1), sigma_f(1), sigma_l(1))';...
                ard_SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_2, noisevar(2), sigma_f(2), sigma_l(2))';...
                ard_SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_3, noisevar(3), sigma_f(3), sigma_l(3))';...
                ard_SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_4, noisevar(4), sigma_f(4), sigma_l(4))'
                ];
            A_k_moreData = AB_k_moreData(:,1:4);
            B_k_moreData = [zeros(4,1), AB_k_moreData(:,5)];
            d_k_moreData = X_Fault_GP_moreData(:,k) - (A_k_moreData*X_Fault_GP_moreData(:,k-1) + B_k_moreData*U_Fault_GP_moreData(:,k-1));

            error_mat_GP_moreData(k) = norm(A_k_moreData-A);


            [U_Fault_GP_moreData(:,k),infeasible] = ctrl_broken{X_Fault_GP_moreData(:,k),[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k-1)], A_k_moreData, B_k_moreData, d_k_moreData};
            X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k)];
%             else
%                 disp("NOT Validated (more data)")
%                 %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
%                 [U_Fault_GP_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_moreData(:,k),[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k-1)]};
%                 X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k)];
%             end
        end
    end     
    
   %% ----------fault with GPse-ard-------------
        %At the moment of the fault, we use a GPse linearisation of the
        %identified systme
    if do_Fault_GPard ==true
        if k==1
            [U_Fault_GPard(:,k),infeasible] = ctrl_nominal{X_Fault_GPard(:,k),U0};
            X_Fault_GPard(:,k+1)=A*X_Fault_GPard(:,k)+B*U_Fault_GPard(:,k);
            
        elseif k < k_fault
            [U_Fault_GPard(:,k),infeasible]=ctrl_nominal{X_Fault_GPard(:,k),U_Fault_GPard(:,k-1)};
            X_Fault_GPard(:,k+1)=A*X_Fault_GPard(:,k)+B*U_Fault_GPard(:,k);
            
        elseif k < k_switch
            [U_Fault_GPard(:,k),infeasible] = ctrl_nominal{X_Fault_GPard(:,k),[U_Fault_GPard(1,k_fault-1);U_Fault_GPard(2,k-1)]};
            X_Fault_GPard(:,k+1)=A*X_Fault_GPard(:,k)+B*[U_Fault_GPard(1,k_fault-1);U_Fault_GPard(2,k)];
        else
            Y_GPard_1 = X_Fault_GPard(1,k_fault+1:k-1)';
            Y_GPard_2 = X_Fault_GPard(2,k_fault+1:k-1)';
            Y_GPard_3 = X_Fault_GPard(3,k_fault+1:k-1)';
            Y_GPard_4 = X_Fault_GPard(4,k_fault+1:k-1)';
            Z_GPard = [X_Fault_GPard(:,k_fault:k-2)', U_Fault_GPard(2,k_fault:k-2)'];

            noisevar =    [0.207011434638091   0.000000000000001   0.000094523372256   0.343759797264307];
            sigma_f = 1.0e+04 *[2.194029301409157   0.000679426763797   0.012952170231516   0.161567940482987];   %size: 1x4
            sigma_l = 1.0e+04 *[0.000043747217450   0.000480105759101   0.000223044213263   0.289137996468431;...
               0.032001534848036   0.012923707339206   0.000932019413174   0.046571952813509;...
               0.882266280783226   0.000801113645661   0.000475033609572   0.000337744140758;...
               1.103686370914297   0.000849685860101   0.000218644426045   0.000178486596270;...
               2.575208850697990   0.000745400599946   0.013843884467876   0.033374538994487]; %size: 5x4

                   
            X_estimation = [ard_SE_mean([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_1, noisevar(1), sigma_f(1), sigma_l(:,1));...
                ard_SE_mean([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_2, noisevar(2), sigma_f(2), sigma_l(:,2));...
                ard_SE_mean([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_3, noisevar(3), sigma_f(3), sigma_l(:,3));...
                ard_SE_mean([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_4, noisevar(4), sigma_f(4), sigma_l(:,4))];


            error_GPard(:,k) = abs(X_estimation-X_Fault_GPard(:,k)) ;

%             if  max(error_GPard(:,k)) < 0.1 %considerated as always
%             validated
    %            disp("Validated")

            AB_k = [ard_SE_deriv([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_1, noisevar(1), sigma_f(1), sigma_l(:,1))';...
                ard_SE_deriv([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_2, noisevar(2), sigma_f(2), sigma_l(:,2))';...
                ard_SE_deriv([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_3, noisevar(3), sigma_f(3), sigma_l(:,3))';...
                ard_SE_deriv([X_Fault_GPard(:,k-1)', U_Fault_GPard(2,k-1)], Z_GPard, Y_GPard_4, noisevar(4), sigma_f(4), sigma_l(:,4))'
                ];
            A_k = AB_k(:,1:4);
            B_k = [zeros(4,1), AB_k(:,5)];
            d_k = X_Fault_GPard(:,k) - (A_k*X_Fault_GPard(:,k-1) + B_k*U_Fault_GPard(:,k-1));

            error_mat_GPard(k) = norm(A_k-A);


            [U_Fault_GPard(:,k),infeasible] = ctrl_broken{X_Fault_GPard(:,k),[U_Fault_GPard(1,k_fault-1);U_Fault_GPard(2,k-1)], A_k, B_k, d_k};
            X_Fault_GPard(:,k+1)=A*X_Fault_GPard(:,k)+B*[U_Fault_GPard(1,k_fault-1);U_Fault_GPard(2,k)];
%             else
%                 disp("NOT Validated ard")
%                 %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
%                 [U_Fault_GPard(:,k),infeasible] = ctrl_nominal{X_Fault_GPard(:,k),[U_Fault_GPard(1,k_fault-1);U_Fault_GPard(2,k-1)]};
%                 X_Fault_GPard(:,k+1)=A*X_Fault_GPard(:,k)+B*[U_Fault_GPard(1,k_fault-1);U_Fault_GPard(2,k)];
%             end
        end
    end
    
   %% ----------GP seard with more data---------------------------------
    %At the moment of the fault, we use a GPse linearisation of the
    %identified system, but we use data from when the system was doing ok
    %at the begining
    if do_Fault_GPard_moreData ==true
        if k==1
            [U_Fault_GPard_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GPard_moreData(:,k),U0};
            X_Fault_GPard_moreData(:,k+1)=A*X_Fault_GPard_moreData(:,k)+B*U_Fault_GPard_moreData(:,k);
        elseif k < k_fault
            [U_Fault_GPard_moreData(:,k),infeasible]=ctrl_nominal{X_Fault_GPard_moreData(:,k),U_Fault_GPard_moreData(:,k-1)};
            X_Fault_GPard_moreData(:,k+1)=A*X_Fault_GPard_moreData(:,k)+B*U_Fault_GPard_moreData(:,k);
        elseif k < k_switch
            [U_Fault_GPard_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GPard_moreData(:,k),[U_Fault_GPard_moreData(1,k_fault-1);U_Fault_GPard_moreData(2,k-1)]};
            X_Fault_GPard_moreData(:,k+1)=A*X_Fault_GPard_moreData(:,k)+B*[U_Fault_GPard_moreData(1,k_fault-1);U_Fault_GPard_moreData(2,k)];
        else
            data_number = 100;
            Y_GPard_1 = [X_Fault_GPard_moreData(1,k_fault+1:k-1)'; Y_1(end-data_number-1+k:end)]; %max(2,k-data_number)
            Y_GPard_2 = [X_Fault_GPard_moreData(2,k_fault+1:k-1)'; Y_2(end-data_number-1+k:end)];
            Y_GPard_3 = [X_Fault_GPard_moreData(3,k_fault+1:k-1)'; Y_3(end-data_number-1+k:end)];
            Y_GPard_4 = [X_Fault_GPard_moreData(4,k_fault+1:k-1)'; Y_4(end-data_number-1+k:end)];
            Z_GPard = [X_Fault_GPard_moreData(:,k_fault:k-2)', U_Fault_GPard_moreData(2,k_fault:k-2)'; Z(end-data_number-1+k:end,:) ];

            
            
            noisevar =    [0.207011434638091   0.000000000000001   0.000094523372256   0.343759797264307];
            sigma_f = 1.0e+04 *[2.194029301409157   0.000679426763797   0.012952170231516   0.161567940482987];   %size: 1x4
            sigma_l = 1.0e+04 *[0.000043747217450   0.000480105759101   0.000223044213263   0.289137996468431;...
                   0.032001534848036   0.012923707339206   0.000932019413174   0.046571952813509;...
                   0.882266280783226   0.000801113645661   0.000475033609572   0.000337744140758;...
                   1.103686370914297   0.000849685860101   0.000218644426045   0.000178486596270;...
                   2.575208850697990   0.000745400599946   0.013843884467876   0.033374538994487]; %size: 5x4  


            X_estimation = [ard_SE_mean([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_1, noisevar(1), sigma_f(1), sigma_l(:,1));...
                ard_SE_mean([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_2, noisevar(2), sigma_f(2), sigma_l(:,2));...
                ard_SE_mean([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_3, noisevar(3), sigma_f(3), sigma_l(:,3));...
                ard_SE_mean([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_4, noisevar(4), sigma_f(4), sigma_l(:,4))];


            error_GPard_moreData(:,k) = abs(X_estimation-X_Fault_GPard_moreData(:,k)) ;

%             if  max(error_GPard_moreData(:,k)) < 0.1 %Considereded always
%             validated
%     %            disp("Validated")

            AB_k_moreData = [ard_SE_deriv([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_1, noisevar(1), sigma_f(1), sigma_l(:,1))';...
                ard_SE_deriv([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_2, noisevar(2), sigma_f(2), sigma_l(:,2))';...
                ard_SE_deriv([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_3, noisevar(3), sigma_f(3), sigma_l(:,3))';...
                ard_SE_deriv([X_Fault_GPard_moreData(:,k-1)', U_Fault_GPard_moreData(2,k-1)], Z_GPard, Y_GPard_4, noisevar(4), sigma_f(4), sigma_l(:,4))'
                ];
            A_k_moreData = AB_k_moreData(:,1:4);
            B_k_moreData = [zeros(4,1), AB_k_moreData(:,5)];
            d_k_moreData = X_Fault_GPard_moreData(:,k) - (A_k_moreData*X_Fault_GPard_moreData(:,k-1) + B_k_moreData*U_Fault_GPard_moreData(:,k-1));

            error_mat_GPard_moreData(k) = norm(A_k_moreData-A);


            [U_Fault_GPard_moreData(:,k),infeasible] = ctrl_broken{X_Fault_GPard_moreData(:,k),[U_Fault_GPard_moreData(1,k_fault-1);U_Fault_GPard_moreData(2,k-1)], A_k_moreData, B_k_moreData, d_k_moreData};
            X_Fault_GPard_moreData(:,k+1)=A*X_Fault_GPard_moreData(:,k)+B*[U_Fault_GPard_moreData(1,k_fault-1);U_Fault_GPard_moreData(2,k)];
%             else
%                 disp("NOT Validated (more data ard)")
%                 %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
%                 [U_Fault_GPard_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GPard_moreData(:,k),[U_Fault_GPard_moreData(1,k_fault-1);U_Fault_GPard_moreData(2,k-1)]};
%                 X_Fault_GPard_moreData(:,k+1)=A*X_Fault_GPard_moreData(:,k)+B*[U_Fault_GPard_moreData(1,k_fault-1);U_Fault_GPard_moreData(2,k)];
%             end
        end
    end
        
   %% -------fault with L2NW and GP---------------------------------------
    %Identifying the system using L2NW, se kernel
    if do_L2NW_GP == true
        if k==1
            [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),U0};
            X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*U_Fault_L2NWGP(:,k);
        elseif k < k_fault
            [U_Fault_L2NWGP(:,k),infeasible]=ctrl_nominal{X_Fault_L2NWGP(:,k),U_Fault_L2NWGP(:,k-1)};
            X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*U_Fault_L2NWGP(:,k);
        elseif k < k_switch
            [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k-1)]};
            X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
        else
 
             lambda =      1.0e+02 *[0.000003423669559   0.000000000000000   3.655141778377240   0.000000000005348];
             sigma_f2 = 1.0e+03 *[1.582306003922267   1.582428374667688   1.581677379931159   1.582324525677518];
             sigma_l2 = [0.351945588681686   0.212271192005340   0.490086358910314   0.241244856635137];

            param = [lambda(1), sigma_f2(1), sigma_l2(1);...
                lambda(2), sigma_f2(2), sigma_l2(2);...
                lambda(3), sigma_f2(3), sigma_l2(3);...
                lambda(4), sigma_f2(4), sigma_l2(4)];


            data_number_l2 = 100;
            Y_L2NW_1GP = [X_Fault_L2NWGP(1,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_2GP = [X_Fault_L2NWGP(2,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_3GP = [X_Fault_L2NWGP(3,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_4GP = [X_Fault_L2NWGP(4,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Z_L2NW_GP = [X_Fault_L2NWGP(:,k_fault:k-2)', U_Fault_L2NWGP(2,k_fault:k-2)'];%;  Z( end -data_number_l2-1+k:end,:)];
            


            X_estimation_L2NWGP = [ard_L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)], Z_L2NW_GP, Y_L2NW_1GP, param(1,:));...
                ard_L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)], Z_L2NW_GP, Y_L2NW_2GP, param(2,:));...
                ard_L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)], Z_L2NW_GP, Y_L2NW_3GP, param(3,:));...
                ard_L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)], Z_L2NW_GP, Y_L2NW_4GP, param(4,:))...
                ];


            error_L2NWGP(:,k) = abs(X_Fault_L2NWGP(:,k) -  X_estimation_L2NWGP);

%             if  max(error_L2NWGP(:,k)) < 300 %considered always validate
            z_old = z_new;
            options = optimoptions('fmincon','MaxFunctionEvaluations',9000);
            z_new = fmincon(@(z) ard_MPC_cost(z, horizon, Q_huge, R_huge, R_delta_huge),...  %@(z) MPC_cost(z, horizon, Q_huge, R_huge, R_delta_huge),
                z_old , C_ineq , c_ineq , [], [], [], [],...
                @(z) ard_MPC_constraitn(z, horizon, X_Fault_L2NWGP(:,k), U_Fault_L2NWGP(2,k-1),...
            Z_L2NW_GP, Y_L2NW_1GP, Y_L2NW_2GP, Y_L2NW_3GP, Y_L2NW_4GP, lambda, sigma_f2, sigma_l2),...
            options); 

            U_Fault_L2NWGP(2,k) = z_new(2,horizon+1);
            U_Fault_L2NWGP(1,k) = U_Fault_L2NWGP(1,k-1);        %there is a fault
            %end FminconProblem



            X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
%             else
%                 disp('L2NW2  NOT validated')
%                 [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k-1)]};
%                 X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
%             end
        end
    end
    
   %% -------fault with L2NW and GPard---------------------------------------
    %Identifying the system using L2NW, se kernel
    if do_L2NW_GPard == true
        if k==1
            [U_Fault_L2NWGP_ard(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP_ard(:,k),U0};
            X_Fault_L2NWGP_ard(:,k+1)=A*X_Fault_L2NWGP_ard(:,k)+B*U_Fault_L2NWGP_ard(:,k);
        elseif k < k_fault
            [U_Fault_L2NWGP_ard(:,k),infeasible]=ctrl_nominal{X_Fault_L2NWGP_ard(:,k),U_Fault_L2NWGP_ard(:,k-1)};
            X_Fault_L2NWGP_ard(:,k+1)=A*X_Fault_L2NWGP_ard(:,k)+B*U_Fault_L2NWGP_ard(:,k);
        elseif k < k_switch
            [U_Fault_L2NWGP_ard(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP_ard(:,k),[U_Fault_L2NWGP_ard(1,k_fault-1);U_Fault_L2NWGP_ard(2,k-1)]};
            X_Fault_L2NWGP_ard(:,k+1)=A*X_Fault_L2NWGP_ard(:,k)+B*[U_Fault_L2NWGP_ard(1,k_fault-1);U_Fault_L2NWGP_ard(2,k)];
        else

            lambda =     1.0e+03 *[1.455543128007772   0.000000000000001
                0.865249329706008   0.007700284142739];
            sigma_f2 =     1.0e+03 *[ 1.584160310465225   1.582760479010113   1.582919812162524   1.582541347721900]; %size: 1x4
            sigma_l2 =  1.0e+03 *[1.454206567885842   0.000331156020930   0.590599646645230   0.006098649717990;...
               0.000964946692309   0.000226992358221   0.000561657338271   0.000472731065963;...
               0.175444781224309   0.007867196986228   0.044388834761815  -0.002571193069446;...
               0.000057030705522  -0.000155566993686   0.000106365542096  -0.000063536431170;...
              -0.000040030459786   0.001292079840159   0.000035331193726  -0.000035900590580]; %size: 5x4
       

            param = [lambda(1), sigma_f2(1), sigma_l2(:,1)';...
                lambda(2), sigma_f2(2), sigma_l2(:,2)';...
                lambda(3), sigma_f2(3), sigma_l2(:,3)';...
                lambda(4), sigma_f2(4), sigma_l2(:,4)'];


            data_number_l2 = 100;
            Y_L2NW_1GP_ard = [X_Fault_L2NWGP_ard(1,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_2GP_ard = [X_Fault_L2NWGP_ard(2,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_3GP_ard = [X_Fault_L2NWGP_ard(3,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_4GP_ard = [X_Fault_L2NWGP_ard(4,k_fault+1:k-1)'];%; Y_1(end -data_number_l2-1+k:end)];
            Z_L2NW_GP_ard = [X_Fault_L2NWGP_ard(:,k_fault:k-2)', U_Fault_L2NWGP_ard(2,k_fault:k-2)'];%;  Z( end -data_number_l2-1+k:end,:)];
            


            X_estimation_L2NWGP_ard = [ard_L2NW_mean([X_Fault_L2NWGP_ard(:,k-1)', U_Fault_L2NWGP_ard(2,k-1)],Z_L2NW_GP_ard,Y_L2NW_1GP_ard,param(1,:));...
                ard_L2NW_mean([X_Fault_L2NWGP_ard(:,k-1)', U_Fault_L2NWGP_ard(2,k-1)],Z_L2NW_GP_ard,Y_L2NW_2GP_ard,param(2,:));...
                ard_L2NW_mean([X_Fault_L2NWGP_ard(:,k-1)', U_Fault_L2NWGP_ard(2,k-1)],Z_L2NW_GP_ard,Y_L2NW_3GP_ard,param(3,:));...
                ard_L2NW_mean([X_Fault_L2NWGP_ard(:,k-1)', U_Fault_L2NWGP_ard(2,k-1)],Z_L2NW_GP_ard,Y_L2NW_4GP_ard,param(4,:))...
                ];


            error_L2NWGP_ard(:,k) = abs(X_Fault_L2NWGP_ard(:,k) -  X_estimation_L2NWGP_ard);

            %if  (max(error_L2NWGP_ard(:,k)) < 300) %considered as always
            %validated
                
                %define Fmincon problem
          z_old = z_new_ard;
            options = optimoptions('fmincon','MaxFunctionEvaluations',30000);
            tic
            z_new_ard = fmincon(@(z) ard_MPC_cost(z, horizon, Q_huge, R_huge, R_delta_huge),...  %@(z) MPC_cost(z, horizon, Q_huge, R_huge, R_delta_huge),
                z_old , C_ineq , c_ineq , [], [], [], [],...
                @(z) ard_MPC_constraitn(z, horizon, X_Fault_L2NWGP_ard(:,k), U_Fault_L2NWGP_ard(2,k-1),...
                Z_L2NW_GP_ard, Y_L2NW_1GP_ard, Y_L2NW_2GP_ard, Y_L2NW_3GP_ard, Y_L2NW_4GP_ard,...
                lambda, sigma_f2, sigma_l2), options ); 
            toc
            
            U_Fault_L2NWGP_ard(2,k) = z_new_ard(2,horizon+1);
            U_Fault_L2NWGP_ard(1,k) = U_Fault_L2NWGP_ard(1,k-1);        %there is a fault
            %end FminconProblem



            X_Fault_L2NWGP_ard(:,k+1)=A*X_Fault_L2NWGP_ard(:,k)+B*[U_Fault_L2NWGP_ard(1,k_fault-1);U_Fault_L2NWGP_ard(2,k)];
%             
%             else
%                 disp('L2NW2  NOT validated')
%                 [U_Fault_L2NWGP_ard(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP_ard(:,k),[U_Fault_L2NWGP_ard(1,k_fault-1);U_Fault_L2NWGP_ard(2,k-1)]};
%                 X_Fault_L2NWGP_ard(:,k+1)=A*X_Fault_L2NWGP_ard(:,k)+B*[U_Fault_L2NWGP_ard(1,k_fault-1);U_Fault_L2NWGP_ard(2,k)];
% 
% 
%             end
        end
    end
    
   %% ---------------ALL----------------------------------
    if mod(k,5) == 0
        save("result_29_L2NWard_bis","X_Fault_L2NWGP_ard","X_estimation_L2NWGP_ard", "U_Fault_L2NWGP_ard", "error_L2NWGP_ard")
        disp(k)
    end
    
    
end    
    

%% Plotting ================================================================
t=Ts*(0:stepNumber-1);
%% States
figure
hold on
leg = [];

subplot(4,1,1) %State 1
hold on
title('Pitch rate [rad/s]')
xlabel('time s')
if do_noFault == true  
    plot(t,X_noFault(1,:),'k-.')
end
if do_Fault == true
    plot(t,X_Fault(1,:),'k:')
end
if do_Fault_perfect == true
    plot(t,X_Fault_Perfect(1,:),'k-','LineWidth',2)
end
if do_Fault_GP == true
    plot(t,X_Fault_GP(1,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t,X_Fault_GP_moreData(1,:),'b-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(1,:),'b:')
end
if do_Fault_GPard == true
    plot(t,X_Fault_GPard(1,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t,X_Fault_GPard_moreData(1,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t,X_Fault_L2NWGP_ard(1,:),'r:')
end

subplot(4,1,2) %State 2-----------------
hold on
title('Air speed [m/s]')
xlabel('time s')
if do_noFault == true  
    plot(t,X_noFault(2,:),'k-.')
end
if do_Fault == true
    plot(t,X_Fault(2,:),'k:')
end
if do_Fault_perfect == true
    plot(t,X_Fault_Perfect(2,:),'k-','LineWidth',2)
end
if do_Fault_GP == true
    plot(t,X_Fault_GP(2,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t,X_Fault_GP_moreData(2,:),'b-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(2,:),'b:')
end
if do_Fault_GPard == true
    plot(t,X_Fault_GPard(2,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t,X_Fault_GPard_moreData(2,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t,X_Fault_L2NWGP_ard(2,:),'r:')
end


subplot(4,1,3) %State 3---------------
hold on
title('Pitch angle [rad]')
xlabel('time s')
if do_noFault == true  
    plot(t,X_noFault(3,:),'k-.')
end
if do_Fault == true
    plot(t,X_Fault(3,:),'k:')
end
if do_Fault_perfect == true
    plot(t,X_Fault_Perfect(3,:),'k-','LineWidth',2)
end
if do_Fault_GP == true
    plot(t,X_Fault_GP(3,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t,X_Fault_GP_moreData(3,:),'b-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(3,:),'b:')
end
if do_Fault_GPard == true
    plot(t,X_Fault_GPard(3,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t,X_Fault_GPard_moreData(3,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t,X_Fault_L2NWGP_ard(3,:),'r:')
end

    
subplot(4,1,4) %State 4---------
hold on
title('Angle of attack [rad]')
xlabel('time s')
if do_noFault == true  
    plot(t,X_noFault(4,:),'k-.')
end
if do_Fault == true
    plot(t,X_Fault(4,:),'k:')
end
if do_Fault_perfect == true
    plot(t,X_Fault_Perfect(4,:),'k-','LineWidth',2)
end
if do_Fault_GP == true
    plot(t,X_Fault_GP(4,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t,X_Fault_GP_moreData(4,:),'b-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(4,:),'b:')
end
if do_Fault_GPard == true
    plot(t,X_Fault_GPard(4,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t,X_Fault_GPard_moreData(4,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t,X_Fault_L2NWGP_ard(4,:),'r:')
end


%legend:
if do_noFault == true  
    leg = [leg, "no Fault scenario"];
end
if do_Fault == true
    leg = [leg, "Fault with nominal model"];
end
if do_Fault_perfect == true
    leg = [leg, "Fault with perfect model"];
end
if do_Fault_GP == true
    leg = [leg, "Fault with GP model"];    
end
if do_Fault_GP_moreData == true
    leg = [leg, "Fault with GP model (with more data)"];
end
if do_L2NW_GP == true
    leg = [leg, "Fault with gaussian L2NW model"];
end
if do_Fault_GPard == true
    leg = [leg, "Fault with GPard model"];    
end
if do_Fault_GPard_moreData == true
    leg = [leg, "Fault with GPard model (with more data)"];
end
if do_L2NW_GPard == true
    leg = [leg, "Fault with gaussian seard L2NW model"];
end
legend(leg)
%% error --------------------------------------------------------------------------
figure
hold on
leg = [];

subplot(4,1,1) %State 1
hold on
title('Error on pitch rate [rad/s]')
xlabel('time s')
ylim([0,0.00001])
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(1,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(1,:),'b-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(1,:),'b:')
end
if do_Fault_GPard == true
    plot(t(1:end-1),error_GPard(1,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t(1:end-1),error_GPard_moreData(1,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t(1:end-1),error_L2NWGP_ard(1,:),'r:')
end

subplot(4,1,2) %State 2-----------------
hold on
title('Error on air speed [m/s]')
xlabel('time s')
ylim([0,0.0001])
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(2,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(2,:),'b-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(2,:),'b:')
end
if do_Fault_GPard == true
    plot(t(1:end-1),error_GPard(2,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t(1:end-1),error_GPard_moreData(2,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t(1:end-1),error_L2NWGP_ard(2,:),'r:')
end

subplot(4,1,3) %State 3---------------
hold on
title('Error on pitch angle [rad]')
xlabel('time s')
ylim([0,0.001])
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(3,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(3,:),'b-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(3,:),'b:')
end
if do_Fault_GPard == true
    plot(t(1:end-1),error_GPard(3,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t(1:end-1),error_GPard_moreData(3,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t(1:end-1),error_L2NWGP_ard(3,:),'r:')
end
    
subplot(4,1,4) %State 4---------
hold on
title('Error on angle of attack [rad]')
xlabel('time s')
ylim([0,0.001])
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(4,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(4,:),'b-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(4,:),'b:')
end
if do_Fault_GPard == true
    plot(t(1:end-1),error_GPard(4,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t(1:end-1),error_GPard_moreData(4,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t(1:end-1),error_L2NWGP_ard(4,:),'r:')
end

%legend:
leg = [];
if do_Fault_GP == true
    leg = [leg, "Fault with GP model"];    
end
if do_Fault_GP_moreData == true
    leg = [leg, "Fault with GP model (with more data)"];
end
if do_L2NW_GP == true
    leg = [leg, "Fault with gaussian L2NW model"];
end
if do_Fault_GPard == true
    leg = [leg, "Fault with GPard model"];    
end
if do_Fault_GPard_moreData == true
    leg = [leg, "Fault with GPard model (with more data)"];
end
if do_L2NW_GPard == true
    leg = [leg, "Fault with gaussian seard L2NW model"];
end
legend(leg)

%% Inputs
figure
hold on
leg = [];

subplot(1,2,1) %State 1
hold on
title('Elevator')
xlabel('time s')
if do_noFault == true  
    plot(t(1:end-1),U_noFault(1,:),'k-.')
end
if do_Fault == true
    plot(t(1:end-1),U_Fault(1,:),'k:')
end
if do_Fault_perfect == true
    plot(t(1:end-1),U_Fault_Perfect(1,:),'k-','LineWidth',2)
end
if do_Fault_GP == true
    plot(t(1:end-1),U_Fault_GP(1,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),U_Fault_GP_moreData(1,:),'b-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),U_Fault_L2NWGP(1,:),'b:')
end
if do_Fault_GPard == true
    plot(t(1:end-1),U_Fault_GPard(1,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t(1:end-1),U_Fault_GPard_moreData(1,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t(1:end-1),U_Fault_L2NWGP_ard(1,:),'r:')
end

subplot(1,2,2) %State 2-----------------
hold on
title('THS')
xlabel('time s')
if do_noFault == true  
    plot(t(1:end-1),U_noFault(2,:),'k-.')
end
if do_Fault == true
    plot(t(1:end-1),U_Fault(2,:),'k:')
end
if do_Fault_perfect == true
    plot(t(1:end-1),U_Fault_Perfect(2,:),'k-','LineWidth',2)
end
if do_Fault_GP == true
    plot(t(1:end-1),U_Fault_GP(2,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),U_Fault_GP_moreData(2,:),'b-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),U_Fault_L2NWGP(2,:),'b:')
end
if do_Fault_GPard == true
    plot(t(1:end-1),U_Fault_GPard(2,:),'r-')    
end
if do_Fault_GPard_moreData == true
    plot(t(1:end-1),U_Fault_GPard_moreData(2,:),'r-.')
end
if do_L2NW_GPard == true
    plot(t(1:end-1),U_Fault_L2NWGP_ard(2,:),'r:')
end


%legend:
if do_noFault == true  
    leg = [leg, "no Fault scenario"];
end
if do_Fault == true
    leg = [leg, "Fault with nominal model"];
end
if do_Fault_perfect == true
    leg = [leg, "Fault with perfect model"];
end
if do_Fault_GP == true
    leg = [leg, "Fault with GP model"];    
end
if do_Fault_GP_moreData == true
    leg = [leg, "Fault with GP model (with more data)"];
end
if do_L2NW_GP == true
    leg = [leg, "Fault with gaussian L2NW model"];
end
if do_Fault_GPard == true
    leg = [leg, "Fault with GPard model"];    
end
if do_Fault_GPard_moreData == true
    leg = [leg, "Fault with GPard model (with more data)"];
end
if do_L2NW_GPard == true
    leg = [leg, "Fault with gaussian seard L2NW model"];
end
legend(leg)

