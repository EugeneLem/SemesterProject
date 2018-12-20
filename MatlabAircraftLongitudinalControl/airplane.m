format long
load('Data_noFault.mat') %load data from when there are no test
Z = Z(:, [1,2,3,4,6]);   %For faulty actuator, we do not want the first input   

import casadi.*
opti = casadi.Opti();   


do_noFault = true;                 %list of calculation we want to do:
do_Fault = true;
do_Fault_perfect = true;
do_Fault_GP = true;
do_Fault_GP_moreData = true;
do_Fault_GP_allData = true;
do_L2NW_1 = true;
do_L2NW_GP = true;



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
R_huge = blkdiag(kron(eye(horizon-1), R(2,2)),P_k(6,6));

%constraint
C_u = [1,0 ; 0,1 ; -1,0 ; 0,-1];        %Constraint on input
c_u = [15;4.6;25;10.4]; 
C_u_delta = [-1 0   1 0;...
            0 -1  0 1;...
            1  0 -1 0;...
            0  1  0 -1];
c_u_delta = [37;0.236;37;0.236]*Ts;% constraint on input change rate


C_u_redu = [ 1 ; -1];
c_u_redu = [4.6;10.4];
C_u_delta_redu = [-1  1; 1  -1];
c_u_delta_redu = [0.236; 0.236]*Ts;

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



%%
%% Simulation constant =====================================================
stepNumber = 20;       %HERE
k_fault=2;               %We initialize everything to zero
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
error_mat_GP = zeros(1,stepNumber-1);

X_Fault_GP_moreData = zeros(4,stepNumber); 
U_Fault_GP_moreData = zeros(2,stepNumber-1); 
X_Fault_GP_moreData(:,1) = X0;
u_Fault_GP_moreData(:,1) = U0;
error_GP_moreData = zeros(4,stepNumber-1);
error_mat_GP_moreData = zeros(1,stepNumber-1);

X_Fault_GP_allData = zeros(4,stepNumber); 
U_Fault_GP_allData = zeros(2,stepNumber-1); 
X_Fault_GP_allData(:,1) = X0;
u_Fault_GP_allData(:,1) = U0;
error_GP_allData = zeros(4,stepNumber-1);
error_mat_GP_allData = zeros(1,stepNumber-1);

X_Fault_L2NW = zeros(4,stepNumber); 
U_Fault_L2NW = zeros(2,stepNumber-1); 
X_Fault_L2NW(:,1) = X0;
u_Fault_L2NW(:,1) = U0;
error_L2NW = zeros(4,stepNumber-1);
error_mat_L2NW = zeros(1,stepNumber-1);


X_Fault_L2NWGP = zeros(4,stepNumber); 
U_Fault_L2NWGP = zeros(2,stepNumber-1); 
X_Fault_L2NWGP(:,1) = X0;
u_Fault_L2NWGP(:,1) = U0;
error_L2NWGP = zeros(4,stepNumber-1);

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

    %         noisevar = 0.000001;
    %         sigma_f = 2.17*10^(3)*ones(1,4);
    %         sigma_l = 1.14*10^(4)*ones(1,4);
    %         sigma_f = [1.677582380126192, 80.087113811947248, 1.371605612665715, 6.573873522600369];
    %         sigma_l = [10.695781461698118, 99.935986868827754,  9.046354763354064, 32.505687565655741];
            noisevar = 0.00000001;
%             sigma_f =[  20.575591520093276,  17.878167302000445,  16.199466399148097,  17.938357317647856];
%             sigma_l =[  99.966954902978543,  34.241347755939984,  78.546996930435938,  88.147249188528974];
            sigma_f =   1.0e+02 *[ 2.639696625733967   0.170294093647814   0.085880808896654   0.101553751261396];
            sigma_l =[ 99.991949107543803  25.640539532091616  34.284915666906144  25.862479178807501];
            %         
            X_estimation = [SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1));...
                SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2));...
                SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3));...
                SE_mean([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))];


            error_GP(:,k) = abs(X_estimation-X_Fault_GP(:,k)) ;

            if  max(error_GP(:,k)) < 0.1
    %            disp("Validated")

                AB_k = [SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1))';...
                    SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2))';...
                    SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3))';...
                    SE_deriv([X_Fault_GP(:,k-1)', U_Fault_GP(2,k-1)], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))'
                    ];
                A_k = AB_k(:,1:4);
                B_k = [zeros(4,1), AB_k(:,5)];
                d_k = X_Fault_GP(:,k) - (A_k*X_Fault_GP(:,k-1) + B_k*U_Fault_GP(:,k-1));

                error_mat_GP(k) = norm(A_k-A);


                [U_Fault_GP(:,k),infeasible] = ctrl_broken{X_Fault_GP(:,k),[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k-1)], A_k, B_k, d_k};
                X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k)];
            else
                disp("NOT Validated")
                %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
                [U_Fault_GP(:,k),infeasible] = ctrl_nominal{X_Fault_GP(:,k),[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k-1)]};
                X_Fault_GP(:,k+1)=A*X_Fault_GP(:,k)+B*[U_Fault_GP(1,k_fault-1);U_Fault_GP(2,k)];
            end
        end
    end
    
    %% --------------------GP with more data---------------------------------
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
            Y_GP_1 = [X_Fault_GP_moreData(1,max(2,k-data_number):k-1)'; Y_1(end-data_number-1+k:end)];
            Y_GP_2 = [X_Fault_GP_moreData(2,max(2,k-data_number):k-1)'; Y_2(end-data_number-1+k:end)];
            Y_GP_3 = [X_Fault_GP_moreData(3,max(2,k-data_number):k-1)'; Y_3(end-data_number-1+k:end)];
            Y_GP_4 = [X_Fault_GP_moreData(4,max(2,k-data_number):k-1)'; Y_4(end-data_number-1+k:end)];
            Z_GP = [X_Fault_GP_moreData(:,max(1,k-data_number-1):k-2)', U_Fault_GP_moreData(2,max(1,k-data_number-1):k-2)'; Z(end-data_number-1+k:end,:) ];
%             Y_GP_1 = [X_Fault_GP_moreData(1,k_fault+1:k-1)'; Y_1(end-data_number-1+k:end)];
%             Y_GP_2 = [X_Fault_GP_moreData(2,k_fault+1:k-1)'; Y_2(end-data_number-1+k:end)];
%             Y_GP_3 = [X_Fault_GP_moreData(3,k_fault+1:k-1)'; Y_3(end-data_number-1+k:end)];
%             Y_GP_4 = [X_Fault_GP_moreData(4,k_fault+1:k-1)'; Y_4(end-data_number-1+k:end)];
%             Z_GP = [X_Fault_GP_moreData(:,k_fault:k-2)', U_Fault_GP_moreData(2,k_fault:k-2)'; Z( end-data_number-1+k:end,[true,true,true,true,false,true]) ];

            
            
    %         noisevar = 0.000001;
    %         sigma_f = 2.17*10^(3)*ones(1,4);
    %         sigma_l = 1.14*10^(4)*ones(1,4);
    %         sigma_f = [1.677582380126192, 80.087113811947248, 1.371605612665715, 6.573873522600369];
    %         sigma_l = [10.695781461698118, 99.935986868827754,  9.046354763354064, 32.505687565655741];
            noisevar = 0.00000001;
%             sigma_f =[  20.575591520093276,  17.878167302000445,  16.199466399148097,  17.938357317647856];
%             sigma_l =[  99.966954902978543,  34.241347755939984,  78.546996930435938,  88.147249188528974];
            sigma_f =   1.0e+02 *[ 2.639696625733967   0.170294093647814   0.085880808896654   0.101553751261396];
            sigma_l =[ 99.991949107543803  25.640539532091616  34.284915666906144  25.862479178807501];
            %         
            X_estimation = [SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1));...
                SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2));...
                SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3));...
                SE_mean([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))];


            error_GP_moreData(:,k) = abs(X_estimation-X_Fault_GP_moreData(:,k)) ;

            if  max(error_GP_moreData(:,k)) < 0.1
    %            disp("Validated")

                AB_k_moreData = [SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1))';...
                    SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2))';...
                    SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3))';...
                    SE_deriv([X_Fault_GP_moreData(:,k-1)', U_Fault_GP_moreData(2,k-1)], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))'
                    ];
                A_k_moreData = AB_k_moreData(:,1:4);
                B_k_moreData = [zeros(4,1), AB_k_moreData(:,5)];
                d_k_moreData = X_Fault_GP_moreData(:,k) - (A_k_moreData*X_Fault_GP_moreData(:,k-1) + B_k_moreData*U_Fault_GP_moreData(:,k-1));

                error_mat_GP_moreData(k) = norm(A_k_moreData-A);


                [U_Fault_GP_moreData(:,k),infeasible] = ctrl_broken{X_Fault_GP_moreData(:,k),[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k-1)], A_k_moreData, B_k_moreData, d_k_moreData};
                X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k)];
            else
                disp("NOT Validated (more data)")
                %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
                [U_Fault_GP_moreData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_moreData(:,k),[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k-1)]};
                X_Fault_GP_moreData(:,k+1)=A*X_Fault_GP_moreData(:,k)+B*[U_Fault_GP_moreData(1,k_fault-1);U_Fault_GP_moreData(2,k)];
            end
        end
    end
        
        
    %% --------------------GP with all data---------------------------------
    %At the moment of the fault, we use a GPse linearisation of the
    %identified system, but we use data from when the system was doing ok
    %at the begining, WITHOUT discardind it... is shit

    if do_Fault_GP_allData == true
        if k==1
            [U_Fault_GP_allData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_allData(:,k),U0};
            X_Fault_GP_allData(:,k+1)=A*X_Fault_GP_allData(:,k)+B*U_Fault_GP_allData(:,k);
        elseif k < k_fault
            [U_Fault_GP_allData(:,k),infeasible]=ctrl_nominal{X_Fault_GP_allData(:,k),U_Fault_GP_allData(:,k-1)};
            X_Fault_GP_allData(:,k+1)=A*X_Fault_GP_allData(:,k)+B*U_Fault_GP_allData(:,k);
        elseif k < k_switch
            [U_Fault_GP_allData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_allData(:,k),[U_Fault_GP_allData(1,k_fault-1);U_Fault_GP_allData(2,k-1)]};
            X_Fault_GP_allData(:,k+1)=A*X_Fault_GP_allData(:,k)+B*[U_Fault_GP_allData(1,k_fault-1);U_Fault_GP_allData(2,k)];
        else
            Y_GP_1 = [X_Fault_GP_allData(1,2:k-1)'; Y_1(1:end)];
            Y_GP_2 = [X_Fault_GP_allData(2,2:k-1)'; Y_2(1:end)];
            Y_GP_3 = [X_Fault_GP_allData(3,2:k-1)'; Y_3(1:end)];
            Y_GP_4 = [X_Fault_GP_allData(4,2:k-1)'; Y_4(1:end)];
            Z_GP = [X_Fault_GP_allData(:,1:k-2)', U_Fault_GP_allData(2,1:k-2)'; Z ];

    %         noisevar = 0.000001;
    %         sigma_f = 2.17*10^(3)*ones(1,4);
    %         sigma_l = 1.14*10^(4)*ones(1,4);
    %         sigma_f = [1.677582380126192, 80.087113811947248, 1.371605612665715, 6.573873522600369];
    %         sigma_l = [10.695781461698118, 99.935986868827754,  9.046354763354064, 32.505687565655741];
            noisevar = 0.00000001;
            sigma_f =[  20.575591520093276,  17.878167302000445,  16.199466399148097,  17.938357317647856];
            sigma_l =[  99.966954902978543,  34.241347755939984,  78.546996930435938,  88.147249188528974];
            %         
            X_estimation = [SE_mean([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1));...
                SE_mean([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2));...
                SE_mean([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3));...
                SE_mean([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))];


            error_GP_allData(:,k) = abs(X_estimation-X_Fault_GP_allData(:,k)) ;

            if  max(error_GP_allData(:,k)) < 0.1
    %            disp("Validated")

                AB_k_allData = [SE_deriv([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_1, noisevar, sigma_f(1), sigma_l(1))';...
                    SE_deriv([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_2, noisevar, sigma_f(2), sigma_l(2))';...
                    SE_deriv([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_3, noisevar, sigma_f(3), sigma_l(3))';...
                    SE_deriv([X_Fault_GP_allData(:,k-1)', U_Fault_GP_allData(2,k-1)], Z_GP, Y_GP_4, noisevar, sigma_f(4), sigma_l(4))'
                    ];
                A_k_allData = AB_k_allData(:,1:4);
                B_k_allData = [zeros(4,1), AB_k_allData(:,5)];
                d_k_allData = X_Fault_GP_allData(:,k) - (A_k_allData*X_Fault_GP_allData(:,k-1) + B_k_allData*U_Fault_GP_allData(:,k-1));

                error_mat_allDataGP(k) = norm(A_k_allData-A);


                [U_Fault_GP_allData(:,k),infeasible] = ctrl_broken{X_Fault_GP_allData(:,k),[U_Fault_GP_allData(1,k_fault-1);U_Fault_GP_allData(2,k-1)], A_k_allData, B_k_allData, d_k_allData};
                X_Fault_GP_allData(:,k+1)=A*X_Fault_GP_allData(:,k)+B*[U_Fault_GP_allData(1,k_fault-1);U_Fault_GP_allData(2,k)];
            else
                disp("NOT Validated (more data)")
                %disp(max(abs(X_estimation-X_Fault_GP(:,k))./X_Fault_GP(:,k)))
                [U_Fault_GP_allData(:,k),infeasible] = ctrl_nominal{X_Fault_GP_allData(:,k),[U_Fault_GP_allData(1,k_fault-1);U_Fault_GP_allData(2,k-1)]};
                X_Fault_GP_allData(:,k+1)=A*X_Fault_GP_allData(:,k)+B*[U_Fault_GP_allData(1,k_fault-1);U_Fault_GP_allData(2,k)];
            end
        end     
    end
        
        %% ----------fault with L2NW-------------
        %Identifying and linearizing the system using L2NW, and a linear
        %kernel (not a very good idea)
    if do_L2NW_1 == true
        if k==1
            [U_Fault_L2NW(:,k),infeasible] = ctrl_nominal{X_Fault_L2NW(:,k),U0};
            X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*U_Fault_L2NW(:,k);
        elseif k < k_fault
            [U_Fault_L2NW(:,k),infeasible]=ctrl_nominal{X_Fault_L2NW(:,k),U_Fault_L2NW(:,k-1)};
            X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*U_Fault_L2NW(:,k);
        elseif k < k_switch
            [U_Fault_L2NW(:,k),infeasible] = ctrl_nominal{X_Fault_L2NW(:,k),[U_Fault_L2NW(1,k_fault-1);U_Fault_L2NW(2,k-1)]};
            X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*[U_Fault_L2NW(1,k_fault-1);U_Fault_L2NW(2,k)];
        else
            param = [3, 4; 3, 4; 3, 4; 3, 4]; %[h, lambda]

            Y_L2NW_1 = [X_Fault_L2NW(1,k_fault+1:k-1)';Y_1(end-data_number-1+k:end)];
            Y_L2NW_2 = [X_Fault_L2NW(2,k_fault+1:k-1)';Y_2(end-data_number-1+k:end)];
            Y_L2NW_3 = [X_Fault_L2NW(3,k_fault+1:k-1)';Y_3(end-data_number-1+k:end)];
            Y_L2NW_4 = [X_Fault_L2NW(4,k_fault+1:k-1)';Y_4(end-data_number-1+k:end)];
            Z_L2NW = [X_Fault_L2NW(:,k_fault:k-2)', U_Fault_L2NW(2,k_fault:k-2)'; Z( end-data_number-1+k:end,:) ];

            

            X_estimation_L2NW = [L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_1,param(1,:),1);...
                L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_2,param(2,:),1);...
                L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_3,param(3,:),1);...
                L2NW_mean([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_4,param(4,:),1)...
                ];


            error_L2NW(:,k) = abs([X_Fault_L2NW(:,k) -  X_estimation_L2NW]);

            if  max(error_L2NW(:,k)) < 300

                O_k = [L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_1,param(1,:),1)';...
                    L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_2,param(2,:),1)';...
                    L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_3,param(3,:),1)';...
                    L2NW_deriv([X_Fault_L2NW(:,k-1)', U_Fault_L2NW(2,k-1)],Z_L2NW,Y_L2NW_4,param(4,:),1)';...
                    ];
                O_A_k = O_k(:,1:4);
                O_B_k = [zeros(4,1), O_k(:,5)];
                O_d_k = X_Fault_L2NW(:,k) - (O_A_k*X_Fault_L2NW(:,k-1) + O_B_k*U_Fault_L2NW(:,k-1));

                error_mat_L2NW(k) = norm(A- O_A_k);

                [U_Fault_L2NW(:,k), infeasible] = ctrl_broken(X_Fault_L2NW(:,k),[U_Fault_L2NW(1,k_fault-1);U_Fault_L2NW(2,k-1)], O_A_k, O_B_k, O_d_k);
                if infeasible ~= 0
                    disp('problem')
                end

                X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*[U_Fault_L2NW(1,k_fault-1);U_Fault_L2NW(2,k)];
            else
                disp('L2NW  NOT validated')
                [U_Fault_L2NW(:,k),infeasible] = ctrl_nominal{X_Fault_L2NW(:,k),[U_Fault_L2NW(1,k_fault-1);U_Fault_L2NW(2,k-1)]};
                X_Fault_L2NW(:,k+1)=A*X_Fault_L2NW(:,k)+B*[U_Fault_L2NW(1,k_fault-1);U_Fault_L2NW(2,k)];


            end
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
            lambda = zeros(1,4);
            %lambda =    1.0e+05 *[0.000000451025780   7.979889364022138   0.024365540051243   0.086369157185887];
             %lambda =      1.0e+06 *[0.000000045105253   1.372971685457147   0.001704647349413   0.004564597558218];
    %         sigma_f = 2.17*10^(3)*ones(1,4);
    %         sigma_l = 1.14*10^(4)*ones(1,4);
    %         param = [10.695781461698118, 1.677582380126192, lambda; 99.935986868827754,  80.087113811947248, lambda;...
    %            9.046354763354064, 1.371605612665715,   lambda; 32.505687565655741 , 6.573873522600369, lambda];
    %         param = [sigma_l, sigma_f, lambda;...
    %             sigma_l, sigma_f, lambda;...
    %             sigma_l, sigma_f, lambda;...
    %             sigma_l, sigma_f, lambda;];
                %sigmaL sigmaF lambda
            sigma_f =   1.0e+02 *[ 2.639696625733967   0.170294093647814   0.085880808896654   0.101553751261396];
            sigma_l =[ 99.991949107543803  25.640539532091616  34.284915666906144  25.862479178807501];
%             sigma_f =[  20.575591520093276,  17.878167302000445,  16.199466399148097,  17.938357317647856];
%             sigma_l =[  99.966954902978543,  34.241347755939984,  78.546996930435938,  88.147249188528974];
    %         sigma_f = [1.677582380126192, 80.087113811947248, 1.371605612665715, 6.573873522600369];
    %         sigma_l = [10.695781461698118, 99.935986868827754,  9.046354763354064, 32.505687565655741];

            param = [sigma_l(1), sigma_f(1), lambda(1);...
                sigma_l(2), sigma_f(2), lambda(2);...
                sigma_l(3), sigma_f(3), lambda(3);...
                sigma_l(4), sigma_f(4), lambda(4)];


            data_number_l2 = 100;
            Y_L2NW_1GP = [X_Fault_L2NWGP(1,k_fault+1:k-1)'; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_2GP = [X_Fault_L2NWGP(2,k_fault+1:k-1)'; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_3GP = [X_Fault_L2NWGP(3,k_fault+1:k-1)'; Y_1(end -data_number_l2-1+k:end)];
            Y_L2NW_4GP = [X_Fault_L2NWGP(4,k_fault+1:k-1)'; Y_1(end -data_number_l2-1+k:end)];
            Z_L2NW_GP = [X_Fault_L2NWGP(:,k_fault:k-2)', U_Fault_L2NWGP(2,k_fault:k-2)';  Z( end -data_number_l2-1+k:end,:)];
            


            X_estimation_L2NWGP = [L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_1GP,param(1,:),2);...
                L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_2GP,param(2,:),2);...
                L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_3GP,param(3,:),2);...
                L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_4GP,param(4,:),2)...
                ];


            error_L2NWGP(:,k) = abs(X_Fault_L2NWGP(:,k) -  X_estimation_L2NWGP);

            if  max(error_L2NWGP(:,k)) < 300
                %Casadi problem: 
                opti = casadi.Opti();
                opti.subject_to; %clean the constraint?
                opti.variable;
                
                tic
                x_tild = opti.variable(4, horizon);
                u_tild = opti.variable(1, horizon);

                u_old = opti.parameter(1,1);
                x_0 = opti.parameter(4,1);


                opti.subject_to( x_tild(:,1) == x_0);
                opti.subject_to( C_u_delta_redu*[u_old;u_tild(1)] <= c_u_delta_redu )
                for i = 1:horizon-1
                    opti.subject_to( C_u_redu*u_tild(:,i) <= c_u_redu );
                    opti.subject_to( C_u_delta_redu*[u_tild(i);u_tild(i+1)] <= c_u_delta_redu );
                    opti.subject_to( x_tild(:,i+1) ==  [L2NW_mean([x_tild(:,i)',u_tild(i)] , Z_L2NW_GP, Y_L2NW_1GP, param(1,:), 3);...
                        L2NW_mean([x_tild(:,i)',u_tild(i)], Z_L2NW_GP, Y_L2NW_2GP, param(2,:), 3);...
                        L2NW_mean([x_tild(:,i)',u_tild(i)], Z_L2NW_GP, Y_L2NW_3GP, param(3,:), 3);...
                        L2NW_mean([x_tild(:,i)',u_tild(i)], Z_L2NW_GP, Y_L2NW_4GP, param(4,:), 3)]);

                end
                opti.minimize(x_tild(:)'*Q_huge*x_tild(:) + u_tild(:)'*R_huge*u_tild(:));

                ops = struct;
                ops.ipopt.print_level = 0;
                ops.ipopt.tol = 1e-3;
                opti.solver('ipopt', ops);


                opti.set_value(x_0,X_Fault_L2NWGP(:,k));
                opti.set_value(u_old,U_Fault_L2NWGP(2,k));
                sol = opti.solve();
                U_Fault_L2NWGP(2,k) = sol.value(u_tild(1));
                U_Fault_L2NWGP(1,k) = U_Fault_L2NWGP(1,k-1);        %there is a fault
                toc
                %end cadasi
                
                

                X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
            else
                disp('L2NW2  NOT validated')
                [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k-1)]};
                X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];


            end
        end
    end
      

        
    %% ---------------ALL----------------------------------
    if mod(k,10) == 0
        disp(k)
    end
    
end





%% Plotting ================================================================
t=Ts*(0:stepNumber-1);
%% States
figure
hold on
leg = [];

subplot(2,2,1) %State 1
hold on
title('Pitch Angle [rad]')
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
if do_Fault_GP_allData == true
    plot(t,X_Fault_GP_moreData(1,:),'b:')
end
if do_L2NW_1 == true
    plot(t,X_Fault_L2NW(1,:),'r-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(1,:),'r-')
end

subplot(2,2,2) %State 2-----------------
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
if do_Fault_GP_allData == true
    plot(t,X_Fault_GP_moreData(2,:),'b:')
end
if do_L2NW_1 == true
    plot(t,X_Fault_L2NW(2,:),'r-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(2,:),'r-')
end

subplot(2,2,3) %State 3---------------
hold on
title('Pitch rate [rad/s]')
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
if do_Fault_GP_allData == true
    plot(t,X_Fault_GP_moreData(3,:),'b:')
end
if do_L2NW_1 == true
    plot(t,X_Fault_L2NW(3,:),'r-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(3,:),'r-')
end
    
subplot(2,2,4) %State 4---------
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
if do_Fault_GP_allData == true
    plot(t,X_Fault_GP_moreData(4,:),'b:')
end
if do_L2NW_1 == true
    plot(t,X_Fault_L2NW(4,:),'r-.')
end
if do_L2NW_GP == true
    plot(t,X_Fault_L2NWGP(4,:),'r-')
end

%legend:
if do_noFault == true  
    leg = [leg, "no Fault scenaria"];
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
if do_Fault_GP_allData == true
    leg = [leg, "Fault with GP model (with all data)"];
end
if do_L2NW_1 == true
    leg = [leg, "Fault with L2NW model"];
end
if do_L2NW_GP == true
    leg = [leg, "Fault with gaussian L2NW model"];
end
legend(leg)
%% error --------------------------------------------------------------------------
figure
hold on
leg = [];

subplot(2,2,1) %State 1
hold on
title('Error on pitch Angle [rad]')
xlabel('time s')
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(1,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(1,:),'b-.')
end
if do_Fault_GP_allData == true
    plot(t(1:end-1),error_GP_moreData(1,:),'b:')
end
if do_L2NW_1 == true
    plot(t(1:end-1),error_L2NW(1,:),'r-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(1,:),'r-')
end

subplot(2,2,2) %State 2-----------------
hold on
title('Error on air speed [m/s]')
xlabel('time s')
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(2,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(2,:),'b-.')
end
if do_Fault_GP_allData == true
    plot(t(1:end-1),error_GP_moreData(2,:),'b:')
end
if do_L2NW_1 == true
    plot(t(1:end-1),error_L2NW(2,:),'r-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(2,:),'r-')
end

subplot(2,2,3) %State 3---------------
hold on
title('Error on pitch rate [rad/s]')
xlabel('time s')
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(3,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(3,:),'b-.')
end
if do_Fault_GP_allData == true
    plot(t(1:end-1),error_GP_moreData(3,:),'b:')
end
if do_L2NW_1 == true
    plot(t(1:end-1),error_L2NW(3,:),'r-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(3,:),'r-')
end
    
subplot(2,2,4) %State 4---------
hold on
title('Error on angle of attack [rad]')
xlabel('time s')
if do_Fault_GP == true
    plot(t(1:end-1),error_GP(4,:),'b-')    
end
if do_Fault_GP_moreData == true
    plot(t(1:end-1),error_GP_moreData(4,:),'b-.')
end
if do_Fault_GP_allData == true
    plot(t(1:end-1),error_GP_moreData(4,:),'b:')
end
if do_L2NW_1 == true
    plot(t(1:end-1),error_L2NW(4,:),'r-.')
end
if do_L2NW_GP == true
    plot(t(1:end-1),error_L2NWGP(4,:),'r-')
end

%legend:
leg = [];
if do_Fault_GP == true
    leg = [leg, "Fault with GP model"];    
end
if do_Fault_GP_moreData == true
    leg = [leg, "Fault with GP model (with more data)"];
end
if do_Fault_GP_allData == true
    leg = [leg, "Fault with GP model (with all data)"];
end
if do_L2NW_1 == true
    leg = [leg, "Fault with L2NW model"];
end
if do_L2NW_GP == true
    leg = [leg, "Fault with gaussian L2NW model"];
end
legend(leg)


%%
% figure
% subplot(2,1,1)
% plot(t(1:end-1),U_noFault(1,:),t(1:end-1),U_Fault(1,:),...
%     t(1:end-1),U_Fault_Perfect(1,:),'o',t(1:end-1),U_Fault_GP(1,:),...
%     t(1:end-1),U_Fault_GP_moreData(1,:),t(1:end-1),U_Fault_L2NWGP(1,:))
% title('Elevator [deg]')
% xlabel('time s')
% subplot(2,1,2)
% plot(t(1:end-1),U_noFault(2,:),t(1:end-1),U_Fault(2,:)...
%     ,t(1:end-1),U_Fault_Perfect(2,:),'o',t(1:end-1),U_Fault_GP(2,:),...
%     t(1:end-1),U_Fault_GP_moreData(2,:),t(1:end-1),U_Fault_L2NWGP(2,:))
% title('THS [deg] ')
% xlabel('time s')
% legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP',  'GP with past data')
% 
% % figure
% % subplot(2,1,1)
% % plot(t(1:end-1),matrix_ratio(1,:))
% % title('A_k')
% % subplot(2,1,2)
% % plot(t(1:end-1),matrix_ratio(2,:))
% % title('B_k')
% 
% 
% dU_noFault = (U_noFault(:,2:end)-U_noFault(:,1:end-1))/Ts;
% dU_Fault = (U_Fault(:,2:end)-U_Fault(:,1:end-1))/Ts;
% dU_Fault_Perfect = (U_Fault_Perfect(:,2:end)-U_Fault_Perfect(:,1:end-1))/Ts;
% dU_Fault_GP = (U_Fault_GP(:,2:end)-U_Fault_GP(:,1:end-1))/Ts;
% dU_Fault_L2NW = (U_Fault_L2NW(:,2:end)-U_Fault_L2NW(:,1:end-1))/Ts;
% dU_Fault_L2NWGP = (U_Fault_L2NWGP(:,2:end)-U_Fault_L2NWGP(:,1:end-1))/Ts;
% dU_Fault_GP_moreData = (U_Fault_GP_moreData(:,2:end)-U_Fault_GP_moreData(:,1:end-1))/Ts;
% dU_Fault_GP_allData = (U_Fault_GP_allData(:,2:end)-U_Fault_GP_allData(:,1:end-1))/Ts;
% 
% 
% figure
% subplot(2,1,1)
% plot(t(1:end-2),dU_noFault(1,:),t(1:end-2),dU_Fault(1,:),...
%     t(1:end-2),dU_Fault_Perfect(1,:),'o',t(1:end-2),dU_Fault_L2NWGP(1,:),...
%     t(1:end-2),dU_Fault_GP_moreData(1,:))
% title('Elevator slew rate [deg/s]')
% xlabel('time s')
% subplot(2,1,2)
% plot(t(1:end-2),dU_noFault(2,:),t(1:end-2),dU_Fault(2,:),...
%     t(1:end-2),dU_Fault_GP(2,:),'o',t(1:end-2),dU_Fault_L2NWGP(2,:),...
%     t(1:end-2),dU_Fault_GP_moreData(2,:))
% title('THS slew rate [deg/s]')
% xlabel('time s')
% legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP', 'L2NWGP')




figure
plot(t(1:end-1),error_mat_GP(1,:),t(1:end-1),error_mat_GP_moreData(1,:),t(1:end-1))
title('norm 2 error on matrix')
xlabel('time [s]')
legend('GP',  'GP with past data')

%----------------SAVE-----------
% figure
% subplot(2,2,1)
% plot(t,X_noFault(1,:),t,X_Fault(1,:),t,X_Fault_Perfect(1,:),'o',t,X_Fault_GP(1,:),t,X_Fault_L2NW(1,:)...
%     ,t,X_Fault_L2NWGP(1,:))
% title('Pitch Angle rad/s')
% xlabel('time s')
% subplot(2,2,2)
% plot(t,X_noFault(2,:),t,X_Fault(2,:),t,X_Fault_Perfect(2,:),'o',t,X_Fault_GP(2,:),t,X_Fault_L2NW(2,:)...
%     ,t,X_Fault_L2NWGP(2,:))
% title('Air speed m/s')
% xlabel('time s')
% subplot(2,2,3)
% plot(t,X_noFault(3,:),t,X_Fault(3,:),t,X_Fault_Perfect(3,:),'o',t,X_Fault_GP(3,:),t,X_Fault_L2NW(3,:)...
%     ,t,X_Fault_L2NWGP(3,:))
% title('Pitch Rate rad')
% xlabel('time s')
% subplot(2,2,4)
% plot(t,X_noFault(4,:),t,X_Fault(4,:),t,X_Fault_Perfect(4,:),'o',t,X_Fault_GP(4,:),t,X_Fault_L2NW(4,:)...
%     ,t,X_Fault_L2NWGP(4,:))
% title('AoA rad')
% xlabel('time s')
% legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP', 'L2NW_1', 'L2NW_2')
% 
% figure
% subplot(2,1,1)
% plot(t(1:end-1),U_noFault(1,:),t(1:end-1),U_Fault(1,:),...
%     t(1:end-1),U_Fault_Perfect(1,:),'o',t(1:end-1),U_Fault_GP(1,:),...
%     t(1:end-1),U_Fault_L2NW(1,:),t(1:end-1),U_Fault_L2NWGP(1,:))
% title('Elevator deg')
% xlabel('time s')
% subplot(2,1,2)
% plot(t(1:end-1),U_noFault(2,:),t(1:end-1),U_Fault(2,:)...
%     ,t(1:end-1),U_Fault_Perfect(2,:),'o',t(1:end-1),U_Fault_GP(2,:),...
%     t(1:end-1),U_Fault_L2NW(2,:),t(1:end-1),U_Fault_L2NWGP(2,:))
% title('THS deg ')
% xlabel('time s')
% legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP', 'L2NW_1', 'L2NW_2')
% 
% % figure
% % subplot(2,1,1)
% % plot(t(1:end-1),matrix_ratio(1,:))
% % title('A_k')
% % subplot(2,1,2)
% % plot(t(1:end-1),matrix_ratio(2,:))
% % title('B_k')
% 
% 
% dU_noFault = (U_noFault(:,2:end)-U_noFault(:,1:end-1))/Ts;
% dU_Fault = (U_Fault(:,2:end)-U_Fault(:,1:end-1))/Ts;
% dU_Fault_Perfect = (U_Fault_Perfect(:,2:end)-U_Fault_Perfect(:,1:end-1))/Ts;
% dU_Fault_GP = (U_Fault_GP(:,2:end)-U_Fault_GP(:,1:end-1))/Ts;
% dU_Fault_L2NW = (U_Fault_L2NW(:,2:end)-U_Fault_L2NW(:,1:end-1))/Ts;
% dU_Fault_L2NWGP = (U_Fault_L2NWGP(:,2:end)-U_Fault_L2NWGP(:,1:end-1))/Ts;
% 
% figure
% subplot(2,1,1)
% plot(t(1:end-2),dU_noFault(1,:),t(1:end-2),dU_Fault(1,:),...
%     t(1:end-2),dU_Fault_Perfect(1,:),'o',t(1:end-2),dU_Fault_GP(1,:),...
%     t(1:end-2),dU_Fault_L2NW(1,:),t(1:end-2),dU_Fault_L2NWGP(1,:))
% title('Elevator slew rate def/s')
% xlabel('time s')
% subplot(2,1,2)
% plot(t(1:end-2),dU_noFault(2,:),t(1:end-2),dU_Fault(2,:),...
%     t(1:end-2),dU_Fault_GP(2,:),'o',t(1:end-2),dU_Fault_GP(2,:),...
%     t(1:end-2),dU_Fault_L2NW(2,:),t(1:end-2),dU_Fault_L2NWGP(2,:))
% title('THS slew rate deg/s')
% xlabel('time s')
% legend('No Fault', 'Fault with nominal MPC', 'Fault with Perfect MPC', 'Fault with GP', 'L2NW_1', 'L2NW_2')
% 
% figure
% 
% subplot(2,2,1)
% plot(t(1:end-1),error_GP(1,:),t(1:end-1),error_L2NW(1,:),t(1:end-1),error_L2NWGP(1,:))
% title('Error on Pitch Angle rad/s')
% xlabel('time s')
% %legend('GP', 'L2NW 1')
% subplot(2,2,2)
% plot(t(1:end-1),error_GP(2,:),t(1:end-1),error_L2NW(2,:),t(1:end-1),error_L2NWGP(2,:))
% title('Error on Air speed m/s')
% xlabel('time s')
% %legend('GP', 'L2NW 1')
% subplot(2,2,3)
% plot(t(1:end-1),error_GP(3,:),t(1:end-1),error_L2NW(3,:),t(1:end-1),error_L2NWGP(3,:))
% title('Error on Pitch Rate rad')
% xlabel('time s')
% %legend('GP', 'L2NW 1')
% subplot(2,2,4)
% plot(t(1:end-1),error_GP(4,:),t(1:end-1),error_L2NW(4,:),t(1:end-1),error_L2NWGP(4,:))
% title('Error on AoA rad')
% xlabel('time s')
% legend('GP', 'L2NW_1', 'L2NW_2')
% 
% %%
% figure
% plot(t(1:end-1),error_mat_GP(1,:),t(1:end-1),error_mat_L2NW(1,:),t(1:end-1),error_mat_L2NWGP(1,:))
% title('norm 2 error on matrix')
% xlabel('time s')
% legend('GP', 'L2NW_1', 'L2NW_2')































%=================SAVE=====================================================
%             %% ----------fault with L2NW and GP Linearized-------------
%     if do_L2NW_1 == true
%         if k==1
%             [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),U0};
%             X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*U_Fault_L2NWGP(:,k);
%         elseif k < k_fault
%             [U_Fault_L2NWGP(:,k),infeasible]=ctrl_nominal{X_Fault_L2NWGP(:,k),U_Fault_L2NWGP(:,k-1)};
%             X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*U_Fault_L2NWGP(:,k);
%         elseif k < k_switch
%             [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k-1)]};
%             X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
%         else
% 
%              lambda =    [0.000041399654408   0.004332991324340   0.000013658470911   0.000019700973465];
%     %         sigma_f = 2.17*10^(3)*ones(1,4);
%     %         sigma_l = 1.14*10^(4)*ones(1,4);
%     %         param = [10.695781461698118, 1.677582380126192, lambda; 99.935986868827754,  80.087113811947248, lambda;...
%     %            9.046354763354064, 1.371605612665715,   lambda; 32.505687565655741 , 6.573873522600369, lambda];
%     %         param = [sigma_l, sigma_f, lambda;...
%     %             sigma_l, sigma_f, lambda;...
%     %             sigma_l, sigma_f, lambda;...
%     %             sigma_l, sigma_f, lambda;];
%                 %sigmaL sigmaF lambda
% 
%             sigma_f =[  20.575591520093276,  17.878167302000445,  16.199466399148097,  17.938357317647856];
%             sigma_l =[  99.966954902978543,  34.241347755939984,  78.546996930435938,  88.147249188528974];
%     %         sigma_f = [1.677582380126192, 80.087113811947248, 1.371605612665715, 6.573873522600369];
%     %         sigma_l = [10.695781461698118, 99.935986868827754,  9.046354763354064, 32.505687565655741];
% 
%             param = [sigma_l(1), sigma_f(1), lambda(1);...
%                 sigma_l(2), sigma_f(2), lambda(2);...
%                 sigma_l(3), sigma_f(3), lambda(3);...
%                 sigma_l(4), sigma_f(4), lambda(4)];
% 
% 
% 
%             Y_L2NW_1GP = X_Fault_L2NWGP(1,k_fault+1:k-1)';
%             Y_L2NW_2GP = X_Fault_L2NWGP(2,k_fault+1:k-1)';
%             Y_L2NW_3GP = X_Fault_L2NWGP(3,k_fault+1:k-1)';
%             Y_L2NW_4GP = X_Fault_L2NWGP(4,k_fault+1:k-1)';
%             Z_L2NW_GP = [X_Fault_L2NWGP(:,k_fault:k-2)', U_Fault_L2NWGP(2,k_fault:k-2)'];
% 
%             X_estimation_L2NWGP = [L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_1GP,param(1,:),2);...
%                 L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_2GP,param(2,:),2);...
%                 L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_3GP,param(3,:),2);...
%                 L2NW_mean([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_4GP,param(4,:),2)...
%                 ];
% 
% 
%             error_L2NWGP(:,k) = abs(X_Fault_L2NWGP(:,k) -  X_estimation_L2NWGP);
% 
%             if  max(error_L2NWGP(:,k)) < 300
% 
%                 O_kGP = [L2NW_deriv([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_1GP,param(1,:),2)';...
%                     L2NW_deriv([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_2GP,param(2,:),2)';...
%                     L2NW_deriv([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_3GP,param(3,:),2)';...
%                     L2NW_deriv([X_Fault_L2NWGP(:,k-1)', U_Fault_L2NWGP(2,k-1)],Z_L2NW_GP,Y_L2NW_4GP,param(4,:),2)';...
%                     ];
%                 O_A_kGP = O_kGP(:,1:4);
%                 O_B_kGP = [zeros(4,1), O_kGP(:,5)];
%                 O_d_kGP = X_Fault_L2NWGP(:,k) - (O_A_kGP*X_Fault_L2NWGP(:,k-1) + O_B_kGP*U_Fault_L2NWGP(:,k-1));
% 
%                 error_mat_L2NWGP(k) = norm(A-O_A_kGP);
% 
%                 [U_Fault_L2NWGP(:,k), infeasible] = ctrl_broken(X_Fault_L2NWGP(:,k),[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k-1)], O_A_kGP, O_B_kGP, O_d_kGP);
%                 if infeasible ~= 0
%                     disp('problem2')
%                 end
% 
%                 X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
%             else
%                 disp('L2NW2  NOT validated')
%                 [U_Fault_L2NWGP(:,k),infeasible] = ctrl_nominal{X_Fault_L2NWGP(:,k),[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k-1)]};
%                 X_Fault_L2NWGP(:,k+1)=A*X_Fault_L2NWGP(:,k)+B*[U_Fault_L2NWGP(1,k_fault-1);U_Fault_L2NWGP(2,k)];
% 
% 
%             end
%         end
%     end
%       
