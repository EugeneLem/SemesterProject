format long

%%Constant Setting
%System
A=[0.86719   ,6.6936*10^(-5)     ,-0.19095  ,0; ...
    -0.027773, 0.99895           , 0.89264, -1.9609; ...
      0.20146, -2.1676*10^(-4)     , 0.88379,  0; ....
      0.2    ,  0                  ,  0      , 1];
B=[-3.7758 *10^(-3)  , -9.0408*10^(-3);...
    0                ,0;...
      -1.2629*10^(-3), -3.2794*10^(-4);...
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
C_u=[1 0 ; 0,1 ; -1,0 ; 0,-1];
c_u=[15;4.6;25;10.4];%Need to check is elevator and THS are not inverted
C_u_delta=[-1 0   1 0;...
            0 -1  0 1;...
            1  0 -1 0;...
            0  1  0 -1];
c_u_delta=[37;0.236;37;0.236]*Ts;%Need to check is elevator and THS are not inverted

%%MPC
x=sdpvar(4,horizon);
u=sdpvar(2,horizon);
A_k=sdpvar(4,4);
B_k=sdpvar(4,2);
u_old=sdpvar(2,1);%For the input change constraint, we need the 0 input
ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
con=[C_u_delta*[u_old;u(:,1)] <= c_u_delta];
obj = 0;   
for i = 1:horizon-1 
    con = [con , x(:,i+1) == A_k*x(:,i) + B_k*u(:,i)]; % System dynamics
    con = [con , C_u*u(:,i) <= c_u];                 % Input constraints
    con = [con , C_u_delta*[u(:,i);u(:,i+1)] <= c_u_delta];  % Input change constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i) + (u(:,i+1)-u(:,i))'*R_delta*(u(:,i+1)-u(:,i));  % Cost function

end
obj=obj+[u(:,horizon);x(:,horizon)]'*P_k*[u(:,horizon);x(:,horizon)];
ctrl = optimizer(con, obj,ops, {x(:,1),A_k,B_k,u_old}, u(:,1));      



%%Let's simulate our system
stepNumber=300;
k_fault=40;
k_switch=k_fault+1;  %The time we decide to use GP
X=zeros(4,stepNumber); 
U=zeros(2,stepNumber-1); 
X(:,1)=[1;1;1;1];
U_0=[0;0];

%GP numbers
inver_wid = 1/88;     %Magic number from Harsh
noiseVar = 0;         %No noise here


for k=1:stepNumber-1                    
    if(k<=k_switch) %We use a classic MPC =)
        if k==1
            [U(:,k),infeasible]=ctrl{X(:,k),A,B,U_0};
        else
            [U(:,k),infeasible]=ctrl{X(:,k),A,B,U(:,k-1)};
        end
    
    
    
    else
        % Kernel creation
        yTrain_1 =  X(1,k_fault+1:k-1)';
        yTrain_2 =  X(2,k_fault+1:k-1)';
        yTrain_3 =  X(3,k_fault+1:k-1)';
        yTrain_4 =  X(4,k_fault+1:k-1)';
        zTrain =   [X(:,k_fault+1:k-1);U(:,k_fault+1:k-1)]';  
            
        trueKer_1 = kernCreate(zTrain, 'rbf');
        trueKer_1.inverseWidth = 1/inver_wid;
        %trueKer_1.variance = 39.49;
        
        trueKer_2 = kernCreate(zTrain, 'rbf');
        trueKer_2.inverseWidth = 1/inver_wid;
        %trueKer_2.variance = 39.49;
        
        trueKer_3 = kernCreate(zTrain, 'rbf');
        trueKer_3.inverseWidth = 1/inver_wid;
        %trueKer_3.variance = 39.49;   
        
        trueKer_4 = kernCreate(zTrain, 'rbf');
        trueKer_4.inverseWidth = 1/inver_wid;
        %trueKer_4.variance = 39.49;        
  
        
        %Determine A_k and B_k:
        
        A_k = [GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_1, trueKer_1, noiseVar);... 
            GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_2, trueKer_2, noiseVar);... 
            GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_3, trueKer_3, noiseVar);... 
            GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_4, trueKer_4, noiseVar)];
        
        B_k=[GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_1, trueKer_1, noiseVar);...
            GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_2, trueKer_2, noiseVar);... 
            GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_3, trueKer_3, noiseVar);... 
            GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_4, trueKer_4, noiseVar)];
        
        %Model validation
        error_nominal=mean(vecnorm( X(:,k_fault:k)-A*X(:,k_fault-1:k-1)+B*U(:,k_fault-1:k-1) ));
        error_GP=mean(vecnorm( X(:,k_fault:k)-A_k*X(:,k_fault-1:k-1)+B_k*U(:,k_fault-1:k-1) ));
        
%         error_GP=zeros(1,k-k_fault);
%         for i=k_fault:k
%             estimation=[GP_SE_mean_var( [X(:,i-1)',U(:,i-1)'], zTrain, yTrain_1, trueKer_1, noiseVar);
%                 GP_SE_mean_var( [X(:,i-1)',U(:,i-1)'], zTrain, yTrain_2, trueKer_2, noiseVar);...
%                 GP_SE_mean_var( [X(:,i-1)',U(:,i-1)'], zTrain, yTrain_3, trueKer_3, noiseVar);...
%                 GP_SE_mean_var( [X(:,i-1)',U(:,i-1)'], zTrain, yTrain_4, trueKer_4, noiseVar)];
%             
%             
%             error_GP(i+1-k_fault)=norm( X(:,i) - estimation );
%         end
%         error_GP=mean(error_GP);
        
        
        if error_GP < 10* error_nominal
            disp('Model Validated')
            k
            [U(:,k),infeasible]=ctrl{X(:,k),A_k,B_k,U(:,k-1)};
        else
            [U(:,k),infeasible]=ctrl{X(:,k),A,B,U(:,k-1)};
        end 
    
    end
            
    %disp(k)
    %disp(infeasible)
    if(k>=k_fault)
        U(1,k)=U(1,k-1);
    end
    X(:,k+1)=A*X(:,k)+B*U(:,k);
end


%%
%%PLOT :D
t=Ts*(0:stepNumber-1);
figure
subplot(3,1,1)
plot(t,X(1,:))
title('Pitch Angle')
subplot(3,1,2)
plot(t,X(3,:))
title('Pitch Rate')
subplot(3,1,3)
plot(t,X(4,:))
title('AoA')

figure
subplot(2,1,1)
plot(t(1:end-1),U(1,:))
title('Elevator')
subplot(2,1,2)
plot(t(1:end-1),U(2,:))
title('THS')

% figure
% subplot(2,1,1)
% plot(t(1:end-1),matrix_ratio(1,:))
% title('A_k')
% subplot(2,1,2)
% plot(t(1:end-1),matrix_ratio(2,:))
% title('B_k')


dU=(U(:,2:end)-U(:,1:end-1))/Ts;
figure
subplot(2,1,1)
plot(t(1:end-2),dU(1,:))
title('Elevator slew rate')
subplot(2,1,2)
plot(t(1:end-2),dU(2,:))
title('THS slew rate')
