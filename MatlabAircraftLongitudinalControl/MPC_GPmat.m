format long

%%Constant Setting
%System
A=[0.86719   ,6.6936*10^(-5)     ,-0.19095  ,0; ...
    -0.027773, 0.99895           , 0.89264, -1.9609; ...
      0.20146, -2.1676*10^(-4)     , 0.88379,  0; ...
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
C_u = [1,0 ; 0,1 ; -1,0 ; 0,-1];
c_u = [15;4.6;25;10.4];%Need to check is elevator and THS are not inverted
C_u_delta = [-1 0   1 0;...
            0 -1  0 1;...
            1  0 -1 0;...
            0  1  0 -1];
c_u_delta = [37;0.236;37;0.236]*Ts;%Need to check is elevator and THS are not inverted

%%MPC
x=sdpvar(4,horizon);
u=sdpvar(2,horizon);
A_k=sdpvar(4,4);
B_k=sdpvar(4,2);
d_k=sdpvar(4,1);
u_old=sdpvar(2,1);%For the input change constraint, we need the 0 input
ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
con = C_u_delta*[u_old;u(:,1)] <= c_u_delta;
obj = 0;   
for i = 1:horizon-1 
    con = [con , x(:,i+1) == A_k*x(:,i) + B_k*u(:,i)+d_k]; % System dynamics
    con = [con , C_u*u(:,i) <= c_u];                 % Input constraints
    con = [con , C_u_delta*[u(:,i);u(:,i+1)] <= c_u_delta];  % Input change constraints
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i) + (u(:,i+1)-u(:,i))'*R_delta*(u(:,i+1)-u(:,i));  % Cost function

end
obj = obj+[u(:,horizon);x(:,horizon)]'*P_k*[u(:,horizon);x(:,horizon)];
ctrl = optimizer(con, obj,ops, {x(:,1),A_k,B_k,d_k,u_old}, u(:,1));      

 




%%Let's simulate our system
stepNumber=700;
k_fault=40;
k_switch=k_fault+1;  %The time we decide to use GP
X=zeros(4,stepNumber); 
U=zeros(2,stepNumber-1); 
X(:,1)=[1;0.1;1;1];
U_0=[0;0];

error_nominal = zeros(4,stepNumber);
error_GPmat = zeros(4,stepNumber);
error_GP = zeros(4,stepNumber);
error_GPml = zeros(4,stepNumber);
error_homemade = zeros(4,stepNumber);


%GP numbers
inver_wid = 1.157508543322349*10^4;%1.136884422110553;   %1/88;     %Magic number from Harsh 88
noiseVar = 0.01;         
model_var = (4.599670410677249*10^4);%1.157508543322349*10^4;%0.001569881746567^2;

sigma_l = 1/4.599670410677249*10^4;
sigma_f = 1.157508543322349*10^4;




for k=1:stepNumber-1                    
    if(k<=k_switch) %We use a classic MPC =)
        if k==1
            [U(:,k),infeasible]=ctrl{X(:,k),A,B,zeros(4,1),U_0};
        else
            [U(:,k),infeasible]=ctrl{X(:,k),A,B,zeros(4,1),U(:,k-1)};
        end
    
    
    
    else
        % Kernel creation
        
        %========GPmat=====================================================
        yTrain_1 =  X(1,k_fault+1:k-1)';
        yTrain_2 =  X(2,k_fault+1:k-1)';
        yTrain_3 =  X(3,k_fault+1:k-1)';
        yTrain_4 =  X(4,k_fault+1:k-1)';
        zTrain =   [X(:,k_fault:k-2);U(:,k_fault:k-2)]'; 


%         trueKer_1 = kernCreate(zTrain, 'rbf');
%         trueKer_1.inverseWidth = inver_wid *ones(1,1);
%         trueKer_1.variance = model_var;        %Value found with GP from matlab toolbox
% 
%         trueKer_2 = kernCreate(zTrain, 'rbf');
%         trueKer_2.inverseWidth = inver_wid *ones(1,1);
%         trueKer_2.variance = model_var;
% 
%         trueKer_3 = kernCreate(zTrain, 'rbf');
%         trueKer_3.inverseWidth = inver_wid *ones(1,1);
%         trueKer_3.variance = model_var; 
% 
%         trueKer_4 = kernCreate(zTrain, 'rbf');
%         trueKer_4.inverseWidth = inver_wid *ones(1,1);
%         trueKer_4.variance = model_var;     
% 
%         error_GPmat(1,k) = abs( GP_SE_mean_var([X(:,k-1);U(:,k-1)]', zTrain, yTrain_1, trueKer_1, noiseVar) - X(1,k));
%         error_GPmat(2,k) = abs( GP_SE_mean_var([X(:,k-1);U(:,k-1)]', zTrain, yTrain_2, trueKer_2, noiseVar) - X(2,k));
%         error_GPmat(3,k) = abs( GP_SE_mean_var([X(:,k-1);U(:,k-1)]', zTrain, yTrain_3, trueKer_3, noiseVar) - X(3,k));
%         error_GPmat(4,k) = abs( GP_SE_mean_var([X(:,k-1);U(:,k-1)]', zTrain, yTrain_4, trueKer_4, noiseVar) - X(4,k));
%       

        %========Base toolbase=============================================
%         gprMdl_1 = fitrgp([X(:,k_fault:k-2);U(:,k_fault:k-2)]',X(1,k_fault+1:k-1)','KernelFunction','squaredexponential'); %creation of the GP
%         gprMdl_2 = fitrgp([X(:,k_fault:k-2);U(:,k_fault:k-2)]',X(2,k_fault+1:k-1)','KernelFunction','squaredexponential'); %creation of the GP
%         gprMdl_3 = fitrgp([X(:,k_fault:k-2);U(:,k_fault:k-2)]',X(3,k_fault+1:k-1)','KernelFunction','squaredexponential');%creation of the GP
%         gprMdl_4 = fitrgp([X(:,k_fault:k-2);U(:,k_fault:k-2)]',X(4,k_fault+1:k-1)','KernelFunction','squaredexponential'); %creation of the GP
%         
%         
%         error_GP(1,k) = abs( predict(gprMdl_1, [X(:,k-1);U(:,k-1)]' ) - X(1,k) );
%         error_GP(2,k) = abs( predict(gprMdl_2, [X(:,k-1);U(:,k-1)]' ) - X(2,k) );
%         error_GP(3,k) = abs( predict(gprMdl_3, [X(:,k-1);U(:,k-1)]' ) - X(3,k) );
%         error_GP(4,k) = abs( predict(gprMdl_4, [X(:,k-1);U(:,k-1)]' ) - X(4,k) );
        
        %========GPml=====================================================
        
%         meanfunc = [];                    % empty: don't use a mean function
%         covfunc = @covSEiso;              % Squared Exponental covariance function
%         likfunc = @likGauss;              % Gaussian likelihood
%         
%         hyp_0 = struct('mean', [], 'cov', [9.314588383482750 6.890117306424340], 'lik', -9.267609631095468);
% %         hyp_1 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_1);
% %         hyp_2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_2);
% %         hyp_3 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_3);
% %         hyp_4 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_4);
% 
%         
%         error_GPml(1,k) = abs( gp(hyp_0, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_1, [X(:,k-1);U(:,k-1)]') - X(1,k) );
%         error_GPml(2,k) = abs( gp(hyp_0, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_2, [X(:,k-1);U(:,k-1)]')  - X(2,k) );
%         error_GPml(3,k) = abs( gp(hyp_0, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_3, [X(:,k-1);U(:,k-1)]')  - X(3,k) );
%         error_GPml(4,k) = abs( gp(hyp_0, @infGaussLik, meanfunc, covfunc, likfunc, zTrain, yTrain_4, [X(:,k-1);U(:,k-1)]')  - X(4,k) );
%         
        %========Homemade==================================================

        
        error_homemade(1,k) = abs( meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_1, noiseVar, sigma_f, sigma_l) - X(1,k) );
        error_homemade(2,k) = abs( meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_2, noiseVar, sigma_f, sigma_l) - X(2,k) );
        error_homemade(3,k) = abs( meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_3, noiseVar, sigma_f, sigma_l) - X(3,k) );
        error_homemade(4,k) = abs( meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_4, noiseVar, sigma_f, sigma_l) - X(4,k) );
  
        
        %Determine A_k and B_k: GPmat======================================


%         M = [GP_SE_Derivative([X(:,k-1);U(:,k-1)]', zTrain, yTrain_1, trueKer_1, noiseVar)';...
%             GP_SE_Derivative([X(:,k-1);U(:,k-1)]', zTrain, yTrain_2, trueKer_2, noiseVar)';...
%             GP_SE_Derivative([X(:,k-1);U(:,k-1)]', zTrain, yTrain_3, trueKer_3, noiseVar)';...
%             GP_SE_Derivative([X(:,k-1);U(:,k-1)]', zTrain, yTrain_4, trueKer_4, noiseVar)';...
%         ];
% 
%         A_k=M(:,1:4);
%         B_k=M(:,5:6);

        %Determine A_k and B_k: Homemade======================================


        M = [derivativeHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_1, noiseVar, sigma_f, sigma_l)';...
            derivativeHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_2, noiseVar, sigma_f, sigma_l)';...
            derivativeHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_3, noiseVar, sigma_f, sigma_l)';...
            derivativeHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_4, noiseVar, sigma_f, sigma_l)';...
        ];

        A_k=M(:,1:4);
        B_k=M(:,5:6);
        d_k=[meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_1, noiseVar, sigma_f, sigma_l);...
            meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_2, noiseVar, sigma_f, sigma_l);...
            meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_3, noiseVar, sigma_f, sigma_l);...
            meanHomemade([X(:,k-1);U(:,k-1)]', zTrain, yTrain_4, noiseVar, sigma_f, sigma_l)]...
            -[A_k , B_k ]*[X(:,k-1);U(:,k-1)];



        
        
        if norm(error_GP(2,k)) < 0.0001 % 0.01*norm(error_nominal(:,k)) %10 %*error_nominal(k)
            disp('Model validated')
            [U(:,k),infeasible]=ctrl{X(:,k),A_k,B_k,d_k,U(:,k-1)};
        else
            disp('Model NOT validated')
            [U(:,k),infeasible]=ctrl{X(:,k),A,B,zeros(4,1),U(:,k-1)};
        end     
    end
            
    disp(k)
    disp(infeasible)
    if infeasible ~= 0
        disp('Problem!!')
    end
    if(k>=k_fault)
        X(:,k+1)=A*X(:,k)+B*[U(1,k_fault);U(2,k)];
    else
        X(:,k+1)=A*X(:,k)+B*U(:,k);
    end
end


%%
%%PLOT :D
t=Ts*(0:stepNumber-1);
figure
subplot(2,2,1)
plot(t,X(1,:))
title('Pitch Angle rad/s')
xlabel('time s')
subplot(2,2,2)
plot(t,X(2,:))
title('Air speed m/s')
xlabel('time s')
subplot(2,2,3)
plot(t,X(3,:))
title('Pitch Rate rad')
xlabel('time s')
subplot(2,2,4)
plot(t,X(4,:))
title('AoA rad')
xlabel('time s')

figure
subplot(2,1,1)
plot(t(1:end-1),U(1,:))
title('Elevator deg')
xlabel('time s')
subplot(2,1,2)
plot(t(1:end-1),U(2,:))
title('THS deg ')
xlabel('time s')

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
title('Elevator slew rate def/s')
xlabel('time s')
subplot(2,1,2)
plot(t(1:end-2),dU(2,:))
title('THS slew rate deg/s')
xlabel('time s')

%%

figure
subplot(2,2,1)
plot(t,error_homemade(1,:));
title('error of the model, state 1')
legend('homemade')
xlabel('time s')
subplot(2,2,2)
plot(t,error_homemade(2,:));
title('error of the model, state 2')
legend('homemade')
xlabel('time s')
subplot(2,2,3)
plot(t,error_homemade(3,:));
title('error of the model, state 3')
legend('homemade')
xlabel('time s')
subplot(2,2,4)
plot(t,error_homemade(4,:));
title('error of the model, state 4')
legend('homemade')
xlabel('time s')

% figure
% subplot(2,2,1)
% plot(t,error_GPmat(1,:),t,error_GP(1,:),t,error_GPml(1,:),t,error_homemade(1,:));
% title('error of the model, state 1')
% legend('GPmat','GP','GPml','homemade')
% xlabel('time s')
% subplot(2,2,2)
% plot(t,error_GPmat(2,:),t,error_GP(2,:),t,error_GPml(2,:),t,error_homemade(2,:));
% title('error of the model, state 2')
% legend('GPmat','GP','GPml','homemade')
% xlabel('time s')
% subplot(2,2,3)
% plot(t,error_GPmat(3,:),t,error_GP(3,:),t,error_GPml(3,:),t,error_homemade(3,:));
% title('error of the model, state 3')
% legend('GPmat','GP','GPml','homemade')
% xlabel('time s')
% subplot(2,2,4)
% plot(t,error_GPmat(4,:),t,error_GP(4,:),t,error_GPml(4,:),t,error_homemade(4,:));
% title('error of the model, state 4')
% legend('GPmat','GP','GPml','homemade')
% xlabel('time s')