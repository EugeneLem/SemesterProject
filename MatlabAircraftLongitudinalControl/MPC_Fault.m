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
C_u_delta=[-1 0 1 0;...
            0 -1 0 1;...
           1   0 -1 0;...
           0   1  0  -1];
c_u_delta=[37;0.236;37;0.236]/Ts;%Need to check is elevator and THS are not inverted


%%MPC!!!!!!!!
x=sdpvar(4,horizon);
u=sdpvar(2,horizon);
ops = sdpsettings('verbose',1);
con=[];
obj = 0;
for i = 1:horizon-1 
    con = [con , x(:,i+1) == A*x(:,i) + B*u(:,i)]; % System dynamics
    con = [con , C_u*u(:,i) <= c_u];                 % Input constraints
    con = [con , C_u_delta*[u(:,i);u(:,i+1)] <= c_u_delta];  % Input change constraints
    %obj = obj + x(:,i)'*Q*x(:,i)  + (u(:,i+1)-u(:,i))'*R_delta*(u(:,i+1)-u(:,i));  % Cost function
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i) + (u(:,i+1)-u(:,i))'*R_delta*(u(:,i+1)-u(:,i));  % Cost function

end
obj=obj+[u(:,horizon);x(:,horizon)]'*P_k*[u(:,horizon);x(:,horizon)];
ctrl = optimizer(con, obj,ops, x(:,1), u(:,1));


%%Let's simulate our system
stepNumber=100;
k_fault=2;
X=zeros(4,stepNumber); 
U=zeros(2,stepNumber); 
X(:,1)=[1;1;1;1];

for k=1:stepNumber-1                    
    U(:,k)=ctrl{X(:,k)};
    if(k>=k_fault)
        U(1,k)=U(1,k-1);
    end
    X(:,k+1)=A*X(:,k)+B*U(:,k);
end
U(:,stepNumber)=ctrl{X(:,k)};
U(1,stepNumber)=U(1,stepNumber-1);

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
plot(t,U(1,:))
title('Elevator')
subplot(2,1,2)
plot(t,U(2,:))
title('THS')