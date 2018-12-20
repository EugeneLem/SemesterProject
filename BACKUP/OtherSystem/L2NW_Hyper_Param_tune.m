clear all;
load('data_for_opti.mat') %Data from a simulation with no GP

param_min = [0, 0.1];
param_max = [6 , 3];
dh = 0.5;
dlambda = 0.1;

h = param_min(1):dh:param_max(1);
lambda = param_min(2):dlambda:param_max(2);

[H,LAMBDA] = meshgrid(h,lambda);
error1 = zeros(size(H));
error2 = zeros(size(H));
error3 = zeros(size(H));
error4 = zeros(size(H));

for i = 1:length(h)
    for j = 1:length(lambda)
        for k = 1:size(Z,1)
            error1(j,i) = error1(j,i) + abs(Y_1(k)- L2NW_mean(Z(k,:),Z([1:k-1, k+1:end],:),Y_1([1:k-1, k+1:end],:),[h(i),lambda(j)],1));
        end
    end
    i
end
figure
surf(H,LAMBDA,error1)
title('state 1')
xlabel('h')
ylabel('lambda')
 %---------------------------------------------------------------------------------------------------
 for i = 1:length(h)
    for j = 1:length(lambda)
        for k = 1:size(Z,1)
            error2(j,i) = error2(j,i) + abs(Y_2(k)- L2NW_mean(Z(k,:),Z([1:k-1, k+1:end],:),Y_2([1:k-1, k+1:end],:),[h(i),lambda(j)],1));
        end
    end
    i
end
figure
surf(H,LAMBDA,error2)
title('state 2')
xlabel('h')
ylabel('lambda')
%----------------------------------------------------- 
 for i = 1:length(h)
    for j = 1:length(lambda)
        for k = 1:size(Z,1)
            error3(j,i) = error3(j,i) + abs(Y_3(k)- L2NW_mean(Z(k,:),Z([1:k-1, k+1:end],:),Y_3([1:k-1, k+1:end],:),[h(i),lambda(j)],1));
        end
    end
    i
end
figure
surf(H,LAMBDA,error3)
title('state 3')
xlabel('h')
ylabel('lambda')
 %----------------------------------------------------------------------------
 for i = 1:length(h)
    for j = 1:length(lambda)
        for k = 1:size(Z,1)
            error4(j,i) = error4(j,i) + abs(Y_4(k)- L2NW_mean(Z(k,:),Z([1:k-1, k+1:end],:),Y_4([1:k-1, k+1:end],:),[h(i),lambda(j)],1));
        end
    end
    i
end
figure
surf(H,LAMBDA,error4)
title('state 4')
xlabel('h')
ylabel('lambda')