load('trainData.mat')
A=[0.86719   ,6.6936*10^(-5)     ,-0.19095  ,0; ...
    -0.027773, 0.99895           , 0.89264, -1.9609; ...
      0.20146, -2.1676*10^(-4)     , 0.88379,  0; ....
      0.2    ,  0                  ,  0      , 1];
B=[-3.7758 *10^(-3)  , -9.0408*10^(-3);...
    0                ,0;...
      -1.2629*10^(-3), -3.2794*10^(-4);...
     0              ,0];
 
 

n=700;
inver_wid=linspace(1/88,88,n);
error_GP=zeros(1,n);
noiseVar=0;




for i=1:n
    
    trueKer_1 = kernCreate(zTrain, 'rbf');
    trueKer_1.inverseWidth = 1/inver_wid(i);
    %trueKer_1.variance = 39.49;

    trueKer_2 = kernCreate(zTrain, 'rbf');
    trueKer_2.inverseWidth = 1/inver_wid(i);
    %trueKer_2.variance = 39.49;

    trueKer_3 = kernCreate(zTrain, 'rbf');
    trueKer_3.inverseWidth = 1/inver_wid(i);
    %trueKer_3.variance = 39.49;   

    trueKer_4 = kernCreate(zTrain, 'rbf');
    trueKer_4.inverseWidth = 1/inver_wid(i);
    %trueKer_4.variance = 39.49;  
    
    
    A_k = [GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_1, trueKer_1, noiseVar);... 
        GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_2, trueKer_2, noiseVar);... 
        GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_3, trueKer_3, noiseVar);... 
        GP_SE_mean_var([1,0,0,0,0,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,1,0,0,0,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,0,1,0,0,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,0,0,1,0,0], zTrain, yTrain_4, trueKer_4, noiseVar)];
    

    B_k=[GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_1, trueKer_1, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_1, trueKer_1, noiseVar);...
        GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_2, trueKer_2, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_2, trueKer_2, noiseVar);... 
        GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_3, trueKer_3, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_3, trueKer_3, noiseVar);... 
        GP_SE_mean_var([0,0,0,0,1,0], zTrain, yTrain_4, trueKer_4, noiseVar), GP_SE_mean_var([0,0,0,0,0,1], zTrain, yTrain_4, trueKer_4, noiseVar)];
        
    error_GP(i) = mean(vecnorm( X(:,2:end)-A_k*X(:,1:end-1)-B_k*U(:,1:end) ));

end

error_nominal=mean(vecnorm( X(:,2:end)-A*X(:,1:end-1)-B*U(:,1:end) ));

plot(inver_wid,error_GP)
min(error_GP)