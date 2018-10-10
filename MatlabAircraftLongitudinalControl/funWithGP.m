rng(0,'twister'); % For reproducibility
n = 100;
x = linspace(-10,10,n)';
y = 1 + x*5e-2 + sin(x)./x + 0.01*randn(n,1);


gprMdl = fitrgp(x,y)%,'Basis','linear','FitMethod','exact','PredictMethod','exact');

% Xnew=linspace(-20,0,n)';
% [ypred,~,yci] = predict(gprMdl,x);
% 
% plot(x,y,'b.');
% hold on;
% plot(x,ypred,'r','LineWidth',1.5);
% plot(x,yci(:,1),'k--');
% plot(x,yci(:,2),'k--');
% xlabel('x');
% ylabel('y');
% legend('Data','GPR predictions');
% hold off

Xnew=linspace(-20,20,n)';
[ypred,~,yci] = predict(gprMdl,Xnew);

plot(x,y,'b.');
hold on;
plot(Xnew,ypred,'r','LineWidth',1.5);
plot(Xnew,yci(:,1),'k--');
plot(Xnew,yci(:,2),'k--');
xlabel('x');
ylabel('y');
legend('Data','GPR predictions');
hold off


