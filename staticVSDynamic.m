% clear the workplace
clear;clc;close all
addpath(genpath('utils'));
%%
% compute the expected utility
euStatic=[];
euDynamic=[];


% initla value
s0 = 0.01;
gamma = 0.01;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0=0.6:0.01:1.5;


for i=theta0
    % static
    disp(i)
    optimalQ=solveFirstOrder(i,s0,gamma,kappa,c,F,R,G,b);
    euStatic(end+1)=staticUtilityStar(i,optimalQ,s0,gamma,kappa,c,F,R,G,b);
    
    % dynamic
    result=HMsimPnL(i,0,2,s0,gamma,kappa,c,F,R,G,b);
    result=HMsimPnL(theta0,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(1,:);
    qList=result(2,:);
    euDynamic(end+1)=dynamicUtilityStar(thetaqList,qList);
    
end

figure(1);
plot(theta0,euStatic)
hold on
plot(theta0,euDynamic)
xlabel("theta0")
ylabel("expected utility")
legend("static","dynamic")
%%
% ration of information cost when theta0=1.15
% initial parameter value
s0 = 0.01;
gamma=0.006:0.001:0.014;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0_1=1.15;
theta0_2=0.95;

% store the ration information cost
st1=[];
dy1=[];
st2=[];
dy2=[];

for i=gamma
    % theta0=1.15
    optimalQ=solveFirstOrder(theta0_1,s0,i,kappa,c,F,R,G,b);
    st1(end+1)=(c*optimalQ+F)/omegaStar(theta0_1,optimalQ,i,s0,b,G);
    
    result=HMsimPnL(theta0_1,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(1,:);
    qList=result(2,:);
    dy1(end+1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,i,s0,b,G));
    
    % theta0=0.95
    optimalQ=solveFirstOrder(theta0_2,s0,i,kappa,c,F,R,G,b);
    st2(end+1)=(c*optimalQ+F)/omegaStar(theta0_2,optimalQ,i,s0,b,G);
    
    result=HMsimPnL(theta0_2,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(1,:);
    qList=result(2,:);
    dy2(end+1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,i,s0,b,G));
    
end
figure;
subplot(211)
plot(gamma,st1)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('cq^*/E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
hold on
plot(gamma,dy1)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('cq^*/E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_q=1.15>-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)


subplot(212)
plot(gamma,st2)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('cq^*/E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
hold on
plot(gamma,dy2)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('cq^*/E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_q=0.95<-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)


