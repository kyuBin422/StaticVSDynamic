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
    i
    % static
    disp(i)
    optimalQ=solveFirstOrder(i,s0,gamma,kappa,c,F,R,G,b);
    euStatic(end+1)=staticUtilityStar(i,optimalQ,s0,gamma,kappa,c,F,R,G,b);
    
    % dynamic
    result=HMsimPnL(i,0,2,s0,gamma,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    euDynamic(end+1)=dynamicUtilityStar(thetaqList,qList,gamma,kappa,c,F,R,G,b);
    
end

figure(1);
plot(theta0,euStatic)
hold on
plot(theta0,euDynamic)
xlabel("theta0")
ylabel("expected utility")
legend("static","dynamic")
%%
% ration of information cost when theta0=1.15 and theta0=0.95
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
StaticRation_theta115=[];
DynamicRation_theta115=[];
StaticRation_theta095=[];
DynamicRation_theta095=[];

for i=gamma
    i
    % theta0=1.15
    optimalQ=solveFirstOrder(theta0_1,s0,i,kappa,c,F,R,G,b);
    % check the value information condition
    if (valueInformation(theta0_1,optimalQ,s0,i,kappa,c,F,R,G,b)<0)
        optimalQ=0;
    end
    StaticQ=optimalQ;
    StaticRation_theta115(end+1)=(c*optimalQ+F)/staticExpectedOmegaStar(theta0_1,optimalQ,i,s0,b,G);
    
    result=HMsimPnL(theta0_1,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicRation_theta115(end+1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,i,s0,b,G));
    
    % theta0=0.95
    optimalQ=solveFirstOrder(theta0_2,s0,i,kappa,c,F,R,G,b);
    % check the value information condition
    if (valueInformation(theta0_2,optimalQ,s0,i,kappa,c,F,R,G,b)<0)
        optimalQ=0;
    end
    StaticQ=optimalQ;
    StaticRation_theta095(end+1)=(c*optimalQ+F)/staticExpectedOmegaStar(theta0_2,optimalQ,i,s0,b,G);
    
    result=HMsimPnL(theta0_2,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicRation_theta095(end+1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,i,s0,b,G));
    
end
figure;
subplot(211)
plot(gamma,StaticRation_theta115)
hold on
plot(gamma,DynamicRation_theta115)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('cE(q^*)+F/E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=1.15>-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

subplot(212)
plot(gamma,StaticRation_theta095)
hold on
plot(gamma,DynamicRation_theta095)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('cE(q^*)+F/E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=0.95<-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)
%%
% s0 with repsect to the optimal q and exptect w
s0 = 0.001:0.0001:0.01;
gamma=0.01;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0_1=1.15;
theta0_2=0.95;

% initial list to restore variable
StaticQ_theta115=[];
DynamicQ_theta115=[];
StaticQ_theta095=[];
DynamicQ_theta095=[];

StaticW_theta115=[];
DynamicW_theta115=[];
StaticW_theta095=[];
DynamicW_theta095=[];


for i =s0
    i
    % theta0=1.15
    [ StaticQ_theta115(end+1),StaticW_theta115(end+1),DynamicQ_theta115(end+1),DynamicW_theta115(end+1)]=...
        ComputeStaticVSDynamic(theta0_1,i,gamma,kappa,c,F,R,G,b);
    
    [ StaticQ_theta095(end+1),StaticW_theta095(end+1),DynamicQ_theta095(end+1),DynamicW_theta095(end+1)]=...
        ComputeStaticVSDynamic(theta0_2,i,gamma,kappa,c,F,R,G,b);
end
figure;
subplot(211)
plot(s0,StaticQ_theta115)
hold on
plot(s0,DynamicQ_theta115)
xlabel('s0','FontWeight','bold','FontSize',16)
ylabel('E(q^*)','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=1.15>-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

subplot(212)
plot(s0,StaticQ_theta095)
hold on
plot(s0,DynamicQ_theta095)
xlabel('s0','FontWeight','bold','FontSize',16)
ylabel('E(q^*)','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=0.95<-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

figure;
subplot(211)
plot(s0,StaticW_theta115)
hold on
plot(s0,DynamicW_theta115)
xlabel('s0','FontWeight','bold','FontSize',16)
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=1.15>-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

subplot(212)
plot(s0,StaticW_theta095)
hold on
plot(s0,DynamicW_theta095)
xlabel('s0','FontWeight','bold','FontSize',16)
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=0.95<-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

%%
% gamma with repsect to the optimal q and exptect w
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

% initial list to restore variable
StaticQ_theta115=[];
DynamicQ_theta115=[];
StaticQ_theta095=[];
DynamicQ_theta095=[];

StaticW_theta115=[];
DynamicW_theta115=[];
StaticW_theta095=[];
DynamicW_theta095=[];


for i =gamma
    i
    % theta0=1.15
    [ StaticQ_theta115(end+1),StaticW_theta115(end+1),DynamicQ_theta115(end+1),DynamicW_theta115(end+1)]=...
        ComputeStaticVSDynamic(theta0_1,s0,i,kappa,c,F,R,G,b);
    
    [ StaticQ_theta095(end+1),StaticW_theta095(end+1),DynamicQ_theta095(end+1),DynamicW_theta095(end+1)]=...
        ComputeStaticVSDynamic(theta0_2,s0,i,kappa,c,F,R,G,b);
end
figure;
subplot(211)
plot(gamma,StaticQ_theta115)
hold on
plot(gamma,DynamicQ_theta115)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('E(q^*)','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=1.15>-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

subplot(212)
plot(gamma,StaticQ_theta095)
hold on
plot(gamma,DynamicQ_theta095)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('E(q^*)','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=0.95<-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

figure;
subplot(211)
plot(gamma,StaticW_theta115)
hold on
plot(gamma,DynamicW_theta115)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=1.15>-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

subplot(212)
plot(gamma,StaticW_theta095)
hold on
plot(gamma,DynamicW_theta095)
xlabel('\gamma','FontWeight','bold','FontSize',16)
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',16)
title('$\hat{\theta}_0=0.95<-G/b$','Interpreter','latex','FontWeight','bold','FontSize',16)
legend("static","dynamic",'FontSize',13)

%%
function [StaticQ,StaticW,DynamicQ,DynamicW]=ComputeStaticVSDynamic(theta0,s0,gamma,kappa,c,F,R,G,b)
optimalQ=solveFirstOrder(theta0,s0,gamma,kappa,c,F,R,G,b);
% check the value information condition
if (valueInformation(theta0,optimalQ,s0,gamma,kappa,c,F,R,G,b)<0)
    optimalQ=0;
end

StaticQ=optimalQ;
StaticW=staticExpectedOmegaStar(theta0,optimalQ,gamma,s0,b,G);

result=HMsimPnL(theta0,0,2,s0,gamma,kappa,c,F,R,G,b);
thetaqList=result(:,1);
qList=result(:,2);
DynamicQ=mean(qList);
DynamicW=mean(omegaStar(thetaqList,qList,gamma,s0,b,G));

end
