clear;clc;close all
addpath(genpath('utils')); 
%%
% compute expected utility, ration, expected q, expected w with respect theta0
s0 =linspace(0.001,0.01,100);
gamma = 0.01;
kappa = 1000;
c = 0.01;
F = 100;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0=1.06;

StaticQ=[];
DynamicQ=[];
StaticW=[];
DynamicW=[];

para=[];
for i=s0
    disp(i)
    for j=gamma
    optimalQ=solveFirstOrder(theta0,i,j,kappa,c,F,R,G,b);
    % check the value information condition
    if valueInformation(theta0,optimalQ,i,j,kappa,c,F,R,G,b)<0
        optimalQ=0;
    end
    StaticQ(end+1)=optimalQ;
    StaticW(end+1)=staticExpectedOmegaStar(theta0,optimalQ,j,i,b,G);

    result=HMsimPnL(theta0,0,2,i,j,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicQ(end+1)=mean(qList);
    DynamicW(end+1)=mean(omegaStar(thetaqList,qList,j,i,b,G));
    
    para(end+1,:)=[i,j];
    end
end