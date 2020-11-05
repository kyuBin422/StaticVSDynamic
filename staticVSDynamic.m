% clear the workplace
clear;clc;close all
addpath(genpath('utils'));
format long

%%
% store the expected utility
euStatic=[];
euDynamic=[];


% initla value
s0=0.01;
theta0=0.6:0.01:1.5;


for i=theta0
    % static
    disp(i)
    optimalQ=solveFirstOrder(i);
    euStatic(end+1)=staticUtilityStar(i,optimalQ);
    
    % dynamic
    euDynamic(end+1)=HMsimPnL(i,0);

end

figure(1);
plot(theta0,euStatic)
hold on
plot(theta0,euDynamic)
xlabel("theta0")
ylabel("expected utility")
legend("static","dynamic")

