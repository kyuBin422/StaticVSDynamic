% clear the workplace
clear;clc;close all
addpath(genpath('utils'));
%%
% compute expected utility, ration, expected q, expected w with respect theta0
s0 =0.01;
gamma = 0.01;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0=0.6:0.01:3.5;

StaticQ=[];
DynamicQ=[];
StaticExpectedUtility=[];
DynamicExpectedUtility=[];
StaticRatio=[];
DynamicRatio=[];
StaticW=[];
DynamicW=[];
for i=theta0
    disp(i)
    optimalQ=solveFirstOrder(i,s0,gamma,kappa,c,F,R,G,b);
    % check the value information condition
    if valueInformation(i,optimalQ,s0,gamma,kappa,c,F,R,G,b)<0
        optimalQ=0;
    end
    StaticQ(end+1)=optimalQ;
    StaticW(end+1)=staticExpectedOmegaStar(i,optimalQ,gamma,s0,b,G);
    StaticExpectedUtility(end+1)=staticUtilityStar(i,optimalQ,s0,gamma,kappa,c,F,R,G,b);
    StaticRatio(end+1)=(c*optimalQ+F)/staticExpectedOmegaStar(i,optimalQ,gamma,s0,b,G);
    
    result=HMsimPnL(i,0,2,s0,gamma,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicQ(end+1)=mean(qList);
    DynamicW(end+1)=mean(omegaStar(thetaqList,qList,gamma,s0,b,G));
    DynamicExpectedUtility(end+1)=dynamicUtilityStar(thetaqList,qList,s0,gamma,kappa,c,F,R,G,b);
    DynamicRatio(end+1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,gamma,s0,b,G));
end

vector2File(theta0,StaticQ,DynamicQ,'\theta_0','E(q^*)','',false,"theta_0 vs Eq")

vector2File(theta0,StaticW,DynamicW,'\theta_0','E(q^*)','',false,"theta_0 vs omega")

vector2File(theta0,StaticExpectedUtility,DynamicExpectedUtility,'\theta_0','Expected Utility','',true,"theta_0 vs Expected Utility")

vector2File(theta0,StaticRatio,DynamicRatio,'\theta_0','cE(q^*)+F/E(\omega^*(q^*))','',true,"theta_0 vs ratio information")

%%
% compute expected utility, ration, expected q, expected w with respect gamma
s0 = 0.01;
gamma=0.006:0.001:0.014;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0_tmp1=1.15;
theta0_tmp2=0.95;

StaticQ=[];
DynamicQ=[];
StaticExpectedUtility=[];
DynamicExpectedUtility=[];
StaticRatio=[];
DynamicRatio=[];
StaticW=[];
DynamicW=[];
for i=gamma
    disp(i)
    % theta0=1.15
    optimalQ=solveFirstOrder(theta0_tmp1,s0,i,kappa,c,F,R,G,b);
    % check the value information condition
    if valueInformation(theta0_tmp1,optimalQ,s0,i,kappa,c,F,R,G,b)<0
        optimalQ=0;
    end
    StaticQ(end+1,1)=optimalQ;
    StaticW(end+1,1)=staticExpectedOmegaStar(theta0_tmp1,optimalQ,i,s0,b,G);
    StaticExpectedUtility(end+1,1)=staticUtilityStar(theta0_tmp1,optimalQ,s0,i,kappa,c,F,R,G,b);
    StaticRatio(end+1,1)=(c*optimalQ+F)/staticExpectedOmegaStar(theta0_tmp1,optimalQ,i,s0,b,G);
    % theta0=0.95
    optimalQ=solveFirstOrder(theta0_tmp2,s0,i,kappa,c,F,R,G,b);
    % check the value information condition
    if valueInformation(theta0_tmp2,optimalQ,s0,i,kappa,c,F,R,G,b)<0
        optimalQ=0;
    end
    StaticQ(end,2)=optimalQ;
    StaticW(end,2)=staticExpectedOmegaStar(theta0_tmp2,optimalQ,i,s0,b,G);
    StaticExpectedUtility(end,2)=staticUtilityStar(theta0_tmp2,optimalQ,s0,i,kappa,c,F,R,G,b);
    StaticRatio(end,2)=(c*optimalQ+F)/staticExpectedOmegaStar(theta0_tmp2,optimalQ,i,s0,b,G);
    
    % theta0=1.15
    result=HMsimPnL(theta0_tmp1,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicQ(end+1,1)=mean(qList);
    DynamicW(end+1,1)=mean(omegaStar(thetaqList,qList,i,s0,b,G));
    DynamicExpectedUtility(end+1,1)=dynamicUtilityStar(thetaqList,qList,s0,i,kappa,c,F,R,G,b);
    DynamicRatio(end+1,1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,i,s0,b,G));
    % theta0=0.95
    result=HMsimPnL(theta0_tmp2,0,2,s0,i,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicQ(end,2)=mean(qList);
    DynamicW(end,2)=mean(omegaStar(thetaqList,qList,i,s0,b,G));
    DynamicExpectedUtility(end,2)=dynamicUtilityStar(thetaqList,qList,s0,i,kappa,c,F,R,G,b);
    DynamicRatio(end,2)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,i,s0,b,G));
end
% plot theta0=1.15
vector2File(gamma,StaticQ(:,1),DynamicQ(:,1),'\gamma','E(q^*)','theta0=1.15',true,"gamma vs Eq theta0=115")

vector2File(gamma,StaticW(:,1),DynamicW(:,1),'\gamma','E(\omega^*(q^*))','theta0=1.15',true,"gamma vs omega theta0=115")

vector2File(gamma,StaticExpectedUtility(:,1),DynamicExpectedUtility(:,1),'\gamma','Expected Utility','theta0=1.15',false,"gamma vs Expected Utility theta0=115")

vector2File(gamma,StaticRatio(:,1),DynamicRatio(:,1),'\gamma','cE(q^*)+F/E(\omega^*(q^*))','theta0=1.15',false,"gamma vs ratio information theta0=115")

% plot theta0=0.95
vector2File(gamma,StaticQ(:,2),DynamicQ(:,2),'\gamma','E(q^*)','theta0=0.95',true,"gamma vs Eq theta0=095")

vector2File(gamma,StaticW(:,2),DynamicW(:,2),'\gamma','E(\omega^*(q^*))','theta0=1.95',true,"gamma vs omega theta0=095")

vector2File(gamma,StaticExpectedUtility(:,2),DynamicExpectedUtility(:,2),'\gamma','Expected Utility','theta0=1.95',false,"gamma vs Expected Utility theta0=095")

vector2File(gamma,StaticRatio(:,2),DynamicRatio(:,2),'\gamma','cE(q^*)+F/E(\omega^*(q^*))','theta0=1.95',false,"gamma vs ratio information theta0=095")

%%
% compute expected utility, ration, expected q, expected w with respect s0
s0 = 0.001:0.001:0.01;
gamma = 0.01;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
theta0_tmp1=1.15;
theta0_tmp2=0.95;

StaticQ=[];
DynamicQ=[];
StaticExpectedUtility=[];
DynamicExpectedUtility=[];
StaticRatio=[];
DynamicRatio=[];
StaticW=[];
DynamicW=[];
for i=s0
    disp(i)
    % theta0=1.15
    optimalQ=solveFirstOrder(theta0_tmp1,i,gamma,kappa,c,F,R,G,b);
    % check the value information condition
    if valueInformation(theta0_tmp1,optimalQ,i,gamma,kappa,c,F,R,G,b)<0
        optimalQ=0;
    end
    StaticQ(end+1,1)=optimalQ;
    StaticW(end+1,1)=staticExpectedOmegaStar(theta0_tmp1,optimalQ,gamma,i,b,G);
    StaticExpectedUtility(end+1,1)=staticUtilityStar(theta0_tmp1,optimalQ,i,gamma,kappa,c,F,R,G,b);
    StaticRatio(end+1,1)=(c*optimalQ+F)/staticExpectedOmegaStar(theta0_tmp1,optimalQ,gamma,i,b,G);
    % theta0=0.95
    optimalQ=solveFirstOrder(theta0_tmp2,i,gamma,kappa,c,F,R,G,b);
    % check the value information condition
    if valueInformation(theta0_tmp2,optimalQ,i,gamma,kappa,c,F,R,G,b)<0
        optimalQ=0;
    end
    StaticQ(end,2)=optimalQ;
    StaticW(end,2)=staticExpectedOmegaStar(theta0_tmp2,optimalQ,gamma,i,b,G);
    StaticExpectedUtility(end,2)=staticUtilityStar(theta0_tmp2,optimalQ,i,gamma,kappa,c,F,R,G,b);
    StaticRatio(end,2)=(c*optimalQ+F)/staticExpectedOmegaStar(theta0_tmp2,optimalQ,gamma,i,b,G);
    
    % theta0=1.15
    result=HMsimPnL(theta0_tmp1,0,2,i,gamma,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicQ(end+1,1)=mean(qList);
    DynamicW(end+1,1)=mean(omegaStar(thetaqList,qList,gamma,i,b,G));
    DynamicExpectedUtility(end+1,1)=dynamicUtilityStar(thetaqList,qList,i,gamma,kappa,c,F,R,G,b);
    DynamicRatio(end+1,1)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,gamma,i,b,G));
    % theta0=0.95
    result=HMsimPnL(theta0_tmp2,0,2,i,gamma,kappa,c,F,R,G,b);
    thetaqList=result(:,1);
    qList=result(:,2);
    DynamicQ(end,2)=mean(qList);
    DynamicW(end,2)=mean(omegaStar(thetaqList,qList,gamma,i,b,G));
    DynamicExpectedUtility(end,2)=dynamicUtilityStar(thetaqList,qList,i,gamma,kappa,c,F,R,G,b);
    DynamicRatio(end,2)=(c*mean(qList)+F)/mean(omegaStar(thetaqList,qList,gamma,i,b,G));
end
% plot theta0=1.15
vector2File(s0,StaticQ(:,1),DynamicQ(:,1),'s0','E(q^*)','theta0=1.15',true,"s0 vs Eq theta0=115")

vector2File(s0,StaticW(:,1),DynamicW(:,1),'s0','E(\omega^*(q^*))','theta0=1.15',true,"s0 vs omega theta0=115")

vector2File(s0,StaticExpectedUtility(:,1),DynamicExpectedUtility(:,1),'s0','Expected Utility','theta0=1.15',false,"s0 vs Expected Utility theta0=115")

vector2File(s0,StaticRatio(:,1),DynamicRatio(:,1),'s0','cE(q^*)+F/E(\omega^*(q^*))','theta0=1.15',false,"s0 vs ratio information theta0=115")

% plot theta0=0.95
vector2File(s0,StaticQ(:,2),DynamicQ(:,2),'s0','E(q^*)','theta0=0.95',true,"s0 vs Eq theta0=095")

vector2File(s0,StaticW(:,2),DynamicW(:,2),'s0','E(\omega^*(q^*))','theta0=0.95',true,"s0 vs omega theta0=095")

vector2File(s0,StaticExpectedUtility(:,2),DynamicExpectedUtility(:,2),'s0','Expected Utility','theta0=0.95',false,"s0 vs Expected Utility theta0=095")

vector2File(s0,StaticRatio(:,2),DynamicRatio(:,2),'s0','cE(q^*)+F/E(\omega^*(q^*))','theta0=0.05',false,"s0 vs ratio information theta0=095")

%%
function vector2File(Para,StaticList,DynamicList,Xlabel,Ylabel,Title,LegendPosition,FileName)
fig=figure;
plot(Para,StaticList,'--')
hold on
plot(Para,DynamicList,'-.')
xlabel(Xlabel,'FontSize',12,'FontWeight','bold')
ylabel(Ylabel,'FontWeight','bold','FontSize',12)
title(Title,'FontSize',12,'FontWeight','bold')
if LegendPosition
    legend("static","dynamic",'FontSize',12)
else
    legend("static","dynamic",'FontSize',12,'Location','southeast')
end
savefig(fig,'image/'+FileName+'.fig')
saveas(fig,'image/imageEPS/'+FileName,'epsc')
close(fig)
end