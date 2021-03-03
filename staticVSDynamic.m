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
figure;
subplot(221)
plot(theta0,StaticQ,'--')
hold on
plot(theta0,DynamicQ,'-.')
ylabel('E(q^*)','FontWeight','bold','FontSize',12)

subplot(222)
plot(theta0,StaticW,'--')
hold on
plot(theta0,DynamicW,'-.')
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',12)

subplot(223)
plot(theta0,StaticExpectedUtility,'--')
hold on
plot(theta0,DynamicExpectedUtility,'-.')
ylabel("Expected Utility",'FontSize',12,'FontWeight','bold')

subplot(224)
plot(theta0,StaticRatio,'--')
hold on
plot(theta0,DynamicRatio,'-.')
xlabel("\theta0",'FontSize',12,'FontWeight','bold')
ylabel('cE(q^*)+F/E(\omega^*(q^*))','FontWeight','bold','FontSize',12)

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,"\theta0",'FontSize',12,'FontWeight','bold')

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
fig=figure;
subplot(221)
plot(gamma,StaticQ(:,1),'LineStyle','--')
hold on
plot(gamma,DynamicQ(:,1),'LineStyle','-.')
ylabel('E(q^*)','FontWeight','bold','FontSize',12)

subplot(222)
plot(gamma,StaticW(:,1),'LineStyle','--')
hold on
plot(gamma,DynamicW(:,1),'LineStyle','-.')
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',12)

subplot(223)
plot(gamma,StaticExpectedUtility(:,1),'LineStyle','--')
hold on
plot(gamma,DynamicExpectedUtility(:,1),'LineStyle','-.')
ylabel("Expected Utility",'FontSize',12,'FontWeight','bold')

subplot(224)
plot(gamma,StaticRatio(:,1),'LineStyle','--')
hold on
plot(gamma,DynamicRatio(:,1),'LineStyle','-.')
ylabel("cE(q^*)+F/E(\omega^*(q^*))",'FontSize',12,'FontWeight','bold')
legend("static","dynamic",'FontSize',12)

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'\gamma','FontSize',12,'FontWeight','bold')
title(han,'theta0=1.15','FontWeight','bold','FontSize',12)

% plot theta0=1.15
fig=figure;
subplot(221)
plot(gamma,StaticQ(:,2),'LineStyle','--')
hold on
plot(gamma,DynamicQ(:,2),'LineStyle','-.')
ylabel('E(q^*)','FontWeight','bold','FontSize',12)

subplot(222)
plot(gamma,StaticW(:,2),'LineStyle','--')
hold on
plot(gamma,DynamicW(:,2),'LineStyle','-.')
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',12)

subplot(223)
plot(gamma,StaticExpectedUtility(:,2),'LineStyle','--')
hold on
plot(gamma,DynamicExpectedUtility(:,2),'LineStyle','-.')
ylabel("Expected Utility",'FontSize',12,'FontWeight','bold')

subplot(224)
plot(gamma,StaticRatio(:,2),'LineStyle','--')
hold on
plot(gamma,DynamicRatio(:,2),'LineStyle','-.')
ylabel("cE(q^*)+F/E(\omega^*(q^*))",'FontSize',12,'FontWeight','bold')
legend("static","dynamic",'FontSize',12)

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'\gamma','FontSize',12,'FontWeight','bold')
title('theta0=0.95','FontWeight','bold','FontSize',12)
%%
% compute expected utility, ration, expected q, expected w with respect s0
s0 = 0.001:0.00005:0.01;
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
fig=figure;
subplot(221)
plot(s0,StaticQ(:,1),'LineStyle','--')
hold on
plot(s0,DynamicQ(:,1),'LineStyle','-.')
ylabel('E(q^*)','FontWeight','bold','FontSize',12)

subplot(222)
plot(s0,StaticW(:,1),'LineStyle','--')
hold on
plot(s0,DynamicW(:,1),'LineStyle','-.')
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',12)

subplot(223)
plot(s0,StaticExpectedUtility(:,1),'LineStyle','--')
hold on
plot(s0,DynamicExpectedUtility(:,1),'LineStyle','-.')
ylabel("Expected Utility",'FontSize',12,'FontWeight','bold')

subplot(224)
plot(s0,StaticRatio(:,1),'LineStyle','--')
hold on
plot(s0,DynamicRatio(:,1),'LineStyle','-.')
ylabel("cE(q^*)+F/E(\omega^*(q^*))",'FontSize',12,'FontWeight','bold')
legend("static","dynamic",'FontSize',12)

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'s0','FontSize',12,'FontWeight','bold')
title(han,'theta0=1.15','FontWeight','bold','FontSize',12)

% plot theta0=0.95
fig=figure;
subplot(221)
plot(s0,StaticQ(:,2),'LineStyle','--')
hold on
plot(s0,DynamicQ(:,2),'LineStyle','-.')
ylabel('E(q^*)','FontWeight','bold','FontSize',12)

subplot(222)
plot(s0,StaticW(:,2),'LineStyle','--')
hold on
plot(s0,DynamicW(:,2),'LineStyle','-.')
ylabel('E(\omega^*(q^*))','FontWeight','bold','FontSize',12)

subplot(223)
plot(s0,StaticExpectedUtility(:,2),'LineStyle','--')
hold on
plot(s0,DynamicExpectedUtility(:,2),'LineStyle','-.')
ylabel("Expected Utility",'FontSize',12,'FontWeight','bold')

subplot(224)
plot(s0,StaticRatio(:,2),'LineStyle','--')
hold on
plot(s0,DynamicRatio(:,2),'LineStyle','-.')
ylabel("cE(q^*)+F/E(\omega^*(q^*))",'FontSize',12,'FontWeight','bold')
legend("static","dynamic",'FontSize',12)

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'s0','FontSize',12,'FontWeight','bold')
title('theta0=0.95','FontWeight','bold','FontSize',12)
