function result = HMsimPnL(X0,p,OutputBool,s0,gamma,kappa,c,F,R,G,b)
addpath(genpath('utils'));
close(figure(1));close(figure(2))
% X0 is the  intial belief
% p = 1, plot the stopping boundaries, else plot only Expected pnl

% s0 = 0.01;
% gamma = 0.01;
% kappa = 100;
% c = 0.01;
% F = 5;
% R = exp(0.05);
% G = -0.5;
% b = 0.5;

T = s0;


rounds = 10000;

X = 1; %this has to correspond to the paramenters in HMdynamicFD
[~, ~,M]=HMdynamicFD(X, s0, gamma, kappa, c, F, R, G, b);

% initial q and thetaq
qList=zeros(rounds,1);
thetaqList=zeros(rounds,1);
% check the condition output figure or not
if p == 1
    hold on
else
    close(1)
end


parfor i = 1:rounds
    
    % create brownian motion with drift = 0, volatility 1
    obj = bm(0,1);
    
    dt = T/1000;
    
    % simulate a sample path for 1000 periods, with dt
    [X1, T1] = simulate(obj, 1000, 'DeltaTime', dt);
    
    % transform paths to plot on other figure axis..
    %s0 = 0.01;
    m=999;
    
    Xmin = X-3*sqrt(s0);  % Minimum theta
    Xmax = X+12*sqrt(s0); % Maximum theta
    dX = (Xmax-Xmin)/(m+1);
    
    %initial belief hat{theta}_0
    %X0 = 1;
    X1 = X1+X0;
    
    
    X1plot = (X1 - Xmin)/dX;
    T1plot = T1/dt+1;
    
    
    % check at which index it hits boundary
    stopping_indx = HMstop(M, T1plot, X1plot);
    
    if p == 1
        plot(T1plot(1:stopping_indx), X1plot(1:stopping_indx));
    end
    %     test1(i)=T1plot(stopping_indx);
    %         test2(i)=X1plot(stopping_indx);
    
    % get actual T = q and X = theta_q values from the stopping_indx
    thetaq = X1(stopping_indx);
    tau = T1(stopping_indx);
    q = tau/(s0^2 - s0*tau);
    qList(i)=q;
    thetaqList(i)=thetaq;
    
end

% meet the value information must be more than one
% and delete the outlier
DeleteFlag= fullInformation(X0,s0,gamma,F,R,G,b)>=0 | dynamicUtilityStar(thetaqList,qList,s0,gamma,kappa,c,F,R,G,b)<0;
thetaqList(DeleteFlag)=[];
qList(DeleteFlag)=[];


if p == 1
    title("X0="+string(X0)+" "+"E(q)="+mean(qList)+" "+"$E(\hat{\theta}_q)=$"+mean(thetaqList),'Interpreter','latex','FontSize',12,'FontWeight','bold')
    figure (2)
    scatter(qList,thetaqList)
    xlabel('$q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
    ylabel('$\hat{\theta}_q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
    title("X0="+string(X0)+" "+"E(q)="+mean(qList)+" "+"$E(\hat{\theta}_q)=$"+mean(thetaqList),'Interpreter','latex','FontSize',12,'FontWeight','bold')
end
result=[];
switch OutputBool
    case 0
        result=thetaqList;
    case 1
        result=qList;
    case 2
        result=[thetaqList,qList];
end
% save('matlab_x0=1point3.mat','qList','thetaqList','-v6')
end

