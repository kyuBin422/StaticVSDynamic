function result = HMsimPnL(X0,p,OutputBool,s0,gamma,kappa,c,F,R,G,b)
addpath(genpath('utils'));
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
    
    % get actual T = q and X = theta_q values from the stopping_indx
    thetaq = X1(stopping_indx);
    tau = T1(stopping_indx);
    q = tau/(s0^2 - s0*tau);
    qList(i)=q;
    thetaqList(i)=thetaq;
    
end

% meet the value information must be more than one

qList(fullInformation(X0,s0,gamma,F,R,G,b)>=0)=0;
qList(dynamicUtilityStar(thetaqList,qList,s0,gamma,kappa,c,F,R,G,b)<=0)=0;


if p == 1
    figure (4)
    scatter(qList,thetaqList)
    xlabel('$q$','Interpreter','latex');
    ylabel('$\hat{\theta}_q$','Interpreter','latex');
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
save('matlab_x0=1point3.mat','qList','thetaqList')
end





function PNL = S1pnl(S1,q,X,s0,gamma,kappa,c,F,R,G,b)

w = max( (G + b*X)/( gamma*(b^2)*( s0/(1+s0*q) ) ),0 );

if q <= 0
    F = 0;
end

PNL = (G + b*S1)* w + (kappa- c*q - F  )*(R)-kappa;

end


function PNL = S2pnl(S2,q,X,s0,gamma,kappa,c,F,R,G,b,theta0)

w = max( (G + b*X)/( gamma*(b^2)*( s0/(1+s0*q) ) ),0 );

if q <= 0
    F = 0;
end

PNL = (G + b*(theta0 - S2* sqrt(s0/(1+s0*q) ) ) )* w + (kappa- c*q - F  )*(R)-kappa;

end

