function dynamicUtility = HMsimPnL(X0,p)
addpath(genpath('utils'));
% X0 is the  intial belief
% p = 1, plot the stopping boundaries, else plot only Expected pnl

s0 = 0.01;
gamma = 0.01;
kappa = 100;
c = 0.01;
F = 5;
R = exp(0.05);
G = -0.5;
b = 0.5;
T = s0; 
% initial q and thetaq
qList=[];
thetaqList=[];


rounds = 10000;

X = 1; %this has to correspond to the paramenters in HMdynamicFD
[~, ~,M]=HMdynamicFD(X, s0, gamma, kappa, c, F, R, G, b);


if p == 1
    hold on
else
    close(1)
end

%initialzie pnl
exp_PNL = zeros(1,100);
exp_PNL2 = zeros(1,100);


for i = 1:rounds
    
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
    qList(end+1)=q;
    thetaqList(end+1)=thetaq;
    
%     % discretize Shock 1 for plotting
%     S1 = linspace(0.6,1.05,100);
%     S2 = linspace(1,4,100);
%     
%     for j = 1:100
%         exp_PNL(j) = exp_PNL(j) + S1pnl(S1(j),q,thetaq,s0,gamma,kappa,c,F,R,G,b);
%         exp_PNL2(j) = exp_PNL2(j) + S2pnl(S2(j),q,thetaq,s0,gamma,kappa,c,F,R,G,b,X0);
%     end
    
end

% exp_PNL = exp_PNL/rounds;
% exp_PNL2 = exp_PNL2/rounds;
% 
% figure (2)
% plot(S1, exp_PNL)
% %plot(S1, exp_PNL, '--')
% 
% figure (3)
% plot(S2, exp_PNL2)
% %plot(S2, exp_PNL2, '--')

figure (4)
scatter(qList,thetaqList)
xlabel('$q$','Interpreter','latex');
ylabel('$\hat{\theta}_q$','Interpreter','latex');
% ultity function

if p == 1
    hold on
else
    close(4)
end
dynamicUtility=dynamicUtilityStar(thetaqList,qList);
% save('matlab_x0=1point3.mat','qList','thetaqList')
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

