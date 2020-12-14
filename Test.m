clc;clear;close all
load matlab_x0=1point3.mat
addpath(genpath('utils'));
addpath(genpath('FlexDist'));


x0=[-0.8093787,0.5,0.5,930.9278351,601.8286083 2873199923];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
up=[];
nonlcon=[];
GHLike(qList,x0)

options = optimoptions(@fmincon,'Display','iter','PlotFcns',@optimplotfval);
params=fmincon(@(x) GHLike(qList,x),x0,A,b,Aeq,beq,lb,up,nonlcon,options);

x=linspace(min(qList),max(qList),length(qList));
y=generalizedHyperbolicDistrbution(x,params(1),params(2),params(3),params(4),params(5),params(6));
figure;
[counts,centers] = hist(thetaqList,100);
bar(centers,counts/trapz(centers,counts))
hold on
plot(x,y,'LineStyle','-','Color','r')
xlabel('$\hat{\theta}_q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
ylabel("density",'FontSize',12,'FontWeight','bold')


function LogL = GHLike(x,Params)

lambda=Params(1);
chi=Params(2);
psi=Params(3);
mu=Params(4);
% here 'Sigma'='sigma'^2
sigma=Params(5);
gamma=Params(6);

if lambda >0
    chi(chi<0)=0;
    psi(chi<=0)=1e-11;
elseif lambda ==0
    chi(chi<=0)=1e-11;
    psi(psi<=0)=1e-11;
elseif lambda <0
    chi(chi<=0)=1e-11;
    psi(psi<0)=1e-11;
end

a= -lambda/2.*log(chi.*psi)+lambda.*log(psi)+(1/2-lambda).*log(psi+(gamma.^2./sigma))-...
    1/2*log(2*pi)-log(sigma)/2-log(besselk(lambda,sqrt(psi.*chi)));

% compute the generalized hyperbolic distrbution
LogL=a+log(besselk(lambda-0.5,sqrt((chi+(x-mu).^2./sigma)...
    .*(psi+gamma.^2./sigma))))+gamma.*(x-mu)./sigma ...
    +(lambda./2-1/4)*(log(chi+(x-mu).^2./sigma)+log(psi+gamma.^2./sigma));

LogL=-sum(LogL);

LogL(~isfinite(LogL))  =  1.0e-20;
LogL(~(~imag(LogL)))   =  1.0e-20;

end























% close all;
% 
% 
% HMsimPnL(0.9,1,2,0.01,0.01,100,0.01,5,exp(0.05),-0.5,0.5)
% figure(3);
% fig3=subplot(221);
% copyobj(allchild(get(figure(1),'CurrentAxes')),fig3)
% figure(4);
% fig3=subplot(221);
% copyobj(allchild(get(figure(2),'CurrentAxes')),fig3)
% 
% 
% HMsimPnL(0.95,1,2,0.01,0.01,100,0.01,5,exp(0.05),-0.5,0.5)
% figure(3);
% fig3=subplot(222);
% copyobj(allchild(get(figure(1),'CurrentAxes')),fig3)
% figure(4);
% fig3=subplot(222);
% copyobj(allchild(get(figure(2),'CurrentAxes')),fig3)
% 
% HMsimPnL(1.1,1,2,0.01,0.01,100,0.01,5,exp(0.05),-0.5,0.5)
% figure(3);
% fig3=subplot(223);
% copyobj(allchild(get(figure(1),'CurrentAxes')),fig3)
% figure(4);
% fig3=subplot(223);
% copyobj(allchild(get(figure(2),'CurrentAxes')),fig3)
% 
% 
% HMsimPnL(1.15,1,2,0.01,0.01,100,0.01,5,exp(0.05),-0.5,0.5)
% figure(3);
% fig3=subplot(224);
% copyobj(allchild(get(figure(1),'CurrentAxes')),fig3)
% figure(4);
% fig3=subplot(224);
% copyobj(allchild(get(figure(2),'CurrentAxes')),fig3)
