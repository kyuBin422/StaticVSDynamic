clc;clear;close all
load matlab_x0=1point3.mat
addpath(genpath('utils'));
addpath(genpath('FlexDist'));
%%
% marginal probility density funcion
% optima q distrbution
[TrainMetric,TestMetric]=generatePDF(qList);
title('q* pdf','Interpreter','latex','FontSize',12,'FontWeight','bold')
xlabel('$q$','Interpreter','latex','FontSize',12,'FontWeight','bold');

% optima hat_theta_q distrbution
[TrainMetric1,TestMetric1]=generatePDF(thetaqList);
title('$\hat{\theta}_q pdf $','Interpreter','latex','FontSize',12,'FontWeight','bold')
xlabel('$\hat{\theta}_q$','Interpreter','latex','FontSize',12,'FontWeight','bold');

figure;
subplot(211)
plot([struct2array(TrainMetric),-8.2352],'--')
hold on
plot([struct2array(TestMetric),-8.2534],'-.')
xticks([1:5])
xticklabels(["Weibull","Gamma","Lognormal","Normal","Generalized Hyperbolic"])
ylabel('log likelihood','FontSize',12,'FontWeight','bold')
title('q* log likelihood  value')
legend('TrainSet','TestSet')

subplot(212)
hold on 
plot([struct2array(TrainMetric1),0.9509],'--')
hold on
plot([struct2array(TestMetric1),0.9604],'-.')
xticks([1:5])
xticklabels(["Weibull","Gamma","Lognormal","Normal","Generalized Hyperbolic"])
ylabel('log likelihood','FontSize',12,'FontWeight','bold')
title('thetaq log likelihood  value ','FontSize',12,'FontWeight','bold')
legend('TrainSet','TestSet')
%%
% joint distrbution
jointDistrubution(thetaqList.',qList.')
%%
% x0=[0.05,1.374,0.9583,1,1,0.3451];
x0=[0.05,1.374,0.9583,1,1,0.3451];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[-Inf,-Inf,-Inf,0,0,-Inf];
up=[100,100,100,100,100,100];
nonlcon=[];

options = optimoptions(@fmincon,'Display','iter','PlotFcns',@optimplotfval);
params=fmincon(@(x) GHLike(thetaqList,x),x0,A,b,Aeq,beq,lb,up,nonlcon,options);



x=linspace(min(thetaqList),max(thetaqList),length(thetaqList));
y=generalizedHyperbolicDistrbution(x,params(1),params(2),params(3),params(4),params(5),params(6));
figure;
[counts,centers] = hist(thetaqList,100);
bar(centers,counts/trapz(centers,counts))
hold on
plot(x,y,'LineStyle','-','Color','r')
xlabel('$\hat{\theta}_q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
ylabel("density",'FontSize',12,'FontWeight','bold')

GHLike(thetaqList,params)


GHLike(thetaqList,params)
function [TrainMetric,TestMetric]=generatePDF(data)
%split the set as train set and test set

flag=randsample(length(data),fix(length(data)*0.8));
trainX=data(flag);
data(flag)=[];
testX=data;
clear flag
% fit the  distribution
Wei=fitdist(trainX,'Weibull');
Gamma=fitdist(trainX,'Gamma');
Lognormal=fitdist(trainX,'Lognormal');
Normal=fitdist(trainX,'Normal');
% output the negativate log likelihood
TrainMetric=struct('Weibull',-Wei.negloglik/length(trainX),...
    'Gamma',-Gamma.negloglik/length(trainX),...
    'Lognormal',-Lognormal.negloglik/length(trainX),...
    'Normal',-Normal.negloglik/length(trainX));
% generate the pdf
x=linspace(min(trainX),max(trainX),length(trainX));

pdf_Wei = pdf(Wei,x);
pdf_Gamma = pdf(Gamma,x);
pdf_Lognormal = pdf(Lognormal,x);
pdf_Normal=pdf(Normal,x);


% y=cdf(Normal,data);
% y=[data,y];
% [~,Normal_p] = kstest(data,'CDF',y);

figure;
histogram(trainX,'Normalization','pdf')
line(x,pdf_Wei,'LineStyle','-','Color','r')
line(x,pdf_Gamma,'LineStyle','-.','Color','b')
line(x,pdf_Lognormal,'LineStyle','--','Color','g')
line(x,pdf_Normal,'LineStyle','-','Color','black')


legend('Data','Weibull','Gamma','Lognormal','Normal')


ylabel('density','FontSize',12,'FontWeight','bold')
% generate the test pdf

pdf_Wei = pdf(Wei,testX);
pdf_Wei(pdf_Wei<=0)=0.0001;
pdf_Gamma = pdf(Gamma,testX);
pdf_Gamma(pdf_Gamma<=0)=0.0001;
pdf_Lognormal = pdf(Lognormal,testX);
pdf_Lognormal(pdf_Lognormal<=0)=0.0001;
pdf_Normal=pdf(Normal,testX);
pdf_Normal(pdf_Normal<=0)=0.0001;

TestMetric=struct(...
    'Weibull',sum(log(pdf_Wei))/length(testX),...
    'Gamma',sum(log(pdf_Gamma))/length(testX),...
    'Lognormal',sum(log(pdf_Lognormal))/length(testX),...
    'Normal',sum(log(pdf_Normal))/length(testX)...
    );
end

function jointDistrubution(theta,q)
figure;
scatterhist(q,theta,'Kernel','overlay')
xlabel('$q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
ylabel('$\hat{\theta}_q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
figure;
histogram2(q, theta,20,'Normalization','pdf','FaceColor','flat','ShowEmptyBins','on')
colorbar
xlabel('$q$','Interpreter','latex','FontSize',12,'FontWeight','bold');
ylabel('$\hat{\theta}_q$','Interpreter','latex','FontSize',12,'FontWeight','bold');

%
%
% x=theta;
% y=q;
% figure;
% scatterhist(x,y)
%
% u = ksdensity(x,x,'function','cdf');
% v = ksdensity(y,y,'function','cdf');
%
% figure;
% scatterhist(u,v)
% xlabel('u')
% ylabel('v')
%
% rng default  % For reproducibility
% [Rho,nu] = copulafit('t',[u v],'Method','ApproximateML')
% r = copularnd('t',Rho,nu,1000);
% u1 = r(:,1);
% v1 = r(:,2);
%
% figure;
% scatterhist(u1,v1)
% xlabel('u')
% ylabel('v')
% set(get(gca,'children'),'marker','.')
% x1 = ksdensity(x,u1,'function','icdf');
% y1 = ksdensity(y,v1,'function','icdf');
%
% figure;
% scatterhist(x1,y1)
% set(get(gca,'children'),'marker','.')
end

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

% LogL(~isfinite(LogL))  =  1.0e-20;
% LogL(~(~imag(LogL)))   =  1.0e-20;

end


