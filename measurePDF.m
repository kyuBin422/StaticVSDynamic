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
% f=generalizedHyperbolicDistrbution(x,lambda,chi,psi,mu,sigma,gamma);
% [parmhat, se, parmci, output]=ghfit(thetaqList);
% GHFit(thetaqList)

x0=[10,10];

% options = optimoptions(@fmincon,'Algorithm','sqp-legacy','Display','iter','PlotFcns',@optimplotfval);

options = optimset('Display','iter','PlotFcns',@optimplotfval);
params=fminsearch(@(x) NormalLike(thetaqList,x),x0,options);

x=linspace(min(thetaqList),max(thetaqList),length(thetaqList));
y=PDF(x,params);
figure;
[counts,centers] = hist(thetaqList,100);
bar(centers,counts/trapz(centers,counts))
hold on
plot(x,y,'LineStyle','-','Color','r')
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

function validationMLE(data)
% divergence
flag=randsample(length(data),fix(length(data)*0.8));
trainX=data(flag);
data(flag)=[];
testX=data;

% empirical distribution function
meanEmpirical=mean(trainX);
varianceEmpirical=var(trainX);

% P true distrbution of data
% Q theory model

[Wei,Gamma,Lognormal,Normal]=generatePDF(trainX);

[counts,centers]=hist(testX,20);

PProbability=counts./length(testX);
histBinWidth=(max(testX)-min(testX))/20;

for QDistrbution=[Wei,Gamma,Lognormal,Normal]
    QProbability=cdf(QDistrbution,centers+histBinWidth/2)-cdf(QDistrbution,centers-histBinWidth/2);
    DKL_P_Q=sum(PProbability.*log(PProbability./QProbability));
    DKL_Q_P=sum(QProbability.*log(QProbability./PProbability));
    disp("divergence KL(P|Q) ="+ DKL_P_Q)
    disp("divergence KL(Q|K) ="+ DKL_Q_P)
    
end

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



function param=gradientDecent(data)
GHLikelihood=@ (x) abs(sum(log(generalizedHyperbolicDistrbution(data,x(1),x(2),x(3),x(4),x(5),x(6)))));
leaningRate=1e-3;
x0=[-0.5547,1.0373,0.9583,mean(data),std(data),skewness(data)];
param=struct('lambda',x0(1),'chi',x0(2),'psi',x0(3),'mu',x0(4),'sigma',x0(5),'gamma',x0(6));

while true
    StartValue=GHLikelihood(x0);
    [x1,result]=gradient(GHLikelihood,x0,'lambda');
    param.lambda=param.lambda-leaningRate*result;
    x0=x1;
    
    [x1,result]=gradient(GHLikelihood,x0,'chi');
    param.chi=param.chi-leaningRate*result;
    x0=x1;
    
    [x1,result]=gradient(GHLikelihood,x0,'psi');
    param.psi=param.psi-leaningRate*result;
    x0=x1;
    
    [x1,result]=gradient(GHLikelihood,x0,'mu');
    param.mu=param.mu-leaningRate*result;
    x0=x1;
    
    [x1,result]=gradient(GHLikelihood,x0,'sigma');
    param.sigma=param.sigma-leaningRate*result;
    x0=x1;
    
    [x1,result]=gradient(GHLikelihood,x0,'gamma');
    param.gamma=param.gammaleaningRate*result;
    x0=x1;
    EndValue=GHLikelihood(x0);
    param
    if (EndValue-StartValue)<1e-5
        break
    end
end
end

function [x1,result]=gradient(Likelihood,x0, flag)
result=0;
x1=x0;
switch flag
    case 'lambda'
        x1(1)=x0(1)+1e-4;
        result=(Likelihood(x1)-Likelihood(x0))/1e-4;
    case 'chi'
        x1(2)=x0(2)+1e-4;
        result=(Likelihood(x1)-Likelihood(x0))/1e-4;
    case 'psi'
        x1(3)=x0(3)+1e-4;
        result=(Likelihood(x1)-Likelihood(x0))/1e-4;
    case 'mu'
        x1(4)=x0(4)+1e-4;
        result=(Likelihood(x1)-Likelihood(x0))/1e-4;
    case  'sigma'
        x1(5)=x0(5)+1e-4;
        result=(Likelihood(x1)-Likelihood(x0))/1e-4;
    case  'gamma'
        x1(6)=x0(6)+1e-4;
        result=(Likelihood(x1)-Likelihood(x0))/1e-4;
end
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

function Logl=NormalLike(x, Params)
mu=Params(1);
sigma=Params(2);

LogL=-log(sqrt(2.*pi).*sigma)-((x-mu)./sigma).^2./2;

Logl=-sum(LogL)/length(x);
end

function f=PDF(x,Params)
mu=Params(1);
sigma=Params(2);

f=1./(sqrt(2.*pi).*sigma).*exp(-((x-mu)./sigma).^2./2);
end

