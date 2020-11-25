load matlab_x0=1point3.mat
clc;close all
qList=qList.';
thetaqList=thetaqList.';
addpath(genpath('utils'));
addpath(genpath('FlexDist'));
%%
% marginal probility density funcion
% optima q distrbution
[TrainMetric,TestMetric]=generatePDF(qList)
title('q* pdf','Interpreter','latex')
xlabel('$q$','Interpreter','latex');

% optima hat_theta_q distrbution
[TrainMetric,TestMetric]=generatePDF(thetaqList)
title('$\hat{\theta}_q$','Interpreter','latex')
xlabel('$\hat{\theta}_q$','Interpreter','latex');


%%
% joint distrbution
jointDistrubution(thetaqList.',qList.')
%%
% f=generalizedHyperbolicDistrbution(x,lambda,chi,psi,mu,sigma,gamma);
% [parmhat, se, parmci, output]=ghfit(thetaqList);
% GHFit(thetaqList)

x0=[-0.5547,1.0373,0.9583,0.5,0.5,0.5];

options = optimset('Display','iter','PlotFcns',@optimplotfval);
params=fminsearch(@(x) GHLike(thetaqList,x),x0,options);

y=generalizedHyperbolicDistrbution(0.5:0.001:2,params(1),params(2),params(3),params(4),params(5),params(6));
figure;
histogram(thetaqList,'Normalization','pdf')
line(0.5:0.001:2,y,'LineStyle','-','Color','r')
%%
param=gradientDecent(thetaqList);
y = generalizedHyperbolicDistrbution(sort(thetaqList),param.lambda,param.chi,param.psi,param.mu,param.sigma,param.gamma);
histogram(thetaqList,'Normalization','pdf')
line(sort(thetaqList),y,'LineStyle','-','Color','r')
%%
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
TrainMetric=struct('Weibull',Wei.negloglik/length(trainX),...
    'Gamma',Gamma.negloglik/length(trainX),...
    'Lognormal',Lognormal.negloglik/length(trainX),...
    'Normal',Normal.negloglik/length(trainX));
% generate the pdf
x=linspace(min(trainX),max(trainX),length(trainX));

pdf_Wei = pdf(Wei,x);
pdf_Gamma = pdf(Gamma,x);
pdf_Lognormal = pdf(Lognormal,x);
pdf_Normal=pdf(Normal,x);


% y=cdf(Normal,data);
% y=[data,y];
% [~,Normal_p] = kstest(data,'CDF',y);

figure
histogram(trainX,'Normalization','pdf')
line(x,pdf_Wei,'LineStyle','-','Color','r')
line(x,pdf_Gamma,'LineStyle','-.','Color','b')
line(x,pdf_Lognormal,'LineStyle','--','Color','g')
line(x,pdf_Normal,'LineStyle','-','Color','black')


legend('Data','Weibull','Gamma','Lognormal','Normal','Location','northwest')


ylabel('probability')
% generate the test pdf

pdf_Wei = pdf(Wei,testX);
pdf_Gamma = pdf(Gamma,testX);
pdf_Lognormal = pdf(Lognormal,testX);
pdf_Normal=pdf(Normal,testX);
TestMetric=struct(...
    'Weibull',-sum(log(pdf_Wei))/length(testX),...
    'Gamma',-sum(log(pdf_Gamma))/length(testX),...
    'Lognormal',-sum(log(pdf_Lognormal))/length(testX),...
    'Normal',-sum(log(pdf_Normal))/length(testX)...
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
xlabel('$q$','Interpreter','latex');
ylabel('$\hat{\theta}_q$','Interpreter','latex');
figure;
histogram2(q, theta,20,'Normalization','pdf','FaceColor','flat','ShowEmptyBins','on')
colorbar
xlabel('$q$','Interpreter','latex');
ylabel('$\hat{\theta}_q$','Interpreter','latex');

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

function phat=GHFit(x)


% generalizedHyperbolicDistrbution(x,lambda,chi,psi,mu,sigma,gamma)

% x0=[-0.5547,1.0373,0.9583,mean(data),std(data),skewness(data)];
%
% ghlike=@(x) sum(log(generalizedHyperbolicDistrbution(data,x(1),x(2),x(3),x(4),x(5),x(6))));
% parameter = fmincon(ghlike, x0);
%
% x=linspace(min(data),max(data),length(data));
% PDF=generalizedHyperbolicDistrbution(x,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6));
% line(x,PDF,'LineStyle','-','Color','r')


PDF=@ (x,lambda,chi,psi,mu,sigma,gamma) sum(log(generalizedHyperbolicDistrbution(x,lambda,chi,psi,mu,sigma,gamma)));

x0=[-0.5547,1.0373,0.9583,mean(x),std(x),skewness(x)];

LB=[-1,0,0,0,0,-5];
UB=[1,10,10,10,100,5];
phat = mle(x,'pdf',PDF,'start',x0,'LowerBound',LB,'UpperBound',UB);
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
a= -lambda/2.*log(chi.*psi)+lambda.*log(psi)+(1/2-lambda).*log(psi+(gamma.^2./sigma.^2))-...
    1/2*log(2*pi)-log(sigma)-log(besselk(lambda,sqrt(psi.*chi)));

% compute the generalized hyperbolic distrbution
LogL=a+log(besselk(lambda-0.5,sqrt((chi+(x-mu).^2./sigma^.2)...
    .*(psi+gamma.^2./sigma.^2))))+gamma.*(x-mu)./sigma.^2 ...
    +(lambda./2-1/4)*(log(chi+(x-mu).^2./sigma.^2)+log(psi+gamma.^2./sigma.^2));

LogL=-sum(LogL)/length(x);

LogL(~isfinite(LogL))  =  1.0e-20;
LogL(~(~imag(LogL)))   =  1.0e-20;

end

