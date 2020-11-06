load matlab_x0=1point3.mat
clc;close all

% % marginal probility density funcion
% % optima q distrbution
% [TrainMetric,TestMetric]=generatePDF(qList.');
% title('q* probability density function','Interpreter','latex')
% xlabel('$q$','Interpreter','latex');
%
% % optima hat_theta_q distrbution
[TrainMetric,TestMetric]=generatePDF(thetaqList.')
title('$\hat{\theta}_q$ probability density function','Interpreter','latex')
xlabel('$\hat{\theta}_q$','Interpreter','latex');
%
% % joint distrbution
% % jointDistrubution(thetaqList.',qList.')







function [TrainMetric,TestMetric]=generatePDF(data)
%split the set as train set and test set

flag=randsample(length(data),fix(length(data)*0.8));
trainX=data(flag);
data(flag)=[];
testX=data;
% fit the  distribution
Wei=fitdist(trainX,'Weibull');
Gamma=fitdist(trainX,'Gamma');
Lognormal=fitdist(trainX,'Lognormal');
Normal=fitdist(trainX,'Normal');
% output the negativate log likelihood
TrainMetric=struct('Weibull',Wei.negloglik/length(trainX),'Gamma',Gamma.negloglik/length(trainX),...
    'Lognormal',Lognormal.negloglik/length(trainX),'Normal',Normal.negloglik/length(trainX));
% generate the pdf
x=linspace(0,max(trainX),length(trainX));

pdf_Wei = pdf(Wei,x);
pdf_Gamma = pdf(Gamma,x);
pdf_Lognormal = pdf(Lognormal,x);
pdf_Normal=pdf(Normal,x);


% y=cdf(Normal,data);
% y=[data,y];
% [~,Normal_p] = kstest(data,'CDF',y);

figure
histogram(trainX,20,'Normalization','pdf')
line(x,pdf_Wei,'LineStyle','-','Color','r')
line(x,pdf_Gamma,'LineStyle','-.','Color','b')
line(x,pdf_Lognormal,'LineStyle','--','Color','g')
line(x,pdf_Normal,'LineStyle','--','Color','black')


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