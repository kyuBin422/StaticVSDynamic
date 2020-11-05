function [S1plot, T1plot, S2plot, T2plot] = HMplotBM

% create brownian motion with drift = 0, volatility 1
obj = bm(0,1);

dt = 0.01/1000;

% simulate a sample path for 1000 periods, with dt 
[S1, T1] = simulate(obj, 1000, 'DeltaTime', dt);

[S2, T2] = simulate(obj, 1000, 'DeltaTime', dt);


% transform paths to plot on other figure axis.. 
X = 1;
s0 = 0.01;
m=999;

Xmin = X-3*sqrt(s0);  % Minimum theta
Xmax = X+12*sqrt(s0); % Maximum theta
dX = (Xmax-Xmin)/(m+1);

S1 = S1+X;
S2 = S2+X;

S1plot = (S1 - Xmin)/dX;
T1plot = T1/dt+1;

S2plot = (S2 - Xmin)/dX;
T2plot = T2/dt+1;

plot(T1plot, S1plot);
hold on
plot(T2plot, S2plot);


% codes for plotting the threshold
%investThreshold = ones(1,length(T1plot));
% Xmin = 1 - 3*sqrt(0.01)
% Xmax = 1+12*sqrt(0.01)
% dX = (Xmax-Xmin)/(1000)
% investThreshold = investThreshold *(1-Xmin)/dX
% plot(T1plot,investThreshold, '--');
% [V, exRegiongamma001]=HMdynamicFD(1, 0.01,0.01,1000,0.1,200,exp(0.05),-0.5,0.5,0.01,1000,999);
% hold on
% plot(T1plot(1:950), S1plot(1:950), ':k')
% plot(T1plot(1:650), S2plot(1:650), ':')
% plot(T1plot,investThreshold, '--k');
