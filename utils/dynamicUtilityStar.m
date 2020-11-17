function f=dynamicUtilityStar(thetaq,q)
% calculate expected utility
gamma=0.01;
R=exp(0.05);
F=5;
kappa=100;
s0=0.01;
b=0.5;
G=-0.5;
c=0.01;
sq=s0./(1+q*s0);

f=-exp(-gamma*((G+b*thetaq).*omegaStar(thetaq,q,gamma,s0,b,G)+R*(kappa-F-c*q))+gamma^2*b^2*omegaStar(thetaq,q,gamma,s0,b,G).^2.*sq/2);
f=mean(f);
 
% paper lemma4
% if  not meet lemma, expected utility equals staticU(theta0,0.01)
% notice q ~=0

% if fullInformation(theta0)>=0
% f=staticUtilityStar(theta0,0.01);
% end
    
end

