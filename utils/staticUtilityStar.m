function f=staticUtilityStar(theta0,q,s0,gamma,kappa,c,F,R,G,b)

% gamma=0.01;
% R=exp(0.05);
% F=5;
% kappa=100;
% s0=0.01;
% b=0.5;
% G=-0.5;
% c=0.01;

d= @(theta0,q) -(G+b*theta0)*sqrt(1+s0*q)./(b*s0*sqrt(q));

f=-exp(-gamma*R*(kappa-c*q-F)).*(normcdf(d(theta0,q))+1./sqrt(1+s0*q).*exp(-(G+b*theta0)^2/(2*b^2*s0)).*(1-normcdf(-(G+b*theta0)./(b*s0*sqrt(q)))));
end

