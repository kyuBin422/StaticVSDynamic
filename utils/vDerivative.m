function f=vDerivative(theta0,q)
gamma=0.01;
R=exp(0.05);
F=5;
kappa=100;
s0=0.01;
b=0.5;
G=-0.5;
c=0.01;

f=  -exp(gamma*R*(c*q+F-kappa))./(2*(s0*q+1).^(3/2)).*(-2*c*gamma*R*(s0*q+1).^(3/2).*normcdf(-sqrt(s0*q+1).*(b*theta0+G)./(b*s0*sqrt(q)))-2*c*gamma*R*(s0*q+1).*exp(-(b*theta0+G)^2/(2*b^2*s0)).*normcdf((b*theta0+G)./(b*s0*sqrt(q)))+s0*exp(-(b*theta0+G)^2/(2*b^2*s0))*normcdf((b*theta0+G)./(b*s0*sqrt(q))));
end