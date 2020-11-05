function f=valueInformation(theta0,q)
gamma=0.01;
R=exp(0.05);
F=5;
kappa=100;
s0=0.01;
b=0.5;
G=-0.5;
c=0.01;
ex=staticUtilityStar(theta0,q);
if theta0>-G/b
    f=ex+exp(-gamma*kappa*R-(b*theta0+G)^2/(2*b^2*s0));
elseif theta0==-G/b
    f=exp(-gamma*kappa*R)-exp(-gamma*R*(kappa-c*q-F)).*(1+1./(sqrt(1+s0*q)))/2;
elseif theta0<-G/b
    f=ex+exp(-gamma*kappa*R);
end

end