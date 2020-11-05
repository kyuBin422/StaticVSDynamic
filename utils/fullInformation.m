function f= fullInformation(theta0)
gamma=0.01;
R=exp(0.05);
F=5;
s0=0.01;
b=0.5;
G=-0.5;

if theta0>-G/b
    f=F-1/(gamma*R)*(-log(normcdf(-(G+b*theta0)/(b*sqrt(s0))))-(b*theta0+G)^2/(2*b^2*s0));
elseif theta0==-G/b
    f=F-log(2)/(gamma*R);
elseif theta0<-G/b
    f=F-1/(gamma*R)*(-log(normcdf(-(G+b*theta0)/(b*sqrt(s0)))));
end

end