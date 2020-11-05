function f=omegaStar(thetaq,q)
gamma=0.01;
s0=0.01;
b=0.5;
G=-0.5;
sq=s0./(1+q*s0);
f =max((G+b*thetaq)./(gamma*b^2*sq),0);
end