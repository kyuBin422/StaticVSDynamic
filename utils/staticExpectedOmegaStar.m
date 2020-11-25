function f=staticExpectedOmegaStar(theta0,q,gamma,s0,b,G)
f=(G+b*theta0)*(1+s0*q)/(s0*b^2*gamma)*normcdf((G+b*theta0)*sqrt(1+s0*q)/(b*s0*sqrt(q)))+...
    sqrt(q+s0*q^2)/(b*gamma)*normpdf(-(G+b*theta0)*sqrt(1+s0*q)/(b*s0*sqrt(q)));
end
