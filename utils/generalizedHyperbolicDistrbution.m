function f=generalizedHyperbolicDistrbution(x,lambda,chi,psi,mu,sigma,gamma)
% ref https://projecteuclid.org/download/pdfview_1/euclid.jam/1425305853
% compute the normalizing constant
a=normalizingConstant(chi,psi,lambda,sigma,gamma);
% compute the generalized hyperbolic distrbution
f=a.*besselk(lambda-0.5,sqrt((chi+(x-mu).^2./sigma^.2)...
    .*(psi+gamma.^2./sigma.^2))).*exp(gamma.*(x-mu)./sigma.^2)...
    .*((chi+(x-mu).^2./sigma.^2).*(psi+gamma.^2./sigma.^2)).^(lambda./2-1/4);
end