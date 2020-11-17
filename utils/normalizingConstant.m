function a=normalizingConstant(chi,psi,lambda,sigma,gamma)
% ref https://projecteuclid.org/download/pdfview_1/euclid.jam/1425305853
% eqution 17
a=(chi.*psi).^(-lambda/2).*psi.^lambda.*(psi+(gamma.^2./sigma.^2)).^(1/2-lambda)./((2*pi)^0.5...
    .*sigma.*besselk(lambda,sqrt(psi.*chi)));
end