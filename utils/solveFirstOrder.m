function solution=solveFirstOrder(theta0)
% the range of q
solution=0.01:0.01:6000;
% work out the value information and v'(q)
fVBar=vDerivative(theta0,solution);
fV=valueInformation(theta0,solution);

fV(fV<=0)=0;
% find the optimal q
[~,flagIndexfV]=max(fV);
[~,flagIndexfVBar]=min(abs(fVBar));
% if the value information is concave
% check v'(q*)==0 and v(q*) is maximal value
if abs(flagIndexfV-flagIndexfVBar)/flagIndexfVBar<0.01
    solution=solution(flagIndexfV);
else
    solution=solution(flagIndexfV);
end
end