function [V,exRegion, M, exV, Vcheck] = HMdynamicFD(X,s0,gamma,kappa,c,F,R,G,b)%,T,n,m)
% Price an american option via crank nicolson finite differences
%
% Inputs:
%
%   X       initial belief hat{theta}_0 (starting point of brownian motion)
%   s0      Initial variance
%   gamma   risk aversion
%   kappa   wealth
%   c       marginal cost of information
%   F       Fixed cost of information
%   R       borrowing and lending payoff multiplier
%   G       expected excess risk free return
%   b       demand sensitivity to quality
%   T       time horizon
%   n+1     Number of points in time grid to use
%   m+2     Number of points in asset price grid to use

%set T = very close to s0
T = s0*9998/10000;
format long;
%set m and n, to use 1000 intervals
n=1000;
m=999;

if T > s0
    error('T>s0');
end

dt = T/(n); % Time step
% create time grid
t = 0:dt:T;
%remove last entry (to deal with limit problem.. set n>5000)
%not doing this anymore to keep n consistent.. else we have n-1 points
%instead.
%t(end) = [];

% Share price grid
Xmin = X-10*sqrt(s0);  % Minimum theta
%Xmin = X-12*sqrt(s0);  % Minimum theta
Xmax = X+30*sqrt(s0); % Maximum theta
dX = (Xmax-Xmin)/(m+1);
x = Xmin:dX:Xmax;

% Set up Grid
V = NaN*ones(n+1,m+2); % Pricing Matrix (t,x)
Vcheck = NaN*ones(n+1,m+2); % Pricing Matrix (t,x)


% Set Terminal condition on Grid
%V(end,:) = -10^20;
for i = 1:m+2
    %V(end,i) = g(x(i),s0,gamma,kappa,c,F,R,G,b,s0-0.0000001);
    %if -( max(G + b*x(i),0) )^2/(2*b^2) + gamma*c*(R) < 0
    %    V(end,i) = 0;
    %elseif -( max(G + b*x(i),0) )^2/(2*b^2) + gamma*c*(R) == 0
    %    V(end,i) = -( G+b*x(i) )^2 /(2*s0*b^2)-gamma*R*(kappa-F);
    %else
    %    V(end,i) = -10^100;
    %end
    V(end,i) = g(x(i),s0,gamma,kappa,c,F,R,G,b,t(end));
end
%V(end,:)

Vcheck(end,:) = V(end,:);

% Set Boundary condition on Grid, not including the last time period, since
% that was set earlier before in the above step.
for i = 1:n
    V(i,1) = g(Xmin,s0,gamma,kappa,c,F,R,G,b,t(i)); % Value of option when stock price is Xmin
    V(i,end) = g(Xmax,s0,gamma,kappa,c,F,R,G,b,t(i)); % Value of option when S = Smax
    %V(i,end) = 0;
end


Vcheck(:,1) = V(:,1);
Vcheck(:,end) = V(:,end);

% Create A matrix for finite difference method
%J = (1:m)';
cc = -1/(dX^2)*ones(m,1);
uu = 1/(2*(dX^2))*ones(m,1);
ll = 1/(2*(dX^2))*ones(m,1);
A = spdiags([ [ll(2:end);0] cc [0; uu(1:end-1)] ],[-1 0 1],m,m);

% Finite difference solver
%theta = 1/2; %crank-nicolson
theta = 1; %implicit

B = eye(m) + (1-theta)*dt*A;
B = sparse(B);
C = eye(m) - (theta)*dt*A ;
C = sparse(C);

% create matrix for exercise value
exV=zeros(n+1,m);
exV(end,:) = V(end,2:end-1);

%finite difference solver with boundary conditions added to matrix B and C
for nn = n+1:-1:2
    old = B*(V(nn,2:end-1)') + (1-theta)* [ll(1)*V(nn,1); zeros(m-2,1); uu(m)*V(nn,end)]+ (theta)*[ll(1)*V(nn-1,1); zeros(m-2,1); uu(m)*V(nn-1,end)];
    new = C\old; % Value of the option
    for i = 1:m
        exV(nn-1,i) = g(Xmin+i*dX,s0,gamma,kappa,c,F,R,G,b,t(nn-1));
    end
    V(nn-1,2:end-1) = max( new', exV(nn-1,:) );
    Vcheck(nn-1,2:end-1) = new';
    % stasify tau is less than s0-2*c*gamma*R
    % 4.656147122503568e+03
%     if t(nn-1)<(s0-2*c*gamma*R)
%         flag=x>(-G/b);
%         flag=flag(2:end-1);
%         V(nn-1,flag)=new(flag)';
%     end
end


% Check for exercise region
exRegion = (exV == V(:,2:end-1));
% fix the numerical error
flag=(x>(-G/b));
flag=flag(2:end-2);
exRegion(t<(s0-2*c*gamma*R),flag)=0;

%plot exercise region
figure (1)
%contourf(exRegion');
M = contour(exRegion','k');

xlabel('$q$','Interpreter','latex');
ylabel('$\hat{\theta}_q$','Interpreter','latex');

%yind = find(x(2:end-1)==X);
%yindm = find(x(2:end-1)==(Xmin+dX));
%yindM = find(x(2:end-1)==(Xmax-dX));

%yticks([yindm yind yindM]);
%yticklabels({Xmin+dX, X, Xmax-dX});

%xticks([1, n+1]);
%xticklabels({'0', '\infty'});

ax = gca;

yTicks = ax.YAxis.TickValues;
yTicks = Xmin+ yTicks*dX;

ax.XAxis.TickValues = ax.XAxis.TickValues+1; %+1 to get axis to plot 1001 point instead of 1000
xTicks = ax.XAxis.TickValues;
xTicks = round( (0+ (xTicks-1)*dt)./(s0^2-s0*(0+ (xTicks-1)*dt))); %-1 here because index now is from 1-1001, instead of 0 to 1000


yticklabels(yTicks);
xticklabels(xTicks);

end


%exercise value function
function value = g(X,s0,gamma,kappa,c,F,R,G,b,tau)

%value =  -exp( - ( (max(G + b*X,0))^2 )/(2*(s0-tau)*b^2) - gamma *( kappa- F - c*(tau)/(s0^2-s0*tau) )*(R)  );

if G+b*X > 0
    value = -exp( - ( (G + b*X)^2 )/(2*(s0-tau)*b^2) - gamma *( kappa- F - c*(tau)/(s0^2-s0*tau) )*(R)  );
else
    value = -exp(  - gamma *( kappa- F - c*(tau)/(s0^2-s0*tau) )*(R)  );
end


%if G+b*X > 0
%    value =  - ( (G + b*X)^2 )/(2*(s0-tau)*b^2) - gamma *( kappa- F - c*(tau)/(s0^2-s0*tau) )*(R)  ;
%else
%    value =  - gamma *( kappa- F - c*(tau)/(s0^2-s0*tau) )*(R)  ;
%end

%value = - value;

%value = - log(-value);
%if X >1.5
%    value = 0;
%end
end