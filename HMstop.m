function stopping_indx = HMstop(M, bmT, bmS)

% M is the contour matrix
% bmT is the time of simulated brownian motion
% bmS is the value of simulated brownian motion


%take relavent boundaries and transpose
Mex = M(:,2:( M(2,1) ) )';

%form matrix for bm path
% path = [bmT, bmS];


% Mex is N*2 matrix
% Path is N*2 matrix

% l =calculate(Mex,path);
if bmS(1)>max(Mex(:,2)) || bmS(1)<min(Mex(:,2))
    stopping_indx=1;
else
    l =calculateStoppingIndex(Mex(:,1),Mex(:,2),bmT,bmS);
    if isempty(l) || l<=0
        stopping_indx=1;
    else
        stopping_indx=fix(l(1));
    end
end

end

% function l=calculate(Mex,path)
% l=[];
% % point in mex
% x3=Mex(1:end-1,1);
% x4=Mex(2:end,1);
% y3=Mex(1:end-1,2);
% y4=Mex(2:end,2);
% for i=1:size(path,1)-1
%     % point in path
%     x1=path(i,1);
%     x2=path(i+1,1);
%     y1=path(i,2);
%     y2=path(i+1,2);
%
%     % calculate the parameter of P D Q
%     P=(x4-x2).*(y4-y3)-(x4-x3).*(y4-y2);
%     D=(x1-x2).*(y4-y3)-(x4-x3).*(y1-y2);
%     Q=(x1-x2).*(y4-y2)-(x4-x2).*(y1-y2);
%
%     if D~=0
%         lam=P./D;
%         eta=Q./D;
%         if ~(lam>=0&lam<=1&eta>=0&eta<=1)
%             continue
%         end
%         lam=lam(lam>=0&lam<=1&eta>=0&eta<=1&D~=0);
%         l=[lam * x1 + (1 - lam) * x2, lam * y1 + (1 - lam) * y2];
%         break
%     end
%     if P~=0 | Q~=0
%         continue
%     end
%     t1=sort(x1,x2);
%     t2=[min(x3,x4),max(x3,x4)];
%     t2=t2(D==0 & P==0 & Q==0);
%
%     if  t1(2)<t2(:,1) | t2(:,2) <t1(1)
%         continue
%     end
%     if t1(1)>t2(:,1)
%         l=t1(1);
%     else
%         l=t2(t1(1)<=t2(:,1),1);
%     end
%     break
% end
% end




