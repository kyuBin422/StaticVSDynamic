function stopping_indx = HMstop(M, bmT, bmS)

% M is the contour matrix
% bmT is the time of simulated brownian motion
% bmS is the value of simulated brownian motion


%take relavent boundaries and transpose
Mex = M(:,2:( M(2,1) ) )';

%form matrix for bm path
path = [bmT, bmS];



% [~, locb] = ismembertol(path,Mex, 3, 'ByRows', true, 'DataScale',1);


% stopping_indx = find(locb,1);
l =calculate(Mex,path);
if isempty(l)
    stopping_indx=1;
else
    stopping_indx=fix(l(1));
end

end

function l=calculate(Mex,path)
l=[];

for i=1:size(path,1)-1
    % point in path
    x1=path(i,1);
    x2=path(i+1,1);
    y1=path(i,2);
    y2=path(i+1,2);
    % point in mex
    x3=Mex(1:end-1,1);
    x4=Mex(2:end,1);
    y3=Mex(1:end-1,2);
    y4=Mex(2:end,2);
    % calculate the parameter of P D Q
    P=(x4-x2).*(y4-y3)-(x4-x3).*(y4-y2);
    D=(x1-x2).*(y4-y3)-(x4-x3).*(y1-y2);
    Q=(x1-x2).*(y4-y2)-(x4-x2).*(y1-y2);
    
    if D~=0
        flag1=D~=0;
        lam=P./D;
        eta=Q./D;
        lam=lam(flag1);
        eta=eta(flag1);
        
        if ~(lam>=0&lam<=1&eta>=0&eta<=1)
            continue
        end
        flag2=lam>=0&lam<=1&eta>=0&eta<=1;
        lam=lam(flag2);
        eta=eta(flag2);
        l=[lam * x1 + (1 - lam) * x2, lam * y1 + (1 - lam) * y2];
        break
    end
    if P~=0 | Q~=0
        continue
    end
    flag3=P~=0 | Q~=0;
    t1=sort(x1,x2);
    t2=[min(x3,x4),max(x3,x4)];
    t2=t2(~flag1 &~flag3);
    
    if  t1(2)<t2(:,1) | t2(:,2) <t1(1)
        continue
    end
    if t1(1)>t2(:,1)
        l=t1(1);
    else
        l=t2(t1(1)<=t2(:,1),1);
    end
    break
end
end

function l=calculateStoppingIndex(Mex,path)
l=[];
for i=1:size(path,1)-1
    for j =1:size(Mex,1)-1
        % point in path
        x1=path(i,1);
        x2=path(i+1,1);
        y1=path(i,2);
        y2=path(i+1,2);
        % point in mex
        x3=Mex(j,1);
        x4=Mex(j+1,1);
        y3=Mex(j,2);
        y4=Mex(j+1,2);
        % calculate the parameter of P D Q
        P=det([x4-x2,x4-x3;y4-y2,y4-y3]);
        D=det([x1-x2,x4-x3;y1-y2,y4-y3]);
        Q=det([x1-x2,x4-x2;y1-y2,y4-y2]);
        
        if D~=0
            lam=P/D;
            eta=Q/D;
            if ~(lam>=0&&lam<=1&&eta>=0&&eta<=1)
                continue
            end
            l(end+1,:)=[lam * x1 + (1 - lam) * x2, lam * y1 + (1 - lam) * y2];
            break
            
        end
        if P~=0 || Q~=0
            continue
        end
        t1=sortrows([x1,y1;x2,y2]);
        t2=sortrows([x3,y3;x4,y4]);
        if t1(2,1)<t2(1,1) || t2(2,1) <t1(1,1)
            continue
        end
        if t1(1,1)>t2(1,1)
            l(end+1,:)=t1(1,:);
        else
            l(end+1,:)=t2(1,:);
        end
        break
    end
    if size(l,1)>=1
        break
    end
end
end





