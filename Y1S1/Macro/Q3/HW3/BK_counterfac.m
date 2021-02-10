function [ga,gg,ghtauL,ghtauI] = BK_counterfac(rhoI,rhoa,rhog,rhoL,a,g,htauL,htauI,c,k,palpha,pdelta,psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
v = [k c];
z = [a g htauL htauI];
[A,B] = LOM(rhoI,rhoa,rhog,rhoL,palpha,pdelta,psigma,pphi  ...
    ,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
[Q,Lambda] = eig(A);
iQ = inv(Q);
if abs(Lambda(1))>1;sel = 1; else; sel = 2; end % which row has explosive eig?
C = iQ*B;
lm = diag(Lambda);
lam = lm(sel);
Theta = (-1/lam)*C(sel,:)*inv(eye(4) -(1/lam)*diag([rhoa rhog rhoL rhoI]) );

for i=1:4
    switch i
        case 1
            z = [a 0*a 0*a 0*a]; 
        case 2
            z = [0*a g 0*a 0*a];
        case 3
            z = [0*a 0*a htauL 0*a];
        case 4
            z = [0*a 0*a 0*a htauI];
    end
    % generate c,k time series counterfactual % I am doing this wrong
    vi = v;
    l = 0*v(:,1);
    y = l;
    
    for j=1:length(vi(:,1)) % k starts at its initial point
        vi(j,2) = (1/iQ(sel,2))*( -1*iQ(sel,1)*vi(j,1) + Theta*z(j,:)');
        l(j) = ((-psigma)/(pphi + palpha)) * vi(j,2) + ((palpha)/(pphi + palpha))* vi(j,1) + ((-1*ptaubarL)/((palpha + pphi)*(1-ptaubarL)))*z(j,3);
        if j<length(vi(:,1))
            vi(j+1,1) = (1-pdelta)*vi(j,1) + pdelta*(Ybar/(pdelta*Kbar))*(z(j,1) + palpha*vi(j,1) + (1-palpha)*l(j) - (Cbar/Ybar)*vi(j,2) - pGbar*z(j,2));
        end
    end
    % from c,k,tau_L calculate l
    %l = ((-psigma)/(pphi + palpha)) * vi(:,2) + ((palpha)/(pphi + palpha))* vi(:,1) + ((ptaubarL)/((palpha + pphi)*(1-ptaubarL)))*z(:,3);
    % from a,k,l calculate y
    y = z(:,1) + palpha*vi(:,1) +(1-palpha)*l;
    switch i
        case 1
            ga = y;
        case 2
            gg = y;
        case 3
            ghtauL = y;
        case 4
            ghtauI = y;
    end
end
end