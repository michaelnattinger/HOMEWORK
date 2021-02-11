function [rhoI,htauI] = BKFP(rhoI0,tol,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
% Solves for fixed point solution to htauT by continually applying
% Blanchard-Kahn.
diff = 999;
maxiter = 1e6;
iter = 1;
while (iter<maxiter)&&(diff>tol)
    [rhoI,htauI] = BK_calcrho(rhoI0,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,psigma,pphi ...
        ,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    diff = abs(rhoI0 - rhoI);
    iter = iter + 1;
    rhoI0 = rhoI;
end
end

function [rhoI,htauI] = BK_calcrho(rhoT0,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,psigma,pphi ,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
    htauI = BK_calchtauI(rhoT0,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,psigma,pphi ,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
    rhoI = htauI(1:end-1)\htauI(2:end);
end

function htauI = BK_calchtauI(rhoI,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,psigma,pphi ,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
[A,B] = LOM(rhoI,rhoa,rhog,rhoL,palpha,pdelta,psigma,pphi  ...
    ,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
[Q,Lambda] = eig(A);
iQ = inv(Q);
if abs(Lambda(1))>1;sel = 1; else; sel = 2; end % which row has explosive eig?
C = iQ*B;
lm = diag(Lambda);
lam = lm(sel);
Theta = (-1/lam)*C(sel,:)*inv(eye(4) -(1/lam)*diag([rhoa rhog rhoL rhoI]) );
zs = [a g htauL];
vs = [k c ];
htauI = Theta(4)^(-1) * ( vs* iQ(sel,:)' - (zs*Theta(1,1:3)')); % follows formula from writeup
end