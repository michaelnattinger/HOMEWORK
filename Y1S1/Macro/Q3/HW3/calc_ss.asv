function [Ybar,Cbar,Kbar,Lbar] = calc_ss(palpha, pbeta, pdelta, psigma, pphi, pAbar, pGbar, ptaubarI, ptaubarL)
syms Y K L C %pa palph pg pph psig ptL ptI pb pd
f1 = Y == pAbar*K^palpha * L^(1-palpha);
f2 = Y == C + pdelta*K + pGbar*Y;
f3 = L^pphi * C^psigma == (1-ptaubarL)*pAbar*(1-palpha)*K^(palpha)*L^(-palpha);
f4 = (1+ptaubarI) == pbeta *( pAbar * palpha * K^(palpha - 1) * L^(1-palpha) +(1-pdelta)*(1+ptaubarI));
soln = solve([f1 f2 f3 f4],[Y K L C]);
% check Y
validY = zeros(1,length(soln.Y));
validL = validY;
validC = validL;
for j=1:length(soln.Y)
    if (isreal(eval(soln.Y(j,1)))&& eval(soln.Y(j,1))>0)
        validY(j) = 1;
    end
    if (isreal(eval(soln.L(j,1)))&& eval(soln.L(j,1))>0)
        validL(j) = 1;
    end
    if (isreal(eval(soln.C(j,1)))&& eval(soln.C(j,1))>0)
        validC(j) = 1;
    end
end
if (sum(validY.*validL.*validC)==0); error('no valid solution'); end
%if (sum(validL)==0); error('no valid solution'); end
%if (sum(validC)==0); error('no valid solution'); end
%if (sum(valid)>1); warning('multiple valid solution'); end
ind = find((validY>0).*(validL>0) .*(validC>0));
Ybar = eval(soln.Y(ind(1),1));
Cbar = eval(soln.C(ind(1),1));
Kbar = eval(soln.K(ind(1),1));
Lbar = eval(soln.L(ind(1),1));
end