function [A,B] = LOM(rhoI0,rhoa,rhog,rhoL,palpha,pdelta,psigma,pphi,pGbar,pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar)
% calculates law of motion EX_p = AX +BZ for log linearized model
syms C K Cp Kp I Ip L Lp g a htauL htauI gp ap htauLp htauIp
% these equations form the law of motion of the endogenous variables
f1 = a + palpha*K + (1-palpha)*L == (Cbar/Ybar)*C + (pdelta*Kbar/Ybar)*I + pGbar*g; % two formulations for Y
f2 = pphi*L + psigma*C == (-1*ptaubarL/(1-ptaubarL))*htauL + palpha*K - palpha *L; % labor supply this period
f3 = Kp == (1-pdelta)*K + pdelta*I;
f4 = pphi*Lp + psigma*Cp == (-1*ptaubarL/(1-ptaubarL))*htauLp + palpha*Kp - palpha *Lp; % labor supply next period
f5 = psigma*(Cp - C) + (ptaubarI/(1+ptaubarI))*htauI == pbeta*((1+ptaubarI)^(-1)*palpha*pAbar*Kbar^(palpha - 1)*Lbar^(1-palpha)*(ap + (1-palpha)*(Lp - Kp)) +(1-pdelta)*(ptaubarI/(1+ptaubarI))*htauIp);
f6 = ap + palpha*Kp + (1-palpha)*Lp == (Cbar/Ybar)*Cp + (pdelta*Kbar/Ybar)*Ip + pGbar*gp; % two formulations for Y
% these equations describe the shock process
f7 = ap == rhoa*a;
f8 = gp == rhog*g;
f9 = htauLp == rhoL*htauL;
f10 = htauIp == rhoI0*htauI;
f = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10];
v = [Kp Cp L Lp I Ip ap gp htauLp htauIp];
obj = solve(f,v);
EXp = [obj.Kp; obj.Cp];
X = [K C];
Z = [a g htauL htauI];
A = eval(jacobian(EXp,X));
B = eval(jacobian(EXp,Z));
end