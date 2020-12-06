function [Kd,diff,iter] = find_Kd_slowiter(tune,k0,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter,ktol,kmiter)
% Solves for equilibrium
% Guesses Kd, calculates Ks
% Updates guess to convex combination of Kd and Ks, with weight on Kd given
% by parameter tune, and initial Kd guess k0
% Note: fast but potentially unstable search with low tune, slow but more
% stable search with high tune
% Michael Nattinger, 2020
iter=1;
diff = 999;
Kd = k0;
while (iter<kmiter)&&(diff>ktol)
Ks= find_Ks(Kd,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter);
diff = abs(Kd-Ks);
Kd = tune*Kd+(1-tune)*Ks;
iter=iter+1;
end