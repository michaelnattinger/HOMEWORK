function [Kd,diff,iter] = find_Kd_slowiter(tune,k0,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter,ktol,kmiter)
iter=1;
diff = 999;
Kd = k0;
while (iter<kmiter)&&(diff>ktol)
Ks= find_Ks(Kd,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter);
diff = abs(Kd-Ks);
Kd = tune*Kd+(1-tune)*Ks;
iter=iter+1;
end