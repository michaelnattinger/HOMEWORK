function Kd = find_Kd(k0,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter)
% Sets up capital demand search as nonlinear optimization problem,
% and solves via fmincon.
% Nattinger, 2020
obj = @(x)abs(find_Ks(x,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter) -x);
Kd = fmincon(obj,k0,[-1;1],[-1e-14;assets(end)]);
end