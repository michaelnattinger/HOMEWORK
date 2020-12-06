function [Kd,diffm] = find_Kd_gridcon(assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter,ktol,kmiter)
center = 5;
ng = 5;
for i=0:5
    grid = center  + linspace(center-ng*(10^(-1*i)),center-ng*(10^(-1*i)),2*ng+1);
    Ks = 0*grid;
    for tt=1:length(Ks)
    Ks(i)= find_Ks(grid(i),assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter);
    end
    diff = abs(Ks - grid);
    [diffm,ind] = min(diff);
    center = grid(ind);
end
Kd = center;