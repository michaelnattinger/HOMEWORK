function [Ktraj,Ctraj,diff] = eq_traj(kss,css,D,palpha,psigma,pbeta,pdelta)
C0 = 3.82:0.0000001:css; % see also comment about laziness in calc_saddle.m
C = C0;                  % but also there are some numerical stability
K = kss+0*C;             % issues with really high horizons so recalculating
nd = length(D);          % a converged starting point can result in
C = repmat(C,nd,1);      % nonconvergence, which this avoids
K = repmat(K,nd,1);
for i=2:nd
    [K(i,:),C(i,:)] = dprop(K(i-1,:),C(i-1,:),D(i),palpha,psigma,pbeta,pdelta);
end
[diff,ind] = min(abs(C(end,:)-css)+abs(K(end,:)-kss));
Ctraj = C(:,ind)';
Ktraj = K(:,ind)';
end