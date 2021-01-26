function [Ktraj,Ctraj,diff] = eq_traj(kss,css,D,palpha,psigma,pbeta,pdelta)
C0 = 3.75:0.0000001:css;
C = C0;
K = kss+0*C;
nd = length(D);
C = repmat(C,nd,1);
K = repmat(K,nd,1);
for i=2:nd
    [K(i,:),C(i,:)] = dprop(K(i-1,:),C(i-1,:),D(i),palpha,psigma,pbeta,pdelta);
end
[diff,ind] = min(abs(C(end,:)-css)+abs(K(end,:)-kss));
Ctraj = C(:,ind)';
Ktraj = K(:,ind)';
% Ctraj = 0*D;
% Ktraj = 0*D;
% Ktraj(1) = kss;
% Ctraj(1) = C0(ind);
% for i=2:nd
%     [Ktraj(i),Ctraj(i)] = dprop(Ktraj(i-1),Ctraj(i-1),D(i),palpha,psigma,pbeta,pdelta);
% end