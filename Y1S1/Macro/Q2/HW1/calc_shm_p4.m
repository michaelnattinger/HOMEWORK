function [cshm,ktrj,ctrj] = calc_shm_p4(k,ctarg,ktarg,l,u,n,delt,bet,gamm,eta,g)
cg = linspace(l,u,n);
k0 = k+cg*0;
c0 = cg;
for tt=1:50
k1 = (1-delt)*k0 + log(k0) - c0;
c1 = c0.*(g^(eta*(1-gamm))*bet*(1-delt + (1./k1))).^(1/gamm);
k0 = k1;
c0 = c1;
end
[~,ind] = min(abs(c1 - ctarg)+abs(k1 - ktarg));
cshm = cg(ind);

ktrj = zeros(1,50);
ctrj = ktrj;
ktrj(1) = k;
ctrj(1) = cshm;
for tt=2:50
    ktrj(tt) = (1-delt)*ktrj(tt-1) + log(ktrj(tt-1)) - ctrj(tt-1);
    ctrj(tt) = ctrj(tt-1).*(g^(eta*(1-gamm))*bet*(1-delt + (1./ktrj(tt)))).^(1/gamm);
end
