function [c,aprime,V] = iter_meth_sto(k,pbeta,pdelta,z,P,pgamma0,ppsi)
nk=length(k);
v0 = [log(k') log(k')]; % initial guess of value function approx geometric
diff = 1e6;
maxiter = 1e6;
tol=1e-6;
iter=1;
posval = true(nk,nk,2);
aprime = 0*v0;
aind = aprime;
for zi = 1:2
for i=1:nk % defines what is possible
posval(i,k>z(zi)*k(i).^(ppsi) + (1-pdelta)*k(i),zi) = 0; %k' legal: c=0: max is pz*k^(ppsi) + (1-pdelta)*k
end
end
while (iter<maxiter)&&(diff>tol)
    V = v0;
for zi = 1:2
   for ki = 1:nk
       val = (z(zi)*k(ki).^(ppsi) + (1-pdelta).*k(ki) - k(posval(ki,:,zi)')).^(1-pgamma0)./(1-pgamma0) + pbeta*(P(zi,1)*v0(posval(ki,:,zi),zi)' + P(zi,2)*v0(posval(ki,:,zi),zi)');
       [V(ki,zi),ind] = max(val);
       aprime(ki,zi) = k(ind);
       aind(ki,zi) = ind;
   end
end
diff = sum(abs(V(:) - v0(:)));
iter = iter+1;
v0 = V;
c = repmat(z,nk,1).*repmat(k',1,2).^ppsi +(1-pdelta).*repmat(k',1,2) - aprime;
end

