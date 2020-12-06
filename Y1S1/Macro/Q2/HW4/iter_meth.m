function [kss,css,c,kprime,v] = iter_meth(k,pbeta,pdelta,pz0,pgamma0,ppsi)
% Iterative method of solving value function

kss = (((1/pbeta) - 1 + pdelta)/(ppsi*pz0))^(1/(ppsi - 1));
css = pz0*kss^(ppsi) - pdelta * kss;

nk=length(k);
v0 = (1/((1-pgamma0)*(1-pbeta)))*((1-pbeta)*k).^(1-pgamma0); % initial guess of value function approx geometric
diff = 1e6;
maxiter = 1e6;
tol=1e-6;
iter=1;
posval = ones(nk); %prealloc
for i=1:nk % defines what is possible
posval(i,k>pz0*k(i).^(ppsi) + (1-pdelta)*k(i)) = 0; %k' legal: c=0: max is pz*k^(ppsi) + (1-pdelta)*k
end
posval = logical(posval);
v=0*v0;
kprime = 0*v0;
% iterate until converged
while (diff>tol&&iter<maxiter)
    if mod(iter,100)==0; disp(['iteration: ' num2str(iter) ', diff = ' num2str(diff) ]); end
    for i=1:nk
        %k0=k(i); %If I start the period with k0 capital
        kpopt = k(posval(i,:)); % I have these options for k'
        vi = (pz0*k(i).^(ppsi) + (1-pdelta)*k(i) - kpopt).^(1-pgamma0)./(1-pgamma0) + pbeta*v0(posval(i,:));   %which yield me this much utility
        [v(i),ind] = max(vi); % value function is the best option
        kprime(i) = kpopt(ind);
    end
    diff = norm(v-v0);
    iter=iter+1;
    v0=v;
end
c = pz0*k.^(ppsi) + (1-pdelta)*k - kprime;
end

