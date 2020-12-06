function [Ks,dist,V,aprime] = find_Ks(kD0,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter)
% Given capital demand guess kD0, calculates Ks
na = length(assets);
pw = (1-palpha)*(kD0)^(palpha);
pr = palpha*(kD0)^(palpha - 1);
legal = true(na,na,2);
for nl = 1:2
   for nap = 1:na
       legal(nap,:,nl) = (pw*l(nl) + (pr+1-pdelta)*assets(nap) - assets')>0;
   end
end

V=0*V0;
aprime = V;
aind = V;
diff=999;
iter=1;
while ((iter<maxiter)&&(diff>tol))
    for nl = 1:2
    for ai = 1:na
        Val = (pw*l(nl) + (1+pr-pdelta)*assets(ai) - assets(legal(ai,:,nl))' ).^(1-pgamma)./(1-pgamma) + pbeta*(V0(legal(ai,:,nl),1)*Q(nl,1) + V0(legal(ai,:,nl),2)*Q(nl,2));
        if isempty(Val); Val = -1e9; end % infeasibility is rewarded with -infinity value - solves strange numerical problem that can arise far away from convergence
        [V(ai,nl),ind] = max(Val);
        aprime(ai,nl) = assets(ind);
        aind(ai,nl) = ind;
    end
    end
    iter = iter+1;
    diff = sum(abs(V(:)-V0(:)));
    V0 = V;
end

% stationary distribution
dist0 = ones(size(V));
dist0 = dist0./(sum(dist0(:)));
diff=999;
iter=1;
while ((iter<maxiter)&&(diff>tol))
    dist = 0*dist0;
    for nl = 1:2
       for ai = 1:na
           if dist0(ai,nl)>1e-15 % any weight? Otherwise skip to save computational time
           targ = aind(ai,nl);
           dist(targ,1) = dist(targ,1)+dist0(ai,nl)*Q(nl,1);
           dist(targ,2) = dist(targ,2)+dist0(ai,nl)*Q(nl,2);
           end
       end
    end
    iter = iter+1;
    diff = sum(abs((dist0(:)-dist(:))));
    dist0 = dist;
end
amarg = sum(dist,2);
if amarg(end)>1e-15; warning('Upper bound binding'); end
Ks = sum(amarg.*assets');
end