function [kss,css,sad,diff] = nlopt_saddle(k,pbeta,pdelta,pz,pgamma,ppsi)
% Calculates saddle path and steady state for growth model

kss = fpinv((1/pz)*((1/pbeta) - 1 + pdelta),ppsi);
css = pz*kss^(ppsi) - pdelta * kss;
[chk,chc] = LOM(kss,css,pz,ppsi,pdelta,pbeta,pgamma);
if abs(chk-kss)>1e-10 || abs(chc-css)>1e-10
    error('SS is not steady state')
end

sad = 0*k;
diff = sad;
for i=1:length(k)
    of = @(x) objfun(x,k(i),[kss css],50,pz,ppsi,pdelta,pbeta,pgamma);
    [sad(i),diff(i)] = fminunc(of,max([(((pz*k(i)^(ppsi) - pdelta*k(i))+(pz*k(i)^ppsi + (1-pdelta)*k(i) - (1/pbeta)+ 1 - pdelta))/2) 0.1]));
end
end

function absdif = objfun(x0,ki,xf,niter,pz,ppsi,pdelta,pbeta,pgamma)
x = [ki x0];
    for i=1:niter; [x(1), x(2)] = LOM(x(1),x(2),pz,ppsi,pdelta,pbeta,pgamma); end
    absdif = sum(abs(x-xf));
    if isnan(absdif); absdif=1e10+rand(); end
    if absdif>1e10; absdif=1e10+rand(); end
end

function y = fpinv(x,ppsi)
y = (x/ppsi)^(1/(1-ppsi));
end

function [kk,cc] = LOM(k,c,pz,ppsi,pdelta,pbeta,pgamma)
kk = pz*k.^ppsi - c + (1-pdelta)*k;
cc = c*(pbeta*(ppsi*pz*k.^(1-ppsi)+1-pdelta))^(1/pgamma);
end