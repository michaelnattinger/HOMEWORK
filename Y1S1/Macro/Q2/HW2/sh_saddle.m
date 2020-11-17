function [kss,css,sad,diff] = sh_saddle(k,pbeta,pdelta,pz,pgamma,ppsi,dif)
% Calculates saddle path and steady state for growth model
kss = (((1/pbeta) - 1 + pdelta)/(ppsi*pz))^(1/(ppsi - 1));
css = pz*kss^(ppsi) - pdelta * kss;
[chk,chc] = LOM(kss,css,pz,ppsi,pdelta,pbeta,pgamma);
if abs(chk-kss)>1e-10 || abs(chc-css)>1e-10
    error('SS is not steady state')
end
sad = 0*k;
diff = sad;
for i=1:length(k)
cg0=dif:dif:2;
cg = cg0;
kg = 0*cg+k(i);
for tt=1:50
   [kg,cg] = LOM(kg,cg,pz,ppsi,pdelta,pbeta,pgamma);
end
[diff(i),ind] = min(abs(kg-kss)+abs(cg-css));
sad(i)=cg0(ind);
end
end

function [kk,cc] = LOM(k,c,pz,ppsi,pdelta,pbeta,pgamma)
kk = pz*k.^ppsi - c + (1-pdelta)*k;
cc = c.*(pbeta*(ppsi*pz*k.^(ppsi-1)+1-pdelta)).^(1/pgamma);
end