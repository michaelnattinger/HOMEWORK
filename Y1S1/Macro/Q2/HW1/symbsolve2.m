clear; close all; clc

% syms an bn gn a b c del
% andrew soln
% f1 = an == a + (c*an*del)/(c+bn*del);
% f2 = bn == b + (c*bn*del)/(c+bn*del);
% f3 = gn == del*gn + (1/2)*(an^2*del^2)/(c+bn*del);
% 
% solns = solve([f1 f2 f3],[an bn gn]);
% 
% f4 = an == a+(c*del^2*an*bn)/((c+del*bn)^2) + (del*an*c)/((c+del*bn)) - 0.5*(del*an*c)/((c+del*bn)^2);
% f5 = bn == b+(c*del^2*bn^2)/((c+del*bn)^2) + (del*bn*c^2)/((c+del*bn)^2);
% f6 = gn == del*gn - 0.5*((c+1)*del^2*an^2)/((c+del*bn)^2) + (del^2*an^2)/(c+del*bn);
% 
% solns2 = solve([f4 f5 f6],[an bn gn]);

pa = 2:20;
pb = pa;%2*[1 2 4 8];
pc = pa;%2*[1 2 4 8];
pd = [0.3 0.5 0.8 0.9];
T = 500;
alph = zeros(length(pa),length(pb),length(pc),length(pd),T);
bet = alph;
gam = alph;
for a = 1:size(pa)
    for b = 1:length(pb)
        for c=1:length(pc)
            for d= 1:length(pd)
for tt=2:T
    alph(a,b,c,d,tt) = pa(a) + (pc(c)*pd(d)*alph(a,b,c,d,tt-1))/(pc(c)+bet(a,b,c,d,tt-1)*pd(d));
    bet(a,b,c,d,tt)  = pb(b) + (pc(c)*pd(d)*bet(a,b,c,d,tt-1))/(pc(c)+bet(a,b,c,d,tt-1)*pd(d));
    gam(a,b,c,d,tt)  = pd(d)*gam(a,b,c,d,tt-1) + 0.5*(alph(a,b,c,d,tt-1)^2*pd(d)^2)/(pc(c)+bet(a,b,c,d,tt-1)*pd(d));
%     alph(tt) = pa+(pc*pd^2*alph(tt-1)*bet(tt-1))/((pc+pd*bet(tt-1))^2) + (pd*alph(tt-1)*pc)/((pc+pd*bet(tt-1))) - 0.5*(pd*alph(tt-1)*pc)/((pc+pd*bet(tt-1))^2);
%     bet(tt)  = pb+(pc*pd^2*bet(tt-1)^2)/((pc+pd*bet(tt-1))^2) + (pd*bet(tt-1)*pc^2)/((pc+pd*bet(tt-1))^2);
%     gam(tt)  = pd*gam(tt-1) - 0.5*((pc+1)*pd^2*alph(tt-1)^2)/((pc+pd*bet(tt-1))^2) + (pd^2*alph(tt-1)^2)/(pc+pd*bet(tt-1));
end
            end
        end
    end
end


figure
subplot(3,4,1)
plot(squeeze(alph(:,1,1,1,end)))
title('\alpha vs a')
subplot(3,4,2)
plot(squeeze(alph(1,:,1,1,end)))
title('\alpha vs b')
subplot(3,4,3)
plot(squeeze(alph(1,1,:,1,end)))
title('\alpha vs c')
subplot(3,4,4)
plot(squeeze(alph(1,1,1,:,end)))
title('\alpha vs \delta')
subplot(3,4,5)
plot(squeeze(bet(:,1,1,1,end)))
title('\beta vs a')
subplot(3,4,6)
plot(squeeze(bet(1,:,1,1,end)))
title('\beta vs b')
subplot(3,4,7)
plot(squeeze(bet(1,1,:,1,end)))
title('\beta vs c')
subplot(3,4,8)
plot(squeeze(bet(1,1,1,:,end)))
title('\beta vs \delta')
subplot(3,4,9)
plot(squeeze(gam(:,1,1,1,end)))
title('\gamma vs a')
subplot(3,4,10)
plot(squeeze(gam(1,:,1,1,end)))
title('\gamma vs b')
subplot(3,4,11)
plot(squeeze(gam(1,1,:,1,end)))
title('\gamma vs c')
subplot(3,4,12)
plot(squeeze(gam(1,1,1,:,end)))
title('\gamma vs \delta')
return


figure
plot(alph)
title('\alpha')
figure
plot(bet)
title('\beta')
figure
plot(gam)
title('\gamma')
