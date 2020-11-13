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

pa = 1;
pb = 2;
pc = 3;
pd = 4;
T = 50000;
alph = zeros(1,T);
bet = alph;
gam = alph;

for tt=2:T
    alph(tt) = pa + (pc*pd*alph(tt-1))/(pc+bet(tt-1)*pd);
    bet(tt)  = pb + (pc*pd*bet(tt-1))/(pc+bet(tt-1)*pd);
    gam(tt)  = pd*gam(tt-1) + 0.5*(alph(tt-1)^2*pd^2)/(pc+bet(tt-1)*pd);
%     alph(tt) = pa+(pc*pd^2*alph(tt-1)*bet(tt-1))/((pc+pd*bet(tt-1))^2) + (pd*alph(tt-1)*pc)/((pc+pd*bet(tt-1))) - 0.5*(pd*alph(tt-1)*pc)/((pc+pd*bet(tt-1))^2);
%     bet(tt)  = pb+(pc*pd^2*bet(tt-1)^2)/((pc+pd*bet(tt-1))^2) + (pd*bet(tt-1)*pc^2)/((pc+pd*bet(tt-1))^2);
%     gam(tt)  = pd*gam(tt-1) - 0.5*((pc+1)*pd^2*alph(tt-1)^2)/((pc+pd*bet(tt-1))^2) + (pd^2*alph(tt-1)^2)/(pc+pd*bet(tt-1));
end

figure
plot(alph)
figure
plot(bet)
figure
plot(gam)
