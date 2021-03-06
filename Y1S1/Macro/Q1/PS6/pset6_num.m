mkdir('pings')
clear; close all; clc

alpha = 0.01:0.01:0.49;
tg = 0*alpha;
objfn = @(t,a) -2*(1-t+a)^(-1) + 2*(1-t)^(-1) + (-3*a + 1 - 2*t)/(2*a*(1-t) + t*(1-t-a));
for i=1:length(alpha)
    tempobj = @(t)abs(objfn(t,alpha(i)));
    tg(i) = fmincon(tempobj,0.1,[1;-1],[1;0]);
end
% cg = (1-tg).*(.5 - (alpha)./(2*(1-tg)));
lg = 0.5+ (alpha./(2*(1-tg)));
cg = (1-tg).*(1-lg);
gg = tg.*(1-lg);
ug = log(lg) + log(alpha + cg) + log(alpha + gg);
figure
plot(alpha,tg,'k')
hold on
plot(alpha,cg,'b--')
plot(alpha,lg,'r-.')
plot(alpha,gg,'g-x')
hold off
set(gcf,'Color',[1 1 1])
title('Ramsey equilibrium')
legend('\tau','c','l','g','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'ramsey.png')
cd('..')

% calculate spp
csp = (1 - alpha)./3;
gsp = csp;
lsp = (2*alpha + 1)./3;
usp = log(lsp) + log(alpha+csp)+log(alpha+gsp);
tsp = gsp./(1-lsp);

% calculate ne
cne = 0.25 - 0.5*alpha;
gne = cne;
lne = alpha + 0.5;
une = log(lne) + log(alpha+cne) + log(alpha+gne);
tne = gne./(1-lne);

figure
plot(alpha,csp,'k')
hold on
plot(alpha,cg,'b--')
plot(alpha,cne,'r-.')
hold off
set(gcf,'Color',[1 1 1])
title('Consumption')
legend('SPP','Ramsey','Nash','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'cons.png')
cd('..')

figure
plot(alpha,lsp,'k')
hold on
plot(alpha,lg,'b--')
plot(alpha,lne,'r-.')
hold off
set(gcf,'Color',[1 1 1])
title('Leisure')
legend('SPP','Ramsey','Nash','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'leisure.png')
cd('..')

figure
plot(alpha,gsp,'k')
hold on
plot(alpha,gg,'b--')
plot(alpha,gne,'r-.')
hold off
set(gcf,'Color',[1 1 1])
title('Gov')
legend('SPP','Ramsey','Nash','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'gov.png')
cd('..')

figure
plot(alpha,usp,'k')
hold on
plot(alpha,ug,'b--')
plot(alpha,une,'r-.')
hold off
set(gcf,'Color',[1 1 1])
title('Utility')
legend('SPP','Ramsey','Nash','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'util.png')
cd('..')

figure
plot(alpha,tsp,'k')
hold on
plot(alpha,tg,'b--')
plot(alpha,tne,'r-.')
hold off
ylim([0 0.6])
set(gcf,'Color',[1 1 1])
title('\tau')
legend('SPP','Ramsey','Nash','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'tau.png')
cd('..')

%% calculate betabar
betabar = 0*tg;
for i=1:length(alpha)
    if ug(i)>une(i)
        temputil = log(lg(i)) + log(alpha(i) + (0.5)*(1-lg(i))) + log(alpha(i) + 0.5*(1-lg(i)));
        genobj = @(b) abs(ug(i)*(1/(1-b)) - temputil - (b/(1-b)*une(i)));
        betabar(i) = fmincon(genobj,0.1,[1;-1],[1;0]);
    else
        betabar(i) = 1;
    end
end
figure
plot(alpha,betabar)
set(gcf,'Color',[1 1 1])
title('\beta_L')
%legend('SPP','Ramsey','Nash','Location','NorthWest')
xlabel('\alpha')
cd('pings')
saveas(gcf,'beta.png')
cd('..')