clear; close all; clc

kgrid = 0.01:0.01:5;
p = parameters();
[mod,ss,lin] = calcmodss(p,0,[]);
p2 = p;
p2.pz = p.pz + 0.1;
[mod2,ss2,lin2] = calcmodss(p2,0,[]);
f = p.pz*kgrid.^p.palpha;
delk0 = f - kgrid*p.pdelta;
delc0 = f + (1 - p.pdelta)*kgrid - ss.y.k;
f2 = p2.pz*kgrid.^p2.palpha;
delk02 = f2 - kgrid*p2.pdelta;
delc02 = f2 + (1 - p2.pdelta)*kgrid - ss2.y.k;
sad2 = (lin2.sseigvecs(2,2)/lin2.sseigvecs(1,2))*(kgrid - ss2.y.k) + ss2.y.c;
sad = (lin.sseigvecs(2,2)/lin.sseigvecs(1,2))*(kgrid - ss.y.k) + ss.y.c;

figure
plot(kgrid,delk0,'k-.')
hold on
plot(kgrid,delc0,'b-.')
plot(kgrid,delk02,'k-')
plot(kgrid,delc02,'b-')
plot(ss.y.k,ss.y.c,'rx')
plot(ss2.y.k,ss2.y.c,'r+')
plot(ss.y.k,1.1762,'ro')
plot([ss.y.k ss.y.k],[ss.y.c 1.1762],'r:')
plot(kgrid,sad,'m-.')
plot(kgrid,sad2,'m-')
plot([ss.y.k ss2.y.k],[1.1762 ss2.y.c],'r-')
hold off
set(gcf,'Color',[1 1 1])
title('Phase diagram')
ylim([1 1.5])
xlim([3 4])
xlabel('K')
ylabel('C')
annotation('arrow',[0.2 0.25],[0.17 0.17])
annotation('arrow',[0.2 0.2],[0.17 0.22])
annotation('arrow',[0.85 0.8],[0.7 0.7])
annotation('arrow',[0.85 0.85],[0.7 0.65])
annotation('arrow',[0.25 0.2],[0.65 0.65])
annotation('arrow',[0.25 0.25],[0.65 0.7])
annotation('arrow',[0.8 0.85],[0.22 0.22])
annotation('arrow',[0.8 0.8],[0.22 0.17])
cd('pings')
saveas(gcf,'phase.png')
cd('..')