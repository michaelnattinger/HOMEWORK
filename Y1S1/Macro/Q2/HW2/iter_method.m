clear; close all; clc
pbeta = 0.95;
pdelta = 0.1;
pz0 = 1;
pz1 = 1.2;
pgamma0 = 2;
pgamma1 = 1.01;
ppsi = 0.35;

k = 0.001:0.001:5; %k grid

[kss0,css0,c,kprime,v] = iter_meth(k,pbeta,pdelta,pz0,pgamma0,ppsi);

delk0 = pz0 * k.^(ppsi) - pdelta*k; %phase
delc0 = pz0 * k.^ppsi + (1-pdelta)*k - (((1/pbeta) -(1-pdelta))/(0.35*pz0))^(1/(ppsi-1)) ;

figure
plot(k,c,'k')
hold on
plot(k,delk0,'b-')
plot(k,delc0,'r-')
plot(kss0,css0,'ko')
hold off
set(gcf,'Color',[1 1 1])
legend('Saddle path (c(k))','\Delta k = 0','\Delta c = 0','steady state','Location','SouthEast')
title('Computational solution to part (a)')
xlabel('k')
ylabel('c')
annotation('arrow',[0.2 0.2 ],[0.5 0.55])
annotation('arrow',[0.2 0.25 ],[0.5 0.5])
annotation('arrow',[0.87 0.87 ],[0.8 0.75])
annotation('arrow',[0.87 0.82 ],[0.8 0.8])
annotation('arrow',[0.5 0.5 ],[0.85 0.9])
annotation('arrow',[0.5 0.45 ],[0.85 0.85])
annotation('arrow',[0.65 0.65 ],[0.5 0.45])
annotation('arrow',[0.65 0.7 ],[0.5 0.5])
cd('pings')
saveas(gcf,'partA.png')
cd('..')

figure
subplot(2,1,1)
plot(k,v)
title('Value function V(k)')
xlabel('k')
ylabel('V(k)')
subplot(2,1,2)
plot(k,kprime)
title('k prime from value function (argmax)')
ylabel('kprime(k)')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'partA2.png')
cd('..')

delk1 = pz0 * k.^(ppsi) - pdelta*k; %phase
delc1 = pz0 * k.^ppsi + (1-pdelta)*k - (((1/pbeta) -(1-pdelta))/(0.35*pz0))^(1/(ppsi-1)) ;

[kss1,css1,c1,kprime1,v1] = iter_meth(k,pbeta,pdelta,pz0,pgamma1,ppsi);
figure
plot(k,c,'k')
hold on
plot(k,delk0,'b-')
plot(k,delc0,'r-')
plot(kss0,css0,'ko')
plot(k,c1,'k-.')
hold off
set(gcf,'Color',[1 1 1])
legend('Saddle path (c(k))','\Delta k = 0','\Delta c = 0','steady state','new saddle path (c(k))','Location','SouthEast')
title('Computational solution to part (b)')
xlabel('k')
ylabel('c')
cd('pings')
saveas(gcf,'partB.png')
cd('..')

figure
subplot(3,1,1)
plot(k,v)
hold on
%plot(k,v1)
hold off
legend('old','Location','SouthEast')
title('Value function V(k)')
subplot(3,1,2)
plot(k,v1)
hold on
%plot(k,v1)
hold off
legend('new','Location','SouthEast')
title('Value function V(k)')
xlabel('k')
ylabel('V(k)')
subplot(3,1,3)
plot(k,kprime)
hold on
plot(k,kprime1)
hold off
legend('old','new','Location','SouthEast')
title('k prime from value function (argmax)')
ylabel('kprime(k)')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'partB2.png')
cd('..')

delk2 = pz1 * k.^(ppsi) - pdelta*k; %phase
delc2 = pz1 * k.^ppsi + (1-pdelta)*k - (((1/pbeta) -(1-pdelta))/(0.35*pz1))^(1/(ppsi-1)) ;

[kss2,css2,c2,kprime2,v2] = iter_meth(k,pbeta,pdelta,pz1,pgamma0,ppsi);
figure
plot(k,c,'k')
hold on
plot(k,delk0,'b-')
plot(k,delc0,'r-')
plot(kss0,css0,'ko')
plot(k,c2,'k-.')
plot(k,delk2,'b-.')
plot(k,delc2,'r-.')
plot(kss2,css2,'kx')
hold off
set(gcf,'Color',[1 1 1])
legend('Saddle path (c(k))','\Delta k = 0','\Delta c = 0','steady state','new saddle path (c(k))','new \Delta k = 0','new \Delta c = 0','new SS','Location','SouthEast')
title('Computational solution to part (c)')
xlabel('k')
ylabel('c')
cd('pings')
saveas(gcf,'partC.png')
cd('..')

figure
subplot(2,1,1)
plot(k,v)
hold on
plot(k,v2)
hold off
legend('old','new')
title('Value function V(k)')
xlabel('k')
ylabel('V(k)')
subplot(2,1,2)
plot(k,kprime)
hold on
plot(k,kprime2)
hold off
legend('old','new','Location','SouthEast')
title('k prime from value function (argmax)')
ylabel('kprime(k)')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'partC2.png')
cd('..')

nt=50;
%find saddle path for transition dynamics
cc0 = css0:0.00001:css2;
kk = 0*cc0+kss0;
cc=cc0;
for i=1:nt
    [kk,cc] = LOM(kk,cc,pz1,ppsi,pdelta,pbeta,pgamma0);
end
[~,ind] = min(abs(cc - css2)+abs(kk-kss2));
cc0 = cc0(ind);
ctraj = cc0*ones(1,nt);
ktraj = kss0*ones(1,nt);
for i=2:nt
    [ktraj(i),ctraj(i)] = LOM(ktraj(i-1),ctraj(i-1),pz1,ppsi,pdelta,pbeta,pgamma0);
end

ct = [css0 ctraj css2];
kt = [kss0 ktraj kss2];
figure
plot(kt,ct,'k')
hold on
plot(kt(1),ct(1),'ko')
plot(kt(2),ct(2),'k+')
plot(kt(end),ct(end),'kx')
hold off
set(gcf,'Color',[1 1 1])
title('Computational solution to particle trajectory for part (c)')
xlabel('k')
ylabel('c')
legend('Particle trajectory','initial steady state','jump','final steady state','Location','SouthEast')
cd('pings')
saveas(gcf,'traj.png')
cd('..')

function [kk,cc] = LOM(k,c,pz,ppsi,pdelta,pbeta,pgamma)
kk = pz*k.^ppsi - c + (1-pdelta)*k;
cc = c.*(pbeta*(ppsi*pz*k.^(ppsi-1)+1-pdelta)).^(1/pgamma);
end