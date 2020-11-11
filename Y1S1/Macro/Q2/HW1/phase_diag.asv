mkdir('pings')
clear; close all; clc;

del = 0.1;
R = 1.05;
kst = 5;
num = 1000;
kgrid = 0:0.1:10;
delk0 = del.*kgrid;
delI0 = (kst)./(R+del) - kgrid.*(1-del)/(R+del);
ob = @(x) abs(R*del*x - (1-del)*del*x + x - kst);
ssk0 = fminunc(ob,4);
ssI0 = del* ssk0;
saddle0 = calc_saddle_p3(kgrid,kst,num,ssk0,ssI0,del,R);
delk1 = delk0;
R1 = 1.1;
delI1 = (kst)./(R1+del) - kgrid.*(1-del)/(R1+del);
ob = @(x) abs(R1*del*x - (1-del)*del*x + x - kst);
ssk1 = fminunc(ob,4);
ssI1 = del*ssk1;
saddle1 = calc_saddle_p3(kgrid,kst,num,ssk1,ssI1,del,R1);

figure
plot(kgrid,delk0,'k-')
hold on
plot(kgrid,delI0,'b-')
plot(kgrid,saddle0,'r-')
plot(kgrid,0*kgrid,'k')
hold off
xlabel('k')
ylabel('I')
ylim([0 5])
legend('\Delta k = 0','\Delta I = 0','saddle path')
set(gcf,'Color',[1 1 1])
title('Phase diagram')
% draw arrows
annotation('arrow',[6 7]/10,[3 3]/5)
annotation('arrow',[6 6]/10,[3 3.5]/5)
annotation('arrow',[2 3]/10,[1.5 1.5]/5)
annotation('arrow',[2 2]/10,[1.5 1]/5)
annotation('arrow',[8 7]/10,[0.25 0.25]/5)
annotation('arrow',[8 8]/10,[0.25 0.75]/5)
annotation('arrow',[5 4]/10,[0.75 0.75]/5)
annotation('arrow',[5 5]/10,[0.75 0.25]/5)
cd('pings')
saveas(gcf,'phase3.png')
cd('..')

figure
plot(kgrid,delk0,'k-')
hold on
plot(kgrid,delI0,'b-')
plot(kgrid,saddle0,'r-')
plot(kgrid,delI1,'b-.')
plot(kgrid(saddle1>0),saddle1(saddle1>0),'r-.')
plot(kgrid,0*kgrid,'k')
hold off
xlabel('k')
ylabel('I')
ylim([0 5])
legend('\Delta k = 0','\Delta I = 0','saddle path','\Delta I = 0 (high R)','saddle path (high R)')
set(gcf,'Color',[1 1 1])
title('Phase diagram comparison')
cd('pings')
saveas(gcf,'phasecomp3.png')
cd('..')

pathk = [ssk0 ssk0 ssk1];
pathI = [ssI0 interp1(kgrid,saddle1,ssk0) ssI1];

figure
plot(kgrid,delk0,'k-')
hold on
plot(kgrid,delI0,'b-')
plot(kgrid,saddle0,'r-')
plot(kgrid,delI1,'b-.')
plot(kgrid(saddle1>0),saddle1(saddle1>0),'r-.')
plot(pathk,pathI,'m-')
plot(pathk(1),pathI(1),'mo')
plot(pathk(2),pathI(2),'m+')
plot(pathk(3),pathI(3),'mx')
plot(kgrid,0*kgrid,'k')
hold off
xlabel('k')
ylabel('I')
ylim([0.4 0.55])
xlim([4.8 5.1])
legend('\Delta k = 0','\Delta I = 0','saddle path','\Delta I = 0 (high R)','saddle path (high R)','transition')
set(gcf,'Color',[1 1 1])
title('Transition dynamics')
cd('pings')
saveas(gcf,'transition3.png')
cd('..')
