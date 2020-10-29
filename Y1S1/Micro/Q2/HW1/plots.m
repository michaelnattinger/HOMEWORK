mkdir('pings')
clear; close all; clc

grid = 0:0.0001:1;

%plot 1
BR1 = grid<=0.6;
BR2 = grid<=(8/11);
figure
plot(BR1,grid,'r-')
hold on
plot(grid,BR2,'b-')
plot([0 1 8/11],[1 0 0.6],'ko')
hold off
set(gcf,'Color',[1 1 1])
xlim([0 1.1])
ylim([0 1.1])
ylabel('\sigma_2')
xlabel('\sigma_1')
legend('BR1','BR2','Equilibria')
title('Best response curves')
cd('pings')
saveas(gcf,'fig1.png')
cd('..')

%plot 2
R = 10;
r = 2;
c = 5;
BR1 = grid<=(R-c)/(R-r);
BR2 = grid<=(R-c)/(R-r);
figure
plot(BR1,grid,'r-')
hold on
plot(grid,BR2,'b-')
plot([0 1 (R-c)/(R-r)],[1 0 (R-c)/(R-r)],'ko')
hold off
set(gcf,'Color',[1 1 1])
xlim([0 1.1])
ylim([0 1.1])
ylabel('p_2')
xlabel('p_1')
legend('BR1','BR2','Equilibria')
set(gca,'xtick',[],'ytick',[])
title('Best response curves')
cd('pings')
saveas(gcf,'fig2.png')
cd('..')