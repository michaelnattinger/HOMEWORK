mkdir('pings')
clear; close all; clc
grid = 0:0.01:10;
figure
p=plot([1 (1/3) (2/3) 0 0],[0 (1/3) (2/3) 2 100],'r');
p.LineWidth = p.LineWidth*2;
hold on
plot(1/3,1/3,'r+')
plot(2/3,2/3,'rx')
hold off
xlim([0 3])
ylim([0 3])
set(gcf,'Color',[1 1 1])
title('Optimal consumption choices')
legend('optimal consumption choice','(1/3,1/3)','(2/3,2/3)')
cd('pings')
saveas(gcf,'indiff2.png')
cd('..')
figure
plot([1 (1/3) (2/3) 0 0]-1,[0 (1/3) (2/3) 2 100],'r')
hold on
plot([-10 10],[0 0],'k')
plot([0 0],[-10 10],'k')
p=plot([1 (1/3) (2/3) 0 0]-1,[0 (1/3) (2/3) 2 100],'r');
p.LineWidth = p.LineWidth*2;
hold off
xlim([0 3]-1)
ylim([0 3])
set(gcf,'Color',[1 1 1])
title('Offer Curve')
cd('pings')
saveas(gcf,'off2.png')
cd('..')

figure
p=plot([0 0 4 7 21 100],[100 12 4 7 0 0],'r');
p.LineWidth = p.LineWidth*2;
hold on
plot(4,4,'r+')
plot(7,7,'rx')
hold off
xlim([0 25])
ylim([0 15])
set(gcf,'Color',[1 1 1])
title('Optimal consumption choices')
legend('optimal consumption choice','(4,4)','(7,7)')
cd('pings')
saveas(gcf,'indiff3.png')
cd('..')
figure
plot([0 0 4 7 21 100] - 1,[100 12 4 7 0 0] - 10,'r');
hold on
plot([-100 100],[0 0],'k')
plot([0 0],[-100 100],'k')
p=plot([0 0 4 7 21 100] - 1,[100 12 4 7 0 0] - 10,'r');
p.LineWidth = p.LineWidth*2;
hold off
xlim([0 25] - 1)
ylim([0 15] - 10)
set(gcf,'Color',[1 1 1])
title('Offer Curve')
cd('pings')
saveas(gcf,'off3.png')
cd('..')