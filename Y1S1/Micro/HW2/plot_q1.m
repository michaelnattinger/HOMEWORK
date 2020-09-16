mkdir('pings')
clear; close all; clc
y1 = 0:-0.01:-100;
y2 = (-1*y1).^(2/3);
figure
h =fill([y1 y1(end) 0],[y2 -25 -25],[0.95 0.95 0.95]);
hold on
h(1).EdgeColor = [0.95 0.95 0.95];
plot([0 0],[-100 100],'k')
plot([-100 100],[0 0],'k')
plot(y1,y2,'k')
hold off
xlim([-100 25])
ylim([-25 25])
set(gcf,'Color',[1 1 1])
title('Y (shaded)')
xlabel('y_1')
ylabel('y_2')
cd('pings')
saveas(gcf,'Y.png')
cd('..')