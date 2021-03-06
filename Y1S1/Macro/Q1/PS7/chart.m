mkdir('pings')
clear; close all; clc
y=1;
pbeta = 0.999;
palpha = 0.99999;
gspac = 0.01;
pt = 0:gspac:y-gspac;
pt1 =pt./(pbeta*(y - pt)) - palpha/pbeta;

figure
plot(pt,pt1,'b-')
hold on
plot([y y],[-10 100],'r--')
plot([palpha/(1+palpha) palpha/(1+palpha)],[-10 100],'r--')
plot(pt(pt1>0),pt1(pt1>0),'r-')
plot([-100 100],[0 0],'k-')
plot([0 0],[-100 100],'k-')
ylim([-5 20])
xlim([-0.2 1.2])
hold off
title('P_{t+1} vs P_t')
ylabel('P_{t+1}')
xlabel('P_t')
set(gca,'xtick',[],'ytick',[])
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'fig1.png')
cd('..')