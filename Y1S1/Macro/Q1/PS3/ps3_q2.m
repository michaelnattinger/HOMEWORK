mkdir('pings')
clear; close all; clc
c1g = linspace(1.25 - sqrt(41)/4,1.25 + sqrt(41)/4,1000);
c2h = 2+ sqrt((41/4)*(1 - (16/41)*(c1g - (1.25)).^2));
c2l = 2 - (c2h - 2);%
forcepos =1;
if forcepos
    c1g = max(c1g, 0*c1g);
    c2h = max(c2h, 0*c2h);
    c2l = max(c2l, 0*c2l);
end
figure
hold on
plot([0 0],[-10 10],'k')
plot([-10 10],[0 0],'k')
p=patch([c1g fliplr(c1g)], [c2h fliplr(c2l)], 'k');
plot(c1g,c2h,'k')
plot(c1g,c2l,'k')
hold off
p.FaceColor = [0.95 0.95 0.95];
p.EdgeColor = p.FaceColor;
title('U \geq U(w_1,w_2)')
set(gcf,'Color',[1 1 1]);
xlabel('c_1')
ylabel('c_2')
xlim([-0.5 3])
ylim([-2 6])
cd('pings')
saveas(gcf,'p1.png')
cd('..')

grid = -2:0.01:5;
eq1 = 1-grid.*2;
eq2 = (1-grid)/2;
eqm = max(eq1,eq2);
pg = [grid grid(end) grid(1)];
py = [eqm 10 10];
if forcepos
pg = max(pg,0*pg);
py = max(py,0*py);
end
figure
hold on
p=patch(pg,py,'k');
plot(grid,eq1,'k')
p.FaceColor = [0.95 0.95 0.95];
p.EdgeColor = p.FaceColor;
plot([0 0],[-10 10],'k')
plot([-10 10],[0 0],'k')
plot(grid,eq2,'k')
hold off
xlim([-0.5 3])
ylim([-0.5 3])
title('U \geq U(w_1,w_2)')
set(gcf,'Color',[1 1 1]);
xlabel('c_1')
ylabel('c_2')
cd('pings')
saveas(gcf,'p2.png')
cd('..')

grid = -5:0.01:15;
eq1 = 12-grid.*2;
eq2 = (12-grid)/2;
eqm = max(eq1,eq2);
pg = [grid grid(end) grid(1)];
py = [eqm 20 20];
if forcepos
pg = max(pg,0*pg);
py = max(py,0*py);
end
figure
hold on
p=patch(pg,py,'k');
plot(grid,eq1,'k')
p.FaceColor = [0.95 0.95 0.95];
p.EdgeColor = p.FaceColor;
plot([0 0],[-10 20],'k')
plot([-10 20],[0 0],'k')
plot(grid,eq2,'k')
hold off
xlim([-2 15])
ylim([-2 15])
title('U \geq U(w_1,w_2)')
set(gcf,'Color',[1 1 1]);
xlabel('c_1')
ylabel('c_2')
cd('pings')
saveas(gcf,'p3.png')
cd('..')