clear; close all; clc
mkdir('pings')
pts0 = [-20 40; -40 70; -70 90; 0 0];
prices = [7 4; 5 5; 4 8];
[lnp,itsct] = do_lines(pts0(1:end-1,1),pts0(1:end-1,2),prices);
% Slope Check
for i=1:3
%disp(lnp(i,1) +prices(i,1)/prices(i,2))
disp(['Line ' num2str(i) ' eqn is: ' num2str(lnp(i,1)) ' *x + ' num2str(lnp(i,2)) '.'])
disp(['Price vector ' num2str(i) ' is (' num2str(prices(i,1)) ',' num2str(prices(i,2)) ').'])
end

%grid = -100:0;
grid = -100:25;
line1 = lnp(1,1)*grid +lnp(1,2);
line2 = lnp(2,1)*grid +lnp(2,2);
line3 = lnp(3,1)*grid +lnp(3,2);
vecscale = 1;
figure
hold on
h = fill([grid(end) itsct(:,1)' grid(1) grid(1) grid(end)],[ line1(end) itsct(:,2)' line3(1) -25 -25],[0.95 0.95 0.95]);
%h=fill([grid -100 0],[yg -25 -25],'k');
%h(1).FaceColor = [0.95 0.95 0.95];
h(1).EdgeColor = [0.95 0.95 0.95];
plot(grid,line1,'k')
plot(grid,line2,'k')
plot(grid,line3,'k')
for i=1:3
plot([pts0(i,1) pts0(i,1)+vecscale*prices(i,1)],[pts0(i,2) pts0(i,2)+vecscale*prices(i,2)],'b-')
end
%plot(grid,yg,'r');
plot(pts0(1:end-1,1),pts0(1:end-1,2),'rx')
xline(0)
yline(0)
hold off
set(gcf,'Color',[1 1 1])
title('Largest production set that can rationalize the data')
%legend('production function','Location','NorthEast')
xlim([-100 25])
ylim([-25 100])
cd('pings')
saveas(gcf,'largestprod.png') %savefig
cd('..')

figure
%plot(grid,yg,'k');
%hold on
%h=fill([grid -100 0],[yg -25 -25],'k');
%h(1).FaceColor = [0.95 0.95 0.95];
%plot(grid,yg,'k');
plot(pts0(1:end,1),pts0(1:end,2),'rx') % Need to make this convex!
hold on
h1 = fill([0 pts0(1:end-1,1)' -100 -100 0],[0 pts0(1:end-1,2)' pts0(end-1,2) -25 -25],[0.95 0.95 0.95]);
% h2 = fill([pts0(2,1) -100 -100 pts0(2,1)],[pts0(2,2) pts0(2,2) -25 -25],'k');
% h3 = fill([pts0(3,1) -100 -100 pts0(3,1)],[pts0(3,2) pts0(3,2) -25 -25],'k');
% h4 = fill([pts0(4,1) -100 -100 pts0(4,1)],[pts0(4,2) pts0(4,2) -25 -25],'k');
% h1(1).FaceColor = [0.95 0.95 0.95];
% h2(1).FaceColor = [0.95 0.95 0.95];
% h3(1).FaceColor = [0.95 0.95 0.95];
% h4(1).FaceColor = [0.95 0.95 0.95];
h1(1).EdgeColor = [0.95 0.95 0.95];
plot([0; pts0(1:end-1,1)],[0; pts0(1:end-1,2)],'k') % Need to make this convex!
plot(pts0(1:end,1),pts0(1:end,2),'rx')
xline(0)
yline(0)
hold off
set(gcf,'Color',[1 1 1])
title('Free disposal, convex, shutdown property')
%legend('production function','Location','NorthEast')
xlim([-100 25])
ylim([-25 100])
cd('pings')
saveas(gcf,'freedispshutdown.png')
cd('..')