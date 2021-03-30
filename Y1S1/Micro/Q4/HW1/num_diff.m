mkdir('pings')
clear; close all; clc

del = 0.01; % width of cgrid
cg = 0:del:1-del; % cgrid
rgrid = linspace(0,1,1e5)';
nc = length(cg);
r1 = 0*cg; % prealloc
r2 = r1;
p1 = r1;
p2 = r1;
for i=1:nc
    %[p1(i),r1(i)] = min(pr1(rgrid,cg(i)));
    [p2(i),r2(i)] = min(pr2(rgrid,cg(i)));
end
c = 0.05:0.05:1-0.05;
r1 = zeros(1,length(c));
p1 = r1;
for ic=1:length(c)
disp(['iteration ' num2str(ic) ', c = ' num2str(c(ic))])
grid = linspace(c(ic),1,1000);
b10 = grid;
b20 = b10;
tol = 1e-10;
maxiter = 1e8;
tune = 0.05;
ndraws = 5e6;
[r1(ic),p1(ic)] = calc_opt_r(grid,c(ic),b10,b20,tol,maxiter,tune,ndraws);
end

%p1 = -1*p1; % convert from negative of profit to profit
p2 = -1*p2;
%r1 = rgrid(r1);
r2 = rgrid(r2);
figure
plot(c,r1,'b--')
hold on
plot(cg,r2,'r-.')
plot(cg,cg,'k-')
hold off
title('(Expected) Profit Maximizing r')
ylabel('Optimal r')
xlabel('c')
legend('First-price','Second-price','c','Location','SouthEast')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'argmax.png')
cd('..')
figure
plot(c,p1,'b--')
hold on
plot(cg,p2,'r-.')
hold off
title('Expected Profit')
legend('First-price','Second-price','Location','NorthEast')
ylabel('Expected profit')
xlabel('c')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'profit.png')
cd('..')

function f = pr2(r,c) % negative profit so optimum is the minimum - to use with fmincon
f = -1*((r-c).*(r+r.^2-2*r.^3) + (1/3).*(((1/6)*(3*r.^4 - 4*r.^3 + 1) - c).*((r-1).^2.*(2*r+1)) + ((1/4)*(r.^2 - 1).^2 - c).*(r.^3 - 3*r+2))) ;
end