clear; close all; clc
read = 0;
if read
[x,xt] = xlsread('cps09mar.xlsx','Sheet1');
save 'data'
else
load 'data'
end

female = x(:,2);
hisp = x(:,3);
edu = x(:,4);
ie = edu>=11;
I = logical(female.*hisp.*ie);
earnings = x(:,5);
hours = x(:,6);
week = x(:,7);
lw = log(earnings./(hours.*week));
Y = lw(I,:);
X = edu(I,:);
X = [ones(length(Y),1) X];
q = [0.05 0.25 0.5 0.75 0.95];
b=qr_nattinger(Y,X,q);
grid = (min(X(:,2)):0.1:max(X(:,2)))';
xg = [ones(size(grid)) grid];
Yh = xg*b;
colore = {'b--' 'r-.' 'k-' 'r-.' 'b--'};
figure
for i=1:length(q)
plot(grid,Yh(:,i),colore{i}); hold on
end
hold off
set(gcf,'Color',[1 1 1])
title('Quantile estimates')
xlabel('Education')
ylabel('Log wage')
legend('0.05', '0.25', '0.5', '0.75', '0.95','Location','NorthWest')
cd('pings')
saveas(gcf,'qrfig.png')
cd('..')