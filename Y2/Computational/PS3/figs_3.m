mkdir('pings')
clear; close all; clc

dat = dlmread('pfs_.dat');

v_50 = dat(:,1);
a_20_h = dat(:,2);
a_20_l = dat(:,3);
x = linspace(0,75,length(v_50));
cd('pings')
figure
plot(x,v_50,'k')
set(gcf,'Color',[1 1 1])
title('Value function (age 50) vs asset holdings')
saveas(gcf,'value.png')
figure
plot(x,a_20_h,'k')
hold on
plot(x,a_20_l,'b')
hold off
set(gcf,'Color',[1 1 1])
title('A prime (age 20) vs asset holdings')
legend('High prod','Low prod')
figure
plot(x,a_20_h - x','k')
hold on
plot(x,a_20_l - x','b')
hold off
set(gcf,'Color',[1 1 1])
title('savings (age 20) vs asset holdings')
legend('High prod','Low prod')
saveas(gcf,'savings.png')
cd('..')
close all
dat_tab = dlmread('pfs_K.dat');
dat_tab = dat_tab';
rows = {'K' 'L' 'w' 'r' 'b' 'W' 'cv'}';
tab = table(rows,dat_tab(:,1),dat_tab(:,2),dat_tab(:,3),dat_tab(:,4), ...
    dat_tab(:,5),dat_tab(:,6),'VariableNames',{'Variable' 'Bench\_SS' 'Bench\_noSS' 'NoRisk\_SS' 'NoRisk\_noSS' 'Labor\_SS' 'Labor\_noSS'});
table2latex(tab,'table1',4)