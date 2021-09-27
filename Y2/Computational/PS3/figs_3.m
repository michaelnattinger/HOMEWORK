mkdir('pings')
clear; close all; clc

dat = dlmread('pfs_.dat');

v_50 = dat(:,1);
a_20_h = dat(:,2);
a_20_l = dat(:,3);
x = linspace(0,75,length(v_50));
figure
plot(x,v_50,'k')
set(gcf,'Color',[1 1 1])
title('Value function (age 50) vs asset holdings')
figure
plot(x,a_20_h,'k')
hold on
plot(x,a_20_l,'b')
hold off
set(gcf,'Color',[1 1 1])
title('A prime (age 20) vs asset holdings')
legend('High prod','Low prod')

figure
plot(x,v_50,'k')
set(gcf,'Color',[1 1 1])
title('Value function (age 50) vs asset holdings')
figure
plot(x,a_20_h - x','k')
hold on
plot(x,a_20_l - x','b')
hold off
set(gcf,'Color',[1 1 1])
title('savings (age 20) vs asset holdings')
legend('High prod','Low prod')