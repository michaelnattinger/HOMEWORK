mkdir('pings')
clear; close all; clc

dat = dlmread('pfs_K.dat');
NT = size(dat,1);
KK = dat(:,1);
xx = 1:NT;
figure
plot(xx,KK,'k')
hold on
plot(xx,0*xx + dat(:,2),'b--')
plot(xx,0*xx + dat(:,3),'r-.')
hold off
title('Capital transition')
set(gcf,'Color',[1 1 1])

