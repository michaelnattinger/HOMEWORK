mkdir('pings')
clear; close all; clc

dat = dlmread('pfs_K.dat');
NT = size(dat,1);
KK = dat(:,1);
LL = dat(:,4);
WW = dat(:,5);
RR = dat(:,6);
KKK = dat(:,7);
LLL = dat(:,8);
WWW = dat(:,9);
RRR = dat(:,10);
cd('pings')
xx = (1:NT) - 1;
figure
plot(xx,KK,'b')
hold on
plot(xx,KKK,'r-.')
%plot(xx,0*xx + dat(:,2),'k')
%plot(xx,0*xx + dat(:,3),'k')
hold off
title('Capital transition')
legend('anticipated','(somewhat) unanticipated','Location','SouthEast')
set(gcf,'Color',[1 1 1])
saveas(gcf,'k.png')

figure
plot(xx,LL,'b')
hold on
plot(xx,LLL,'r-.')
hold off
title('Labor transition')
legend('anticipated','(somewhat) unanticipated')
set(gcf,'Color',[1 1 1])
saveas(gcf,'l.png')

figure
plot(xx,WW,'b')
hold on
plot(xx,WWW,'r-.')
hold off
title('Wage transition')
legend('anticipated','(somewhat) unanticipated','Location','SouthEast')
set(gcf,'Color',[1 1 1])
saveas(gcf,'w.png')

figure
plot(xx,RR,'b')
hold on
plot(xx,RRR,'r-.')
hold off
title('Rental price transition')
legend('anticipated','(somewhat) unanticipated')
set(gcf,'Color',[1 1 1])
saveas(gcf,'r.png')

cd('..')

% EV
EV = dat(1,11:12);
vote = dat(1,13:14);
mat = [EV;vote];
vars = {'EV'; 'Vote'};
tab=table(vars,mat(:,1),mat(:,2),'VariableNames',{' ' 'Experiment 1' 'Experiment 2'});
table2latex(tab,'tab1');