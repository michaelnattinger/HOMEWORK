mkdir('pings')
clear; close all; clc
A = [4 2; 2 4]; V = [3;5];
soln = A\V; % solution to 5

% graphs for 4

a = 15;
b1 = -1.5;
b2 = 2*b1;
xx = 0:0.01:5;
D = a+b1*xx;
MR = a + b2*xx;
figure
subplot(3,1,1)
plot(xx,D,'b')
hold on
plot(xx,MR,'r')
plot(xx,0*xx,'k')
hold off
title('4 A')
legend('D','MR','Location','SouthWest')

S = xx.^(0.3);
subplot(3,1,2)
plot(xx,S,'k')
title('4 B')
legend('S','Location','SouthEast')

Cutoff = 2;
Resume = 3;
[~,indi] = min(abs(xx-Cutoff));
[~,indf] = min(abs(xx-Resume));
D2 = D;
D2(indi:indf) = D(indi);
MR2 = MR;
MR2(indi:indf) = D(indi);
n = length(D);
m = n-indf;

D2(indf:end) = D(indi:indi+m);
MR2(indf:end) = MR(indi:indi+m);
MR2([indi indf]) = NaN; 

subplot(3,1,3)
plot(xx,D2,'b')
hold on
plot(xx,MR2,'r')
plot(xx,0*xx,'k')
hold off
title('4 C')
legend('D','MR','Location','SouthWest')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'q4.png')
cd('..')