mkdir('pings')
clear; close all; clc
d = 0.05:0.05:10;
%d = exp(d);
count=0;
for a = [0.3 0.6 0.9]
    count=count+1;
APA = d./(2-a);
BPA = APA;
APRA = APA.*(d - APA + a*BPA );%d.^2/(2-a)^2;
BPRA = BPA.*(d - BPA + a*APA );%APRA;
TPRA = APRA + BPRA;

APB = d*(2+a)./(4 - 2*a^2);
BPB = (d./2)*(1+a*(2+a)/(2*(2-a^2)));

APRB = APB.*(d - APB + a*BPB );
BPRB = BPB.*(d - BPB + a*APB );

figure
subplot(2,2,1)
plot(d,APA,'k')
hold on
plot(d,APB,'r--')
hold off
set(gcf,'Color',[1 1 1])
title('Firm A Price')
xlabel('d')
ylabel('Price')
subplot(2,2,2)
plot(d,BPA,'k')
hold on
plot(d,BPB,'r--')
hold off
set(gcf,'Color',[1 1 1])
title('Firm B Price')
xlabel('d')
ylabel('Price')
subplot(2,2,3)
plot(d,APRA,'k')
hold on
plot(d,APRB,'r--')
hold off
set(gcf,'Color',[1 1 1])
title('Firm A Profits')
xlabel('d')
ylabel('Profits')
subplot(2,2,4)
plot(d,BPRA,'k')
hold on
plot(d,BPRB,'r--')
hold off
set(gcf,'Color',[1 1 1])
title('Firm B Profits')
xlabel('d')
ylabel('Profits')
legend('Part A','Part B','Location','NorthWest')
suptitle(['Problem 6 results: \alpha = ' num2str(a,2)])
cd('pings')
saveas(gcf,['fig' num2str(count) '.png'])
cd('..')
end