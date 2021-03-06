mkdir('pings')
clear; close all; clc
%% First section: prepare data
rerunclean = 0;
if rerunclean
reload = 0;
if reload
[x,xt] = xlsread('macropset3q3.xls','matlab');
[I,it] = xlsread('macropset3q3.xls','inv');
save 'raw'
else
load 'raw'
end 
names = xt(1,2:end);
dates = datetime(xt(2:end,1));
dateI = datetime(it(2:end,1));
% Question 1 - plot raw x,I
figure
for i=1:3
subplot(3,1,i)
plot(dates,x(:,i),'k')
title(names{i})
end
set(gcf,'Color',[1 1 1])
suptitle('Raw data')
cd('pings')
saveas(gcf,'raw.png')
cd('..')

figure
plot(dateI,I,'k');
set(gcf,'Color',[1 1 1])
suptitle('Raw investment data')
cd('pings')
saveas(gcf,'rawi.png')
cd('..')

% Question 2 - plot logged x,I and log detrended x,I
x = log(x);
I = log(I);
% hp filter
xhp = hpfilter(x); % note: default smoothing is 1600 which is correct for quarterly data
Ihp = hpfilter(I);

figure
for i=1:3
subplot(3,1,i)
plot(dates,x(:,i),'k')
hold on
plot(dates,xhp(:,i),'r')
hold off
title(names{i})
end
legend('log data','HP trend')
set(gcf,'Color',[1 1 1])
suptitle('Log data with trend')
cd('pings')
saveas(gcf,'log.png')
cd('..')


figure
plot(dateI,I,'k');
hold on
plot(dateI,Ihp,'r');
hold off
set(gcf,'Color',[1 1 1])
suptitle('Log investment data with trend')
cd('pings')
saveas(gcf,'logi.png')
cd('..')


x = x-xhp;
I = I-Ihp;

figure
for i=1:3
subplot(3,1,i)
plot(dates,x(:,i),'k')
title(names{i})
ylabel('log deviation from trend')
end
set(gcf,'Color',[1 1 1])
suptitle('Log detrended data')
cd('pings')
saveas(gcf,'det.png')
cd('..')

figure
plot(dateI,I,'k');
ylabel('log deviation from trend')
set(gcf,'Color',[1 1 1])
suptitle('Log detrended investment data')
cd('pings')
saveas(gcf,'deti.png')
cd('..')

close all
% Question 3 - calculate k and plot
pdelta = 0.025;
k = 0*I;
for tt=2:length(k)
    k(tt) = (1-pdelta)*k(tt-1) + pdelta*I(tt-1);
end
T = size(x,1);
k = k(end-T+1:end);
figure
plot(dates,k,'k')
set(gcf,'Color',[1 1 1])
ylabel('log deviation from trend')
suptitle('Capital approximation')
cd('pings')
saveas(gcf,'detk.png')
cd('..')

close all
I = I(end-T+1:end);
y = x(:,1);
c = x(:,2);
l = x(:,3);
save 'clean'
else
load 'clean'
end

%% Second section: back out a,g,tauL
palpha = (1/3);
psigma = 1;
pphi = 1;
pGbar = (1/3); % as ratio of Y
pAbar = 1;
ptaubarL = 0;
ptaubarI = 0;
pbeta = 0.99; 
a = y - palpha*k - (1-palpha)*l; 
[Ybar,Cbar,Kbar,Lbar] = calc_ss(palpha, pbeta, pdelta, psigma, pphi, pAbar, pGbar, ptaubarI, ptaubarL);
Ibar = pdelta*Kbar;
g = (1/pGbar)*(y - (Cbar/Ybar)*c - (Ibar/Ybar)*I);
htauL = (-1)*(pphi*l + psigma*c - palpha*k + palpha * l);
% calculate shock persistance
rhoa = a(1:end-1)\a(2:end); % just ols ar(1) with no intercept
rhog = g(1:end-1)\g(2:end);
rhoL = htauL(1:end-1)\htauL(2:end);

%% Blanchard-Kahn implementation
rhoI0 = 0;
tol = 1e-6;
rebkfp = 1;
if rebkfp
[rhoI,htauI] = BKFP(rhoI0,tol,rhoa,rhog,rhoL,a,g,htauL,c,k,palpha,pdelta,psigma,pphi,pGbar, ...
    pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
save 'bkfp' 'rhoI' 'htauI'
else
load 'bkfp'
end

% rho table
rho = [rhoa rhog rhoL rhoI]';
tab = table(rho,'rowNames',{'a' 'g' 'tau L' 'tau I'},'VariableNames',{'rho'});
table2latex(tab,'rhotable.tex')
% and figure
figure
plot(dates,a,'k')
hold on
plot(dates,g,'r')
plot(dates,htauL,'b')
plot(dates,htauI,'m')
hold off
set(gcf,'Color',[1 1 1])
title('Wedges since 1980')
legend('a','g','\tau_L','\tau_I')
cd('pings')
saveas(gcf,'wedges.png')
cd('..')

figure
subplot(2,2,1)
plot(dates,a,'k')
title('a')
subplot(2,2,2)
plot(dates,g,'k')
title('g')
subplot(2,2,3)
plot(dates,htauL,'k')
title('\tau_L')
subplot(2,2,4)
plot(dates,htauI,'k')
title('\tau_I')
set(gcf,'Color',[1 1 1])
suptitle('Wedges')
cd('pings')
saveas(gcf,'wedges2.png')
cd('..')
close all
%% Counterfactual
[ga,gg,ghtauL,ghtauI] = BK_counterfac(rhoI,rhoa,rhog,rhoL,a,g,htauL,htauI,c,k,palpha,pdelta,psigma,pphi,pGbar, ...
    pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
% and plots
figure
plot(dates,ga,'k-.')
hold on
plot(dates,gg,'r')
plot(dates,ghtauL,'b--')
plot(dates,ghtauI,'m-.')
plot(dates,y,'g-.')
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
cd('pings')
saveas(gcf,'wedgesg.png')
cd('..')

date2007 = dates(datenum(dates)>datenum('12/31/2007'));
inds = find(datenum(dates)>datenum('12/31/2007'));
date2007 = datetime(date2007(datenum(date2007)<datenum('12/31/2009')));
inds = inds(datenum(date2007)<datenum('12/31/2009'));
figure
plot(dates,ga,'k-.')
hold on
plot(dates,gg,'r')
plot(dates,ghtauL,'b--')
plot(dates,ghtauI,'m-.')
plot(dates,y,'g-.')
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2007(1) date2007(end)])
cd('pings')
saveas(gcf,'wedgesfin.png')
cd('..')

dga = ga(inds,:) - ga(inds(1),:);
dgg = gg(inds,:) - gg(inds(1),:);
dghtauL = ghtauL(inds,:) - ghtauL(inds(1),:);
dghtauI = ghtauI(inds,:) - ghtauI(inds(1),:);
dy = y(inds,:) - y(inds(1),:);
figure
plot(date2007,dga,'k-.')
hold on
plot(date2007,dgg,'r')
plot(date2007,dghtauL,'b--')
plot(date2007,dghtauI,'m-.')
plot(date2007,dy,'g-.')
hold off
set(gcf,'Color',[1 1 1])
ylabel('Net change since 2008')
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2007(1) date2007(end)])
cd('pings')
saveas(gcf,'wedgesfindiff.png')
cd('..')



date2020 = dates(datenum(dates)>datenum('12/31/2019'));
inds = find(datenum(dates)>datenum('12/31/2019'));
date2020 = datetime(date2020(datenum(date2020)<datenum('12/31/2020')));
inds = inds(datenum(date2020)<datenum('12/31/2020'));

figure
plot(dates,ga,'k-.')
hold on
plot(dates,gg,'r')
plot(dates,ghtauL,'b--')
plot(dates,ghtauI,'m-.')
plot(dates,y,'g-.')
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2020(1) date2020(end)])
cd('pings')
saveas(gcf,'wedgescov.png')
cd('..')

dga = ga(inds,:) - ga(inds(1),:);
dgg = gg(inds,:) - gg(inds(1),:);
dghtauL = ghtauL(inds,:) - ghtauL(inds(1),:);
dghtauI = ghtauI(inds,:) - ghtauI(inds(1),:);
dy = y(inds,:) - y(inds(1),:);
figure
plot(date2020,dga,'k-.')
hold on
plot(date2020,dgg,'r')
plot(date2020,dghtauL,'b--')
plot(date2020,dghtauI,'m-.')
plot(date2020,dy,'g-.')
hold off
set(gcf,'Color',[1 1 1])
title('Wedge effects on GDP')
ylabel('Net change since 2020')
legend('a','g','\tau_L','\tau_I','Actual GDP','Location','SouthWest')
xlim([date2020(1) date2020(end)])
cd('pings')
saveas(gcf,'wedgescovdiff.png')
cd('..')

close all

% test to ensure everything is right
[gtest] = BK_counterfac_test(rhoI,rhoa,rhog,rhoL,a,g,htauL,htauI,c,k,palpha,pdelta,psigma,pphi,pGbar, ...
    pAbar,ptaubarI,ptaubarL,pbeta,Ybar,Kbar,Cbar,Lbar);
figure
plot(dates,gtest,'k')
hold on
plot(dates,y,'r')
hold off
title('test comparison')
set(gcf,'Color',[1 1 1])