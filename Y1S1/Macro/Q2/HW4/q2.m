mkdir('pings')
addpath('C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Macro\PS5')
clear; close all; clc
pbeta = 0.95;
pdelta = 0.1;
pz0 = 1;
pgamma0 = 2;
pgamma1 = 1.01;
ppsi = 0.35;
rerun = 0;
if rerun
k = 0.001:0.001:5; %k grid
[~,~,c0,kprime0,v0] = iter_meth(k,pbeta,pdelta,pz0,pgamma0,ppsi);
z = [0.8 1.2];
P = [0.9 0.1;0.1 0.9];
[c,kprime,v] = iter_meth_sto(k,pbeta,pdelta,z,P,pgamma0,ppsi);
save res_q2
else
    load res_q2
end
figure
plot(k,kprime(:,1),'k')
hold on
plot(k,kprime(:,2),'r')
plot(k,kprime0,'b')
hold off
set(gcf,'Color',[1 1 1])
legend('low prod','high prod','deterministic','Location','SouthEast')
title('k prime vs k')
ylabel('k prime')
xlabel('k')
cd('pings')
saveas(gcf,'q2k.png')
cd('..')

figure
plot(k,c(:,1),'k')
hold on
plot(k,c(:,2),'r')
plot(k,c0,'b')
hold off
set(gcf,'Color',[1 1 1])
legend('low prod','high prod','deterministic','Location','SouthEast')
title('c vs k')
ylabel('c')
xlabel('k')
cd('pings')
saveas(gcf,'q2c.png')
cd('..')

% simulate
dosim=0;
if dosim==1
T = 1e5;
burn = floor(T/10);
dist0 = ones(size(k'));
dist0 = dist0./sum(dist0(:));
stat = 1; % initial state
cons = zeros(T,1);
Y = cons;
kap = cons;
I = cons;
w = cons;
r = cons;
% make transition matrices
mat1 = zeros(length(k));
mat2 = mat1;
for i=1:length(k)
    mat1(:,i) = 0+kprime(i,1)==k;
    mat2(:,i) = 0+kprime(i,2)==k;
end
tic
for tt = 1:T+burn
    if mod(tt,burn)==0; toc; disp(['iteration ' num2str(tt)]); end
    draw = rand;
    if draw<0.1; stat = 1-stat; end % z transition with probability 0.1
    switch stat
        case 0
            dist = mat1*dist0;
        otherwise
            dist = mat2*dist0;
    end
    if tt>burn % after burn-in, calculate 
        t = tt-burn;
        % calculate values 
        cons(t) = sum(c(:,stat+1).*dist0);
        kap(t) = sum(k'.*dist0);
        I(t) =sum(kprime(:,stat+1).*dist0 - (1-pdelta)*k'.*dist0);
        Y(t) = kap(t)^(ppsi);
        r(t) = sum((ppsi)*kap(t).^(ppsi-1));
        w(t) = sum(kap(t).^(ppsi) - r(t)*kap(t));
    end
    dist0=dist;
end
save simres
else
    load simres
end
% calculate deterministic steady state
matss = 0*mat1;
for i=1:length(k)
    matss(:,i) = 0+kprime0(i)==k;
end
init=ones(size(k'));
init=init./(sum(init));
tol=1e-6;
diff = 999;
maxiter=1e6;
iter=1;
while (iter<maxiter)&&(diff>tol)
    next = matss*init;
    iter = iter+1;
    diff = sum(abs(next - init));
    init = next;
end
consss = sum(c0(:).*next);
kapss = sum(k'.*next);
Yss = kapss^(ppsi);
Iss =sum(kprime0'.*next - (1-pdelta)*k'.*next);
rss = sum((ppsi)*kapss.^(ppsi-1));
wss = sum(kapss.^(ppsi) - rss*kapss);
        

consa = mean(cons);
kapa = mean(kap);
Ia = mean(I);
ra = mean(r);
wa = mean(w);
Ya = mean(Y);

Deterministic = [ kapss Yss consss Iss wss rss rss-pdelta]';
Simulation = [kapa Ya consa  Ia wa ra ra-pdelta]';
tab = table(Deterministic,Simulation,'RowNames',{'Capital' 'Output' 'Consumption' 'Investment' 'wages' 'rental rate' 'Interest rate (rental net of depr.)'});
table2latex(tab,'q2.tex')

cvol = std(cons);
Yvol = std(Y);
Ivol = std(I);
cycorr = corr(cons,Y);
rycorr = corr(r,Y);

Simulation = [cvol Yvol Ivol cycorr rycorr]';
tab = table(Simulation,'RowNames',{'Consumption Volatility' 'Output Volatility' 'Investment Volatility' 'C-Y Correlation' 'r-Y Correlation'});
table2latex(tab,'q2vols.tex')