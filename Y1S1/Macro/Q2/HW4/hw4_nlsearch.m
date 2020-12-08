addpath('C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q2')
clear; close all; clc
% Finds capital demand guess that clears the market.
% Michael Nattinger, 2020
recalc1 = 0;
recalc2 = 0;
if recalc1>0
    Q = [0.85 0.15; 0.05 0.95];
    tol=1e-5;       % maybe jack these up
    maxiter=1e6;
    pbeta = 0.95; % parameter values
    pgamma = 3;
    palpha = 0.36;
    pdelta = 0.08;
    l = [0.7 1.1];
    abar = 50; % upper bound for asset grid
    agr = 5e-2;
    assets = 0:agr:abar; % asset grid
    na = length(assets);
    %V0 = ones(na,2); % value function initial guess
    V0 = [log(assets'+1) log(assets'+1)];
    k0 = 5.0162;%5.5;
    ktol = 0.01;
    kmiter = 1000;
    tic
    Kd = find_Kd_slowiter(0.99,k0,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter,ktol,kmiter);
    toc
    % calculate results
    [Ks,dist,V,aprime] = find_Ks(Kd,assets,Q,V0,pbeta,pgamma,pdelta,palpha,l,tol,maxiter);
    save results1
else
    load results1
end
if recalc2>0
    tol=1e-4;       % maybe jack these up
    maxiter=1e6;
    abar2 = 50;
    assets2 = -2:agr:abar2;
    na2 = length(assets2);
    V02 = [log(assets2'+3) log(assets2'+3)];
    k02 = 4.97;
    %Kd2 = find_Kd(k02,assets2,Q,V02,pbeta,pgamma,pdelta,palpha,l,tol,maxiter);
    % calculate results
    ktol = 0.01;
    kmiter = 1000;
    [Kd2,diff2,iter] = find_Kd_slowiter(0.99,k02,assets2,Q,V02,pbeta,pgamma,pdelta,palpha,l,tol,maxiter,ktol,kmiter);
    tic
    [Ks2,dist2,V2,aprime2] = find_Ks(k02,assets2,Q,V02,pbeta,pgamma,pdelta,palpha,l,tol,maxiter);
    toc
    save results
else
    load results
end

w = (1-palpha)*(Kd)^(palpha);
r = palpha*(Kd)^(palpha - 1);
Ilev = w*l.*ones(size(V)) + r*repmat(assets',1,2);
Clev  = Ilev - aprime + (1-pdelta)*repmat(assets',1,2);
% calculate cdf
Cvec = sort(unique(Clev));
Ivec = sort(unique(Ilev));
Cden = 0*Cvec;
Iden = 0*Ivec;
for nl = 1:2
   for ai = 1:na
        Cden(Cvec==Clev(ai,nl)) = Cden(Cvec==Clev(ai,nl))+dist(ai,nl);
        Iden(Ivec==Ilev(ai,nl)) = Iden(Ivec==Ilev(ai,nl))+dist(ai,nl);
   end
end
Ccdf = cumsum(Cden);
Icdf = cumsum(Iden);

w2 = (1-palpha)*(Kd2)^(palpha);
r2 = palpha*(Kd2)^(palpha - 1);
Ilev2 = w2*l.*ones(size(V2)) + r2*repmat(assets2',1,2);
Clev2  = Ilev2 - aprime2 + (1-pdelta)*repmat(assets2',1,2);
% calculate cdf
Cvec2 = sort(unique(Clev2));
Ivec2 = sort(unique(Ilev2));
Cden2 = 0*Cvec2;
Iden2 = 0*Ivec2;
for nl = 1:2
   for ai = 1:na2
        Cden2(Cvec2==Clev2(ai,nl)) = Cden2(Cvec2==Clev2(ai,nl))+dist2(ai,nl);
        Iden2(Ivec2==Ilev2(ai,nl)) = Iden2(Ivec2==Ilev2(ai,nl))+dist2(ai,nl);
   end
end
Ccdf2 = cumsum(Cden2);
Icdf2 = cumsum(Iden2);

figure
subplot(2,1,1)
plot(Cvec,Ccdf,'k-')
title('Consumption distribution')
subplot(2,1,2)
plot(Ivec,Icdf,'k-')
title('Income distribution')
set(gcf,'Color',[1 1 1])
suptitle('Consumption and income distributions')
cd('pings')
saveas(gcf,'dis.png')
cd('..')

figure
subplot(2,1,1)
plot(Cvec,Ccdf,'k-')
hold on
plot(Cvec2,Ccdf2,'r-')
hold off
title('Consumption distribution')
subplot(2,1,2)
plot(Ivec,Icdf,'k-')
hold on
plot(Ivec2,Icdf2,'r-')
hold off
title('Income distribution')
set(gcf,'Color',[1 1 1])
legend('debt limit 0','debt limit -2','Location','SouthEast')
suptitle('Consumption and income distributions')
cd('pings')
saveas(gcf,'comp.png')
cd('..')

Equilibrium = [Kd w r]';
DebtEquilibrium = [Kd2 w2 r2]';

tab1 = table(Equilibrium,'RowNames',{'Capital' 'Wage' 'Interest rate'});
tab2 = table(Equilibrium,DebtEquilibrium,'RowNames',{'Capital' 'Wage' 'Interest rate'});
table2latex(tab1,'tab1.tex')
table2latex(tab2,'tab2.tex')