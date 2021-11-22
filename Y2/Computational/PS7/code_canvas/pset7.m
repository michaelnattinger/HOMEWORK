mkdir('pings')
clear; close all; clc
rng(99999)
% Computational pset 7: SMM.
% Michael Nattinger 10/25/2021
%% Parameters (true)
P0.prho = 0.5;
P0.psigma = 1;
P0.px0 = 0;
P0.T = 200;
P0.H = 1;
P.T = P0.T;
P.px0 = 0;

%% Part 0: initial simulation
x0 = simul(P0);
MT = [mean(x0) ((x0 - mean(x0))'*(x0 - mean(x0))/P0.T)]';
P.H = 10;
e = normrnd(0,1,P.T,P.H);

%% Part 1: underidentification
W = eye(2);
ng = 100;
obj = zeros(ng,ng);
g_sig = linspace(0.8,1.2,ng);
g_rho = linspace(0.35,0.65,ng);
for i_sig = 1:ng
   P.psigma = g_sig(i_sig);
   for i_rho = 1:ng
      P.prho = g_rho(i_rho); 
       obj(i_sig,i_rho) = smm_obj(MT,W,e,P);
   end
end
figure
surf(g_sig,g_rho,obj)
set(gcf,'Color',[1 1 1])
title('Objective function (under-id)')
cd('pings')
saveas(gcf,'obj_un.png')
cd('..')

obfm = @(x) smm_fm(x,MT,W,e,P);
b1 = fminsearch(obfm,[0.5 1]);%fmincon(obfm,[0.5 1],[],[],[],[],[-10 0],[10 10]);%
P.prho = b1(1);
P.psigma = b1(2);
% given b1, compute W*
iT = 4; % lag length for NW Cov
x1 = simul(P,e);
SyTH = GammajTH(0,x1,P);
for j=1:iT
    gamm = GammajTH(j,x1,P);
    SyTH = SyTH + (1 - (j/(iT + 1)))*(gamm+gamm');
end
STH = (1+(1/P.H))*SyTH;
Wst = inv(STH);
obfm = @(x) smm_fm(x,MT,Wst,e,P);
[b2,J2] = fminsearch(obfm,[0.5 1]);
P.prho = b2(1);
P.psigma = b2(2);
s = 1e-12;
x1 = simul(P,e);
Jacob = calc_Jacob(x1,e,P,s);
Cov = inv(Jacob'*Wst*Jacob)/P.T;
SE = sqrt(diag(Cov))';
vars = {'rho' 'sigma'}';
tab = table(vars,b1',b2',SE','VariableNames',{' ' 'b1' 'b2' 'SE'});
jtab = table(Jacob(:,1),Jacob(:,2));
table2latex(tab,'tab_un');
table2latex_mid(jtab,'jtab_un')


Jtest = P.T*(P.H/(1+P.H))*J2;
disp(['J test: ' num2str(Jtest)])
disp('(unsurprisingly) this is not zero.')
%% Part 2: exact identification
mx0 = mean(x0);
MT = [((x0 - mx0)'*(x0 - mx0)/P0.T) ((x0(2:end,:) - mx0)'*(x0(1:end-1,:) - mx0))/(P0.T-1)]';
obj_2 = 0*obj;
for i_sig = 1:ng
   P.psigma = g_sig(i_sig);
   for i_rho = 1:ng
      P.prho = g_rho(i_rho); 
       obj_2(i_sig,i_rho) = smm_obj_2(MT,W,e,P);
   end
end
figure
surf(g_sig,g_rho,obj_2)
set(gcf,'Color',[1 1 1])
title('Objective function (exact-id)')
cd('pings')
saveas(gcf,'obj_ex.png')
cd('..')
W=eye(2);
obfm_2 = @(x) smm_fm_2(x,MT,W,e,P);
b1_2 = fminsearch(obfm_2,[0.5 1]); % much better
P.prho = b1_2(1);
P.psigma = b1_2(2);
% given b1, compute W*
iT = 4; % lag length for NW Cov
x1 = simul(P,e);
SyTH = GammajTH_2(0,x1,P);
for j=1:iT
    gamm = GammajTH_2(j,x1,P);
    SyTH = SyTH + (1 - (j/(iT + 1)))*(gamm+gamm');
end
STH = (1+(1/P.H))*SyTH;
Wst = inv(STH);
obfm_2 = @(x) smm_fm_2(x,MT,Wst,e,P);
[b2_2,J2_2] = fminsearch(obfm_2,[0.5 1]);
P.prho = b2_2(1);
P.psigma = b2_2(2);
s = 1e-12;
x1 = simul(P,e);
Jacob = calc_Jacob_2(x1,e,P,s);
Cov = inv(Jacob'*Wst*Jacob)/P.T;
SE2 = sqrt(diag(Cov))';
Jtest = P.T*(P.H/(1+P.H))*J2_2;

tab = table(vars,b1_2',b2_2',SE2','VariableNames',{' ' 'b1' 'b2' 'SE'});
jtab = table(Jacob(:,1),Jacob(:,2));
table2latex(tab,'tab_ex');
table2latex_mid(jtab,'jtab_ex')

%% Part 3: bootstrapping
mx0 = mean(x0);
MT = [mx0 ((x0 - mx0)'*(x0 - mx0)/P0.T) ((x0(2:end,:) - mx0)'*(x0(1:end-1,:) - mx0))/(P0.T-1)]';
W = eye(3);
obj_3 = 0*obj;
for i_sig = 1:ng
   P.psigma = g_sig(i_sig);
   for i_rho = 1:ng
      P.prho = g_rho(i_rho); 
       obj_3(i_sig,i_rho) = smm_obj_3(MT,W,e,P);
   end
end
figure
surf(g_sig,g_rho,obj_3)
set(gcf,'Color',[1 1 1])
title('Objective function (over-id)')
cd('pings')
saveas(gcf,'obj_ov.png')
cd('..')

ns = 5000;
W = eye(3);
b1_3 = zeros(ns,2);
b2_3 = b1_3;
for is = 1:ns
    rng(is) %new seed, I actually don't need to do this but might as well
    [x0] = simul(P0); % simulate the data
    mx = mean(x0);
    MT = [mean(x0) ((x0 - mx)'*(x0 - mx)/P0.T) ((x0(2:end,1) - mx)'*(x0(1:end-1,1) - mx)/(P0.T-1))]';
    e = normrnd(0,1,P.T,P.H);
    
    obfm_3 = @(x) smm_fm_3(x,MT,W,e,P);
    b1_3(is,:) = fminsearch(obfm_3,[0.5 1]);
    P.prho = b1_3(is,1);
    P.psigma = b1_3(is,2);
    x1 = simul(P,e);
    SyTH = GammajTH_3(0,x1,P);
    for j=1:iT
        gamm = GammajTH_3(j,x1,P);
        SyTH = SyTH + (1 - (j/(iT + 1)))*(gamm+gamm');
    end
    STH = (1+(1/P.H))*SyTH;
    Wst = inv(STH);
    obfm_3 = @(x) smm_fm_3(x,MT,Wst,e,P);
    [b2_3(is,:),J2_3] = fminsearch(obfm_3,[0.5 1]);
    if is<1.5 % first draw
        P.prho = b2_3(is,1);
        P.psigma = b2_3(is,2);
        s = 1e-12;
        x1 = simul(P,e);
        Jacob = calc_Jacob_3(x1,e,P,s);
        Cov = inv(Jacob'*Wst*Jacob)/P.T;
        SE3 = sqrt(diag(Cov))';
        Jtest = P.T*(P.H/(1+P.H))*J2_3;

        tab = table(vars,b1_3(is,:)',b2_3(is,:)',SE3','VariableNames',{' ' 'b1' 'b2' 'SE'});
        jtab = table(Jacob(:,1),Jacob(:,2));
        table2latex(tab,'tab_ov');
        table2latex_mid(jtab,'jtab_ov')
    end
end
b2_3mean = mean(b2_3);

figure
subplot(2,1,1)
histogram(b2_3(:,1),50)
title('\rho estimates')
subplot(2,1,2)
histogram(b2_3(:,2),50)
title('\sigma estimates')
set(gcf,'Color',[1 1 1])
suptitle(['Histogram of estimates from ' num2str(ns) ' bootstrap samples'])
cd('pings')
saveas(gcf,'boot.png')
cd('..')