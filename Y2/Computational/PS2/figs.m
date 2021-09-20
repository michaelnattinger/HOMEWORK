mkdir('pings')
clear; close all; clc

dat = dlmread('pfs_.dat');
cd('pings')
[T,N] = size(dat);
x = linspace(-2,5,T);
apr_e = dat(:,1);
apr_u = dat(:,2);
pmf_e = dat(:,3);
pmf_u = dat(:,4);
v_e   = dat(:,5);
v_u   = dat(:,6);
figure
plot(x,apr_e,'b-')
hold on
plot(x,apr_u,'r--')
plot(x,x,'k:')
[~,i_zero] = min(abs(x'-apr_e));
plot([x(i_zero) x(i_zero)],[-2 5],'m-.')
hold off
legend('a prime (e)','a prime (u)','45 degree line','a hat')
set(gcf,'Color',[1 1 1])
xlabel('Asset holdings today')
ylabel('Asset holdings tomorrow')
saveas(gcf,'aprime.png')
rbound = x(i_zero);
wealth_e = x+1;
wealth_u = x+0.5;
wealth = sort(unique([wealth_e wealth_u]));
wealth_mass_e = 0*wealth;
wealth_mass_u = 0*wealth;
for i=1:T
    wealth_mass_e(wealth_e(i)==wealth) = wealth_mass_e(wealth_e(i)==wealth) + pmf_e(i);
    wealth_mass_u(wealth_u(i)==wealth) = wealth_mass_u(wealth_u(i)==wealth) + pmf_u(i);
end
wealth_mass = wealth_mass_e + wealth_mass_u;
figure
plot(wealth,wealth_mass_e,'k')
hold on
plot(wealth,wealth_mass_u,'r')
hold off
legend('Population mass (e)','Population mass (u)')
set(gcf,'Color',[1 1 1])
ylabel('Probability mass')
xlabel('Wealth (a + s)')
saveas(gcf,'wealth.png')
weight_wealth = wealth_mass.*wealth;
lorenz = cumsum(weight_wealth)./sum(weight_wealth');
cs_wm = cumsum(wealth_mass);
figure
plot(cs_wm,lorenz,'b')
hold on
plot(cs_wm,0*cumsum(wealth_mass),'k');
plot(cs_wm,cumsum(wealth_mass),'k:')
hold off
set(gcf,'Color',[1 1 1])
title('Lorenz Curve')
saveas(gcf,'lorenz.png')
% calc_gini: integrate across cumsum(wealth_mass) calcing diffs between
% cumsum(wealth_mass) and lorenz
diff = cs_wm - lorenz;
% integrating here
isum = 0;
for i = 2:length(diff)
    isum = isum +  diff(i) *(cs_wm(i) -cs_wm(i-1)); 
end
gini = isum/(1/2);

% calculate welfare
mass_e = sum(pmf_e);
mass_u = sum(pmf_u);
c_all = mass_e + 0.5*mass_u;
palpha = 1.5;
u_flow = (c_all^(1-palpha) - 1)/(1-palpha);
pbeta = 0.9932;
W_FB = u_flow./(1-pbeta);
lambda_e = 0*v_e;
lambda_u = 0*v_e;
WG = 0;
WINC = 0;
count = 0;
for i = 1:length(v_e)
    lambda_e(i) = ((W_FB + 1/((1-palpha)*(1-pbeta)))/(v_e(i)+1/((1-palpha)*(1-pbeta))))^(1/(1-palpha)) - 1;
    lambda_u(i) = ((W_FB + 1/((1-palpha)*(1-pbeta)))/(v_u(i)+1/((1-palpha)*(1-pbeta))))^(1/(1-palpha)) - 1;
    WG = WG + lambda_e(i)*pmf_e(i);
    WG = WG + lambda_u(i)*pmf_u(i); 
    % calc WINC
    WINC = WINC + v_e(i)*pmf_e(i);
    WINC = WINC + v_u(i)*pmf_u(i);
    % calc frac_pop
    count = count + pmf_e(i)*(lambda_e(i)>0);
    count = count + pmf_u(i)*(lambda_u(i)>0);
end
figure
plot(x,lambda_e,'k')
hold on
plot(x,lambda_u,'b')
hold off
title('\lambda vs asset holdings')
set(gcf,'Color',[1 1 1])
legend('\lambda (e)','\lambda (u)')
xlabel('Asset holdings today')
ylabel('\lambda')
saveas(gcf,'lambda.png')
figure
plot(x,v_e,'k')
hold on
plot(x,v_u,'b')
hold off
title('Value functions')
set(gcf,'Color',[1 1 1])
cd('..')

