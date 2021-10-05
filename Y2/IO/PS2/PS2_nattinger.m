clear; close all; clc
% IO Problem Set 2
% Nattinger
rng(999) % set seed
% This first part of the code is as given to us by the pset
J = 1000;
unif = rand(J,1);
unif2 = rand(J,1);
Num = ceil(10*unif);
collude = [(Num(1:500)>=8); zeros(500,1)];
F = 1;
c0 = 1;
c1 = 0.9;
b0 = 1;
b1 = 0;
z = 0;
h = 0;
nu = 0;
eta = 0;
a0 = 3;
a1 = 1;
%% Question 2
n_exp = 2;
for i_exp = 1:n_exp
    switch i_exp
        case 1 %c
            LernerCN = c1./Num;% <===== My code starts here, eqn is in problem 1 (g)
            LernerM = c1;
        case 2 %d
            LernerCN = (a0+nu-b0-eta)./(a0+nu+Num*(b0+eta));
            LernerM = ones(size(Num)).*(a0+nu-b0-eta)/(a0+nu+b0+eta);
    end
elasticity = 1/c1;
Herf = 1./Num;
n_exp = 2;
Lerner = collude.*LernerM + (1-collude).*LernerCN;
lnLernerobs = log(Lerner) + 0.1*(unif2 - 0.5);
lnHerf = log(Herf);
n_reg = 3;
regcoeff = NaN(2,n_reg); %prealloc
regse = NaN(2,n_reg); %prealloc
regFp = NaN(1,n_reg);
regF = NaN(1,n_reg);
for i_reg = 1:n_reg
    switch i_reg
        case 1
            XX = lnHerf;
            YY = lnLernerobs;% Constant added in regression function
        case 2
            XX = lnHerf(1:500);
            YY = lnLernerobs(1:500);
        case 3
            XX = lnHerf(501:end);
            YY = lnLernerobs(501:end);
    end
    mdl = fitlm(XX,YY);
    n = length(YY);
    regcoeff(:,i_reg) = mdl.Coefficients.Estimate;
    regse(:,i_reg) = mdl.Coefficients.SE;
    % perform F test
    % Null hyp:
    r0 = YY-XX;
    r = r0 - mean(r0); % residuals of restricted regression in this case
    RSS1 = r'*r;
    RSS2 = mdl.SSE;
    F = (RSS1 - RSS2)/((RSS2/(n-2))); %F statistic 
    p0 = fcdf(F,1,n-2);
    regF(i_reg) = F;
    regFp(i_reg) = 1-p0; % p value for F test
end
% make tables
out = summary_coef_se(regcoeff,regse);
vars = {'cons' ' ' 'log(H)' ' '}';
tab = table(vars,out(:,1),out(:,2),out(:,3));
table2latex_mid(tab,['regs_' num2str(i_exp)]) %output to tex table
out = [regF;regFp];
vars = {'F' 'p(f$>$F)'}';
tab = table(vars,out(:,1),out(:,2),out(:,3));
table2latex_mid(tab,['regs_' num2str(i_exp) '_2'])
end


%% Question 3
a0 = 5;
a1 = 1;
F = 1;
b0 = 1;
b1 = 0;
n_exp = 2;
coeffs = NaN(2,n_exp);
SEs = coeffs;
for i_exp = 1:n_exp
    switch i_exp
        case 1
            v = 2*rand(1000,1) - 1;
            eta = 0*v;
        case 2
            eta = 2*rand(1000,1) - 1;
            v = 0*eta;
    end
    % calculate N
    N = (1/sqrt(F*a1))*(a0 + v - b0 - eta)-1;
    Herf = 1./N;
    % calc L
    Lerner = (a0 + v - b0 - eta)./(a0 + v + N.*(b0 + eta));
    % run OLS
    logH = log(Herf);
    logL = log(Lerner);
    noise = unif(1000,1)/10 - 0.05;
    logLobs = logL+noise;
    mdl = fitlm(logH,logLobs);
    coeffs(:,i_exp) = mdl.Coefficients.Estimate;
    SEs(:,i_exp) = mdl.Coefficients.SE;
end

out = summary_coef_se(coeffs,SEs);
vars = {'cons' ' ' 'log(H)' ' '}';
tab = table(vars,out(:,1),out(:,2));
table2latex_mid(tab,'reg3')