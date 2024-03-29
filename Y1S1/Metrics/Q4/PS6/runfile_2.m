clear; close all; clc
read = 0;
if read
    cd('..\PS4')
    [x,xt] = xlsread('cps09mar.xlsx','Sheet1');
    save 'data'
    cd('..\PS6')
else
    cd('..\PS4')
    load 'data'
    cd('..\PS6')
end

hisp = x(:,3);
women = x(:,2);
hiwo = logical(hisp.*women);
mlist = 1:9;
x = x(hiwo,:);
exp = x(:,1) - x(:,4) - 6;
edu = x(:,4);
n = size(x,1);
Y = log(x(:,5)./(x(:,6).*x(:,7)));
married = double(x(:,12)<4);
region = zeros(n,3);
for i=1:3
    region(x(:,10)==i,i) = 1;
end
AIC = 0*mlist;
BIC = AIC;
for mnum = mlist
    if mnum<4 % exp,sq
        X_exp = exp.^(1:2);
    elseif mnum<7 % exp - 4
        X_exp = exp.^(1:4);
    else % exp - six
        X_exp = exp.^(1:6);
    end
    if mod(mnum,3)==1 % college
        educ = double(edu>=16);
    elseif mod(mnum,3)==2 % spline
        educ = [edu zeros(n,1)];
        educ(edu>9,2) = edu(edu>9,1) - 9; 
    else % dummy
        educ = zeros(n,6);
        educ(edu == 12,1) = 1;
        educ(edu == 13,2) = 1;
        educ(edu == 14,3) = 1;
        educ(edu == 16,4) = 1;
        educ(edu == 18,5) = 1;
        educ(edu == 20,6) = 1;
    end
    X = [married region X_exp educ];
    mdl = fitlm(X,Y);
    AIC(mnum) = mdl.ModelCriterion.AIC;
    BIC(mnum) = mdl.ModelCriterion.BIC;
end
mat = [AIC;BIC];
tab = table(mat(:,1),mat(:,2),mat(:,3),mat(:,4),mat(:,5),mat(:,6),mat(:,7),mat(:,8),mat(:,9), ...
    'VariableNames',{'(1)' '(2)' '(3)' '(4)' '(5)' '(6)' '(7)' '(8)' '(9)'}, ...
    'RowNames',{'AIC' 'BIC'});
table2latex(tab,'table2',6)