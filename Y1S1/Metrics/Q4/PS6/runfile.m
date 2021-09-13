clear; close all; clc
exload = 0;
if exload
[x,names] = xlsread('CHJ2004','Sheet1');
save data
else
load data
end
tinkind = x(:,18)./1000;
tinkind(tinkind<0,:) = 0;
income = x(:,1)./1000;
Dincome = (income-1).*(income>1);
X1 = [income Dincome];
mdl1 = fitlm(X1,tinkind);
prop = mean(tinkind==0);
ind2 = (tinkind>0);
Y2 = tinkind(ind2,:);
X2 = X1(ind2,:);
mdl2 = fitlm(X2,Y2);
[mdl3,~,cov3] = TOBIT(tinkind,X1,0,inf,1);
[mdl4,cov4] = CLAD(tinkind,X1,0,inf,1);
tab=table(mdl1.Coefficients.Estimate,mdl2.Coefficients.Estimate,mdl3,mdl4, ...
    'VariableNames',{'A' 'C' 'D' 'E'},'RowNames',{'cons' 'inc' 'Dinc'});
table2latex(tab,'table1',5)