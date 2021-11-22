clear; close all; clc
[x,xt] = xlsread('FHP_Illistration','Sheet1');

c1 = 29;
c2 = 59;
x1 = x(1:c1,:);
x2 = x(c1+1:c2,:);
x3 = x(c2+1:end,:);
yy1 = x1(:,5);
xx1 = [ones(size(yy1,1),1) x1(:,[6 7])];
yy2 = x2(:,5);
xx2 = [ones(size(yy2,1),1) x2(:,[6 7])];
yy3 = x3(:,5);
xx3 = [ones(size(yy3,1),1) x3(:,[6 7])];

b = zeros(3);
n = zeros(3,1);
r2 = n;
for ii=1:3
    switch ii
        case 1
            xx = xx1;
            yy = yy1;
        case 2
            xx = xx2;
            yy = yy2;
        case 3
            xx = xx3;
            yy = yy3;
    end
    nnans = ~isnan(yy.*prod(xx,2));
    xx = xx(nnans,:);
    yy = yy(nnans,:);
    b(:,ii) = xx\yy;
    mdl = fitlm(xx(:,2:end),yy);
    r2(ii) = mdl.Rsquared.Ordinary;
    n(ii) = length(yy);
end
vars = {'intercept' 'Q_{t-1}' 'CF'}';
tab = table(vars,b(:,1),b(:,2),b(:,3),'VariableNames',{' ' '(1)' '(2)' '(3)'});
mat = [r2'; n'];
vars = {'R2' 'N'}';
tab2 = table(vars,mat(:,1),mat(:,2),mat(:,3),'VariableNames',{' ' '(1)' '(2)' '(3)'});
table2latex_mid(tab,'tab1')
table2latex_mid(tab2,'tab2')