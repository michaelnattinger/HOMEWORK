clear; close all; clc
reload = 0;
if reload
x1=readtable('AK91.csv');
y = x1.lwage; ed = x1.educ; n = length(y);
yob = zeros(n,9); sob = zeros(n,50); qob = zeros(n,3);
for t=1:n
    if x1.yob(t)>30&&x1.yob(t)<40
        yob(t,x1.yob(t) - 30) = 1;
    end
    if x1.sob(t)>0 && x1.sob(t)< 51
        sob(t,x1.sob(t)) = 1;
    end
    if x1.qob(t)>1 && x1.qob(t)<5
        qob(t,x1.qob(t)-1) = 1;
    end
end
save 'cleanAK91'
else
load 'cleanAK91'
end
ssob = sum(sob); ssobi = ssob>0;
sob = sob(:,ssobi); % remove columns with no data
x = [ed ones(n,1) yob sob]; z = [qob ones(n,1) yob sob];
% formulas
bhat2sls = (x'*z*inv(z'*z)*z'*x)\(x'*z*inv(z'*z)*z'*y);
e = y-x*bhat2sls; Qzz = z'*z/n; Qxz = x'*z/n;
Om = 0*Qzz;
for i=1:n
    Om = Om + z(i,:)'*z(i,:)*(e(i)^2);
end
Om = Om/n;
minv = inv(Qxz*inv(Qzz)*Qxz');
VhB = minv*(Qxz*inv(Qzz)*Om*inv(Qzz)*Qxz')*minv/n;
tab = table([bhat2sls(1,1); sqrt(VhB(1,1))],'VariableNames', ...
    {'Results'},'RowNames',{'Beta hat' 'SE'});
table2latex(tab,'ps3.tex')

