clear; close all; clc

[data,sh] = xlsread('cps09mar.xlsx','Sheet1');
sh = sh(1,:);
n_y = 'earnings';
i_y = strcmp(n_y,sh);
i_h = strcmp('hours',sh);
i_w = strcmp('week',sh);
y = log(data(:,i_y)./(data(:,i_h).*data(:,i_w)));
t = length(y);
n_x = {'education' 'age'};
x=ones(t,1);
for i=1:length(n_x)
    x = [x data(:,strcmp(n_x{i},sh))];
end
x =[x (x(:,1+find(strcmp(n_x,'age'))) - x(:,1+find(strcmp(n_x,'education'))) -6)];
x(:,1+find(strcmp(n_x,'age'))) = x(:,end);
x(:,end) = x(:,end).^2/100;
n_x = {'education' 'experience' 'experience\^{}2/100' 'constant'};

marital = data(:,strcmp(sh,'marital'));
race = data(:,strcmp(sh,'race'));
female = data(:,strcmp(sh,'female'));
male = ~female;
single = marital==7; %is this right?
asian = race==4;
sam = logical(single.*asian.*male);
x =x(sam,:);
y = y(sam,:);
weirdobs = x(:,3)>=45;
x = x(~weirdobs,:);
y = y(~weirdobs,:);
% reorder to put constant last
x = [x(:,2:end) x(:,1)];
beta = (x'*x)\(x'*y);
disp(['\beta = '] )
disp(beta)
r = y-x*beta;
ybar = mean(y);
r2 = 1-(r'*r)/((y-ybar)'*(y-ybar)); % compute sum of squares as inner products
disp(['R^2 = ' num2str(r2)])
sse  = (r'*r);
disp(['sse = ' num2str(sse)])
% re-estimate
x2i = x(:,[2:4]);
x2f = x(:,1);
x2f = x2f-x2i*(x2i\x2f);
y2 = y - x2i*(x2i\y);
b2 = x2f\y2;
res2 = y2-x2f*b2;
r22 = 1-(res2'*res2)/((y2-mean(y2))'*(y2-mean(y2)));
sse2 = res2'*res2;
disp(['new beta is ' num2str(b2)])
disp(['R2 is ' num2str(r22)])
disp(['sse = ' num2str(sse2)])

a = sum(r);
b = sum(r.*x(:,strcmp('education',n_x)));
c = sum(r.*x(:,strcmp('experience',n_x)));
d = sum(r.*x(:,strcmp('education',n_x)).^2);
e = sum(r.*x(:,strcmp('experience',n_x)).^2);
f = sum(r.*(x*beta));
g = sum(r.*r);

% make tables for tex
table1 = table(beta,'RowNames',n_x);
results = [r2 sse]';
table2 = table(results,'RowNames',{'R\^{}2' 'SSE'});
reestimate = [b2 r22 sse2]';
table3 = table(reestimate,'RowNames',{'coefficient estimate' 'R\^{}2' 'SSE'});
sums = [a b c d e f g]';
table4 = table(sums,'RowNames',{'a' 'b' 'c' 'd' 'e' 'f' 'g'});

table2latex(table1,'table1.tex')
table2latex(table2,'table2.tex')
table2latex(table3,'table3.tex')
table2latex(table4,'table4.tex')