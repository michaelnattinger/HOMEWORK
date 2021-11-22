clear; close all; clc
%% read in data
[x,xt] = xlsread('mpd','Sheet1'); % This part just reads-in and preprocesses data
names = xt(1,:);
i_y = strcmp(names, 'i_close_first_year');
YY = x(:,i_y);
T = length(YY);
xlist = {'i_large_loan' 'i_medium_loan' 'rate_spread' 'i_refinance'  'age_r'  'cltv' 'dti' 'cu' 'first_mort_r' 'score_0' 'score_1'   'i_FHA' 'i_open_year2' 'i_open_year3' 'i_open_year4' 'i_open_year5'}; % 'i_open_year2-i_open_year5' also, need to hardcode this!
nx = length(xlist);
XX = zeros(T,nx);
for i_x = 1:nx
    ii = strcmp(names,xlist{i_x});
    XX(:,i_x) = x(:,ii);
end
%XX = [ones(T,1) XX x(:,25)-x(:,28)];
XX = [ones(T,1) XX];
%xlswrite('data',[YY XX],'Sheet1','A1')
xlist = {'constant' xlist{:,:}};
[~,nx] = size(XX);
for i_x = 1:nx % fixes _ issue with strings and latex
    strng = xlist{i_x};
    nstrng = strng;
    ns = length(strng);
    for ii=ns:-1:2
        if strcmp(strng(ii),'_')
            nstrng = [nstrng(1:ii-1) '\' nstrng(ii:end)];
        end
    end
    xlist{i_x} = nstrng;
end
%% Actual codes are below
B0 = zeros(nx,1); 
B0(1) = -1;
LL0 = LL_mn(B0,XX,YY);
SCR0 = score_mn(B0,XX,YY);
HSS0 = hessian_mn(B0,XX);
s = 1e-6; % delta for numerical derivatives
num_f = num_f_mn(B0,XX,YY,s); 
num_s = num_s_mn(B0,XX,YY,s);
obj = @(x) -1*LL_mn(x,XX,YY); % objective function for optim routines
B0 = [-6 1 0.5 0.5 0 1 0 0 1 0.5 -0.5 0 1 1 1 0.5 0.5]'; % good guess for starting point
tol = 1e-6; % use same tol for all optimization routines
ops = optimset('MaxFunEvals',1e6,'TolX',tol,'TolFun',tol,'Display','none');
tune = 0.5; % for newton's method
tic % start timer
[B_newt,J_newt] = fmin_mn(XX,YY,B0,tol,tune);
t_newt = toc; % end timer
tic
[B_nm,J_nm] = fminsearch(obj,B0,ops); % Nelder-Mead Simplex method
t_nm = toc;
tic
[B_bfgs,J_bfgs] = fminunc(obj,B0,ops); % BFGS method
t_bfgs = toc;

%% Results
xnames = {'computation time' 'loglikelihood' xlist{:,:}}';
tab = table(xnames,[t_newt; J_newt; B_newt],[t_nm; -1*J_nm; B_nm],[t_bfgs; -1*J_bfgs; B_bfgs],'VariableNames',{' ' 'Newton' 'Nelder-Mead' 'BFGS'});
table2latex(tab,'results',5)

tab2 = table(SCR0,num_f);
table2latex_mid(tab2,'score');

str = 'tab3 = table(';
for i=1:nx
    str = [str 'HSS0(:,' num2str(i) '),'];
end
str = str(1:end-1);
str = [str ');'];
eval(str);
%tab3 = table(HSS0);
table2latex_mid(tab3,'hessian');
%tab4 = table(num_s);
str = 'tab4 = table(';
for i=1:nx
    str = [str 'num_s(:,' num2str(i) '),'];
end
str = str(1:end-1);
str = [str ');'];
eval(str);
table2latex_mid(tab4,'num_hess');
