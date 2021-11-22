clear; close all; clc

% read-in data, not very exciting:
[x,xt] = xlsread('mpd','Sheet1'); % This part just reads-in and preprocesses data
names = xt(1,:);

i_y = strcmp(names, 'i_close_first_year'); % this is just temporary. Overwritten later.
YY = x(:,i_y);
T = length(YY); 
xlist = {'score_0' 'rate_spread' 'i_large_loan' 'i_medium_loan' 'i_refinance'  'age_r'  'cltv' 'dti' 'cu' 'first_mort_r' 'i_FHA' 'i_open_year2' 'i_open_year3' 'i_open_year4' 'i_open_year5' 'score_0' 'score_1' 'score_2' 'i_open_0' 'i_open_1' 'i_open_2'}; 
nx = length(xlist);
XX = zeros(T,nx);
for i_x = 1:nx
    ii = strcmp(names,xlist{i_x});
    XX(:,i_x) = x(:,ii);
end

YY = 1-XX(:,end-2:end);
ZZ = XX(:,end-5:end-3);
XX = XX(:,end-6:end);
%XX = [ones(T,1) XX]; % will just hardcode alpha separate in likelihoods
%xlist = {'constant' xlist{:,:}};
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


%% likelihood computation
% initial parameters
ndraw = 100;
params.alpha0 = 0;
params.alpha1 = -1;
params.alpha2 = -1;
params.beta = zeros(nx,1);
params.gamma = 0.3;
params.rho = 0.5;
x0 = [params.alpha0 params.alpha1 params.alpha2 params.gamma params.rho params.beta']';

% pre-draw everything via halton draws (will map through inverse cdf)
p = haltonset(3,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
draws_3 = net(p,ndraw);
p = haltonset(2,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
draws_2 = net(p,ndraw);
p = haltonset(1,'Skip',1e3,'Leap',1e2);
p = scramble(p,'RR2');
draws_1 = net(p,ndraw);
% pre read-in quadrature nodes and weights
[nodes1, weights1] = nwspgr('KPU', 1, 20);
[nodes2, weights2] = nwspgr('KPU', 2, 20);

tic
GHK = LL_GHK(XX,YY,ZZ,params,draws_1,draws_2); 
toc
tic
accrej = LL_accrej(XX,YY,ZZ,params,draws_1,draws_2,draws_3); 
toc
tic
quad = LL_quad(XX,YY,ZZ,params,nodes1,weights1,nodes2,weights2);
toc

ops = optimset('Display','none');
obj = @(x)objfn(x,XX,YY,ZZ,params,draws_1,draws_2);
tic
[xopt,negLL] = fminunc(obj,x0,ops);
toc
params.alpha0 = xopt(1);
params.alpha1 = xopt(2);
params.alpha2 = xopt(3);
params.gamma = xopt(4);
params.rho = xopt(5);
params.beta = xopt(6:end);
GHK_opt = LL_GHK(XX,YY,ZZ,params,draws_1,draws_2); 
accrej_opt = LL_accrej(XX,YY,ZZ,params,draws_1,draws_2,draws_3); 
quad_opt = LL_quad(XX,YY,ZZ,params,nodes1,weights1,nodes2,weights2);

vars = {'$\alpha_0$' '$\alpha_1$' '$\alpha_2$' '$\gamma$' '$\rho$'};

for i_x = 1:nx
    vars = {vars{:,:} [xlist{i_x}]};
end
vars = vars';
tab = table(vars,xopt,'VariableNames',{' ' '$\hat{\theta}$'});
table2latex(tab,'mle.tex',4);

xx = [quad; GHK; accrej];
vars = {'Quadrature' 'GHK' 'Accept/Reject'}';
tab = table(vars,xx,'VariableNames',{'Method' 'Log Likelihood'});
table2latex(tab,'sll.tex',4);

%function negLL = objfn(x,XX,YY,ZZ,params,nodes1,weights1,nodes2,weights2) 
function negLL = objfn(x,XX,YY,ZZ,params,draws_1,draws_2)
% turns vector input x into parameters and computes (negative) log
% likelihood
params.alpha0 = x(1);
params.alpha1 = x(2);
params.alpha2 = x(3);
params.gamma = x(4);
params.rho = x(5);
params.beta = x(6:end); % due to issues in the optimization routine I switched to GHK
%negLL = -1*LL_quad(XX,YY,ZZ,params,nodes1,weights1,nodes2,weights2);
negLL = -1*LL_GHK(XX,YY,ZZ,params,draws_1,draws_2);
end
