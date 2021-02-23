%% parameter definitions and preallocation
clear; close all; clc
rng(999); % for reproducability
alpha0 = 1;%100; %scale down all coefficients by 100 to look nicer in table
delta0 = 1;%100;
beta0 = 0.01;%1;
T = [50 150 250];
rho1 = [0.7 0.9 0.95];
betas = zeros(4,3,3); % preallocation
CIl = zeros(4,3,3);
CIu = zeros(4,3,3);
scl = norminv(0.975); %~1.96 
%% Part (i): generate results once
for i=1:3 
    for j=1:3
        [Y,X] = gen_dat(T(i),rho1(j),alpha0,delta0,beta0);
        % calculate OLS coefficients, and h-r 95% CI
        tX = [ones(T(i),1) (1:T(i))' X(2:end) Y(1:end-1) ]; % rhs matrix
        betas(:,i,j) = tX\Y(2:end,:);
        r = Y(2:end,:) - tX*betas(:,i,j);
        Vhc0 = inv(tX'*tX)*(tX'*diag(r.^2)*tX)*inv(tX'*tX);
        SE = sqrt(diag(Vhc0));
        CIl(:,i,j) = betas(:,i,j)-scl*SE;
        CIu(:,i,j) = betas(:,i,j)+scl*SE;
    end
end

%% Part (2): repeat and calculate simulation results
nsim = 10000;
betasN = zeros(4,3,3,nsim);
coverN = zeros(4,3,3,nsim);
for n=1:nsim
for i=1:3
    for j=1:3
        [Y,X] = gen_dat(T(i),rho1(j),alpha0,delta0,beta0);
        % calculate OLS coefficients, and h-r 95% CI
        tX = [ones(T(i),1) (1:T(i))' X(2:end) Y(1:end-1) ];
        betasN(:,i,j,n) = tX\Y(2:end,:);
        r = Y(2:end,:) - tX*betasN(:,i,j,n);
        Vhc0 = inv(tX'*tX)*(tX'*diag(r.^2)*tX)*inv(tX'*tX);
        SE = sqrt(diag(Vhc0));
        Cl = betasN(:,i,j,n)-scl*SE;
        Cu = betasN(:,i,j,n)+scl*SE;
        if alpha0>Cl(1) && alpha0<Cu(1) % covered true value?
            coverN(1,i,j,n) = 1;
        end
        if beta0>Cl(2) && beta0<Cu(2)
            coverN(2,i,j,n) = 1;
        end
        if delta0>Cl(3) && delta0<Cu(3)
            coverN(3,i,j,n) = 1;
        end
        if rho1(j)>Cl(4) && rho1(j)<Cu(4)
            coverN(4,i,j,n) = 1;
        end
    end
end
end
PcoverN = mean(coverN,4); %average over our simulations
MbetasN = mean(betasN,4);
%% output tex tables
% construct table for single iteration results
mats = zeros(9,9);
for i = 1:3
    for j=1:3
        mats(3*(i-1)+j,[1 4 7]) = betas([1 3 4],i,j)';
        mats(3*(i-1)+j,1+[1 4 7]) = CIl([1 3 4],i,j)';
        mats(3*(i-1)+j,2+[1 4 7]) = CIu([1 3 4],i,j)';
    end
end
tabs = table(mats(:,1),mats(:,2),mats(:,3),mats(:,4),mats(:,5),mats(:,6),mats(:,7),mats(:,8),mats(:,9),...
    'RowNames',{'$(T=50,\rho_1=0.7)$' '$(T=50,\rho_1=0.9)$' '$(T=50,\rho_1=0.95)$' ...
    '$(T=150,\rho_1=0.7)$' '$(T=150,\rho_1=0.9)$' '$(T=150,\rho_1=0.95)$' ...
    '$(T=250,\rho_1=0.7)$' '$(T=250,\rho_1=0.9)$' '$(T=250,\rho_1=0.95)$'}, ...
    'VariableNames',{'$\hat{\alpha}_0$' '$\hat{\alpha}_0$ LB' '$\hat{\alpha}_0$ UB' ...
    '$\hat{\delta}_0$' '$\hat{\delta}_0$ LB' '$\hat{\delta}_0$ UB'  ...
    '$\hat{\rho}_1$' '$\hat{\rho}_1$ LB' '$\hat{\rho}_1$ UB'});
table2latex(tabs,'ps5s.tex')
% construct table for simulation results
mat = zeros(9,6);
for i = 1:3
    for j=1:3
        mat(3*(i-1)+j,[1 3 5]) = MbetasN([1 3 4],i,j)';
        mat(3*(i-1)+j,[2 4 6]) = PcoverN([1 3 4],i,j)'; 
    end
end
tab = table(mat(:,1),mat(:,2),mat(:,3),mat(:,4),mat(:,5),mat(:,6),...
    'RowNames',{'$(T=50,\rho_1=0.7)$' '$(T=50,\rho_1=0.9)$' '$(T=50,\rho_1=0.95)$' ...
    '$(T=150,\rho_1=0.7)$' '$(T=150,\rho_1=0.9)$' '$(T=150,\rho_1=0.95)$' ...
    '$(T=250,\rho_1=0.7)$' '$(T=250,\rho_1=0.9)$' '$(T=250,\rho_1=0.95)$'}, ...
    'VariableNames',{'$\hat{\alpha}_0$ Mean' '$\hat{\alpha}_0$ Coverage' ...
    '$\hat{\delta}_0$ Mean' '$\hat{\delta}_0$ Coverage' '$\hat{\rho}_1$ Mean' '$\hat{\rho}_1$ Coverage'});
table2latex(tab,'ps5.tex')