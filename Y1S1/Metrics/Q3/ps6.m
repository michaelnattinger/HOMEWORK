clear; close all; clc
beta0 = 1;
pdelta = [0 1 1 1]';
philist = [0 0.8];
nlist = [40 70 100];
T = 4;
b_ols = zeros(length(nlist),length(philist),5);
b_fe = zeros(length(nlist),length(philist),4);
CIr = zeros(length(nlist),length(philist),2);
CIc = CIr;
scl = norminv(0.975);
for in = 1:length(nlist)
    n = nlist(in);
    for iphi = 1:length(philist)
        phi = philist(iphi);
        [X,Y] = sim_panel(beta0,pdelta,phi,n,T);
        t = repmat(diag([0 1 1 1]),n,1);
        t = t(:,2:4);
        b_ols(in,iphi,:) = ([ones(T*n,1) X(:) t])\Y(:);
        dX = X - mean(X);
        dY = Y - mean(Y);
        dt = t - mean(t);
        XX = ([dX(:) dt]);
        b_fe(in,iphi,:) = XX\dY(:);
        r = dY(:) - XX*reshape(b_fe(in,iphi,:),4,1);
        Vr = inv(XX'*XX)*(XX'*diag(r.^2)*XX)*inv(XX'*XX);
        clustsum = zeros(4,4);
        for i = 1:n
            clustsum = clustsum + XX(T*(i-1)+1:T*i,:)' *r(T*(i-1)+1:T*i,:) * r(T*(i-1)+1:T*i,:)' * XX(T*(i-1)+1:T*i,:);
        end
        Vc = inv(XX'*XX)*clustsum*inv(XX'*XX);
        CIr(in,iphi,1) = b_fe(in,iphi,1) - scl*sqrt(Vr(1));
        CIr(in,iphi,2) = b_fe(in,iphi,1) + scl*sqrt(Vr(1));
        CIc(in,iphi,1) = b_fe(in,iphi,1) - scl*sqrt(Vc(1));
        CIc(in,iphi,2) = b_fe(in,iphi,1) + scl*sqrt(Vc(1));
    end
end


%% simulation
nsim = 10000;
b_ols_sim = zeros(length(nlist),length(philist),5,nsim);
b_fe_sim = zeros(length(nlist),length(philist),4,nsim);
CRr = zeros(length(nlist),length(philist),nsim);
CRc = CRr;

for ns = 1:nsim
for in = 1:length(nlist)
    n = nlist(in);
    for iphi = 1:length(philist)
        phi = philist(iphi);
        [X,Y] = sim_panel(beta0,pdelta,phi,n,T);
        t = repmat(diag([0 1 1 1]),n,1);
        t = t(:,2:4);
        b_ols_sim(in,iphi,:,ns) = ([ones(T*n,1) X(:) t])\Y(:);
        dX = X - mean(X);
        dY = Y - mean(Y);
        dt = t - mean(t);
        XX = ([dX(:) dt]);
        b_fe_sim(in,iphi,:,ns) = XX\dY(:);
        r = dY(:) - XX*reshape(b_fe_sim(in,iphi,:,ns),4,1);
        Vr = inv(XX'*XX)*(XX'*diag(r.^2)*XX)*inv(XX'*XX);
        clustsum = zeros(4,4);
        for i = 1:n
            clustsum = clustsum + XX(T*(i-1)+1:T*i,:)' *r(T*(i-1)+1:T*i,:) * r(T*(i-1)+1:T*i,:)' * XX(T*(i-1)+1:T*i,:);
        end
        Vc = inv(XX'*XX)*clustsum*inv(XX'*XX);
        CIr_sim(1) = b_fe_sim(in,iphi,1,ns) - scl*sqrt(Vr(1));
        CIr_sim(2) = b_fe_sim(in,iphi,1,ns) + scl*sqrt(Vr(1));
        CIc_sim(1) = b_fe_sim(in,iphi,1,ns) - scl*sqrt(Vc(1));
        CIc_sim(2) = b_fe_sim(in,iphi,1,ns) + scl*sqrt(Vc(1));
        CRr(in,iphi,ns) = (beta0>CIr_sim(1))*(beta0<CIr_sim(2));
        CRc(in,iphi,ns) = (beta0>CIc_sim(1))*(beta0<CIc_sim(2));
    end
end
end
CRr = mean(CRr,3);
CRc = mean(CRc,3);
b_ols_sim = median(b_ols_sim,4);
b_fe_sim = median(b_fe_sim,4);

mat = zeros(6,4);
mat2 = zeros(6,9);

for i=1:length(nlist)
    for j=1:length(philist)
        row = (i-1)*2+j;
        mat(row,:) = [b_fe_sim(i,j,1) b_ols_sim(i,j,2) CRr(i,j) CRc(i,j)];
        mat2(row,:) = [b_fe(i,j,1) b_ols(i,j,2) reshape(b_fe(i,j,2:4),1,3) CIr(i,j,1) CIr(i,j,2) CIc(i,j,1) CIc(i,j,2)];
    end
end


tab = table(mat(:,1),mat(:,2),mat(:,3),mat(:,4), ...
    'RowNames',{'$\phi=0,n=40$' '$\phi=0.8,n=40$' '$\phi=0,n=70$' '$\phi=0.8,n=70$' '$\phi=0,n=100$' '$\phi=0.8,n=100$'}, ...
    'VariableNames',{'FE' 'OLS' 'EW' 'Cluster'});
table2latex(tab,'ps6.tex')
tab2 = table(mat2(:,1),mat2(:,2),mat2(:,3),mat2(:,4),mat2(:,5),mat2(:,6),mat2(:,7),mat2(:,8),mat2(:,9), ...
    'RowNames',{'$\phi=0,n=40$' '$\phi=0.8,n=40$' '$\phi=0,n=70$' '$\phi=0.8,n=70$' '$\phi=0,n=100$' '$\phi=0.8,n=100$'}, ...
    'VariableNames',{'FE' 'OLS' '$\delta_2$' '$\delta_3$' '$\delta_4$' 'EW LB' 'EW UB' 'C LB' 'C UB'});
table2latex(tab2,'ps6p2.tex')