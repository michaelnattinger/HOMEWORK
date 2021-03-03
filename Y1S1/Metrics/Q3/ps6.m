clear; close all; clc
beta0 = 1;
pdelta = [0 1 1 1]';
philist = [0 0.8];
nlist = [40 70 100];
T = 4;
b_ols = zeros(length(nlist),length(philist),3);
b_fe = zeros(length(nlist),length(philist),2);
CIr = zeros(length(nlist),length(philist),2);
CIc = CIr;
scl = norminv(0.975);
for in = 1:length(nlist)
    n = nlist(in);
    for iphi = 1:length(philist)
        phi = philist(iphi);
        [X,Y] = sim_panel(beta0,pdelta,phi,n,T);
        t = repmat((1:T)',1,n);
        b_ols(in,iphi,:) = ([ones(T*n,1) X(:) t(:)])\Y(:);
        dX = X - mean(X);
        dY = Y - mean(Y);
        dt = t - mean(t);
        XX = ([X(:) t(:)]);
        b_fe(in,iphi,:) = XX\Y(:);
        r = Y(:) - XX*reshape(b_fe(in,iphi,:),2,1);
        Vr = inv(XX'*XX)*(XX'*diag(r.^2)*XX)*inv(XX'*XX);
        clustsum = zeros(2,2);
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
nsim = 100;
b_ols_sim = zeros(length(nlist),length(philist),nsim);
b_fe_sim = zeros(length(nlist),length(philist),nsim);
CRr = zeros(length(nlist),length(philist),nsim);
CRc = CRr;

for ns = 1:nsim
for in = 1:length(nlist)
    n = nlist(in);
    for iphi = 1:length(philist)
        phi = philist(iphi);
        [X,Y] = sim_panel(beta0,pdelta,phi,n,T);
        t = repmat((1:T)',1,n);
        b_ols_sim(in,iphi,:,ns) = ([ones(T*n,1) X(:) t(:)])\Y(:);
        dX = X - mean(X);
        dY = Y - mean(Y);
        dt = t - mean(t);
        XX = ([X(:) t(:)]);
        b_fe_sim(in,iphi,:,ns) = XX\Y(:);
        r = Y(:) - XX*reshape(b_fe_sim(in,iphi,:,ns),2,1);
        Vr = inv(XX'*XX)*(XX'*diag(r.^2)*XX)*inv(XX'*XX);
        clustsum = zeros(2,2);
        for i = 1:n
            clustsum = clustsum + XX(T*(i-1)+1:T*i,:)' *r(T*(i-1)+1:T*i,:) * r(T*(i-1)+1:T*i,:)' * XX(T*(i-1)+1:T*i,:);
        end
        Vc = inv(XX'*XX)*clustsum*inv(XX'*XX);
        CIr_sim(1) = b_fe_sim(in,iphi,1,ns) - scl*sqrt(Vr(1));
        CIr_sim(2) = b_fe_sim(in,iphi,1,ns) + scl*sqrt(Vr(1));
        CIc_sim(1) = b_fe_sim(in,iphi,1,ns) - scl*sqrt(Vc(1));
        CIc_sim(2) = b_fe_sim(in,iphi,1,ns) + scl*sqrt(Vc(1));
        CRr(in,iphi,ns) = ;% what goes here?
    end
end
end