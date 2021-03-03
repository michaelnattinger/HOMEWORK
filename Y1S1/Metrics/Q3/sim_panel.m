function [X,Y] = sim_panel(beta0,pdelta,phi,n,T)
% generates panel data
X = zeros(T,n);
alpha = normrnd(0,1,1,n);
ind = alpha>0.6;
X(3:4,ind) = 1;
epsilon = zeros(T+1,n);
epsilon(1,:) = normrnd(0,1,1,n);
u = normrnd(0,1,T,n);
for tt = 1:T
    epsilon(tt+1,:) = phi*epsilon(tt,:) + u(tt,:);
end
Y = X*beta0 + repmat(pdelta,1,n) + repmat(alpha,T,1) + epsilon(2:end,:);