function [Y,X] = gen_dat(T,rho1,alpha0,delta0,beta0)
% generates time series data for HW5
V = normrnd(0,1,T,1);
U = normrnd(0,1,T,1);
Y = zeros(T+1,1);
X = Y;
Y(1) = normrnd(0,1,1,1);
X(1) = normrnd(0,1,1,1);
for i=1:T
    X(i+1) = X(i)*0.3 + V(i);
    Y(i+1) = alpha0 + i*beta0 + X(i+1)*delta0 + Y(i)*rho1 + U(i);
end