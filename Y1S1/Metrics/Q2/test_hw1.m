clear; close all; clc
rng(99) % set seed - any seed is fine
n = 1000; % # obs for sim
e = randn(n,1); %true error
btrue = 1; % true beta
X = randn(n,1); % independent variables draws
Y = X*btrue + e; % dependent variable is X*b+e

bols = inv(X'*X)*(X'*Y);%X\Y; % OLS betahat
r = Y -X*bols; % OLS residuals

sum1 = sum(X.^2.*r); % quantity  we are asked about
sum2 = sum(X.^2.*Y) - sum(X.^3 * inv(X'*X)*(X'*Y)); %rewritten version of quantity - is exactly identical
disp(num2str(sum1));
disp(num2str(sum2));