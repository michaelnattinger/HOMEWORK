clear; close all; clc
syms pbeta pkappa pphi psigma

mat = [1/pbeta -pkappa/pbeta; ...
    pphi - (1/(pbeta*psigma)) 1+pkappa/(pbeta*psigma)];
mat = subs(mat,pbeta,1);
[Q,Lambda] = eig(mat);
iQ = inv(Q);