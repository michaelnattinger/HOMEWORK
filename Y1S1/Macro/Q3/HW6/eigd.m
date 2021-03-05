clear; close all; clc
syms pbeta pkappa pphi psigma

mat = [1/pbeta -pkappa/pbeta; ...
    pphi - (1/(pbeta*psigma)) 1+pkappa/(pbeta*psigma)];
mat = subs(mat,pbeta,1);
[Q,Lambda] = eig(mat);
iQ = inv(Q);
C = [0; -1/psigma];
C = iQ*C;
il1inv = C(1)*Lambda(1)^(-1);
il2inv = C(2)*Lambda(4)^(-1);