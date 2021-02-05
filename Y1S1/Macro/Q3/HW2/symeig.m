clear; close all; clc

syms pbeta pdelta pphi ppsi palpha pgamma ptheta plambda

A = [pbeta^(-1), (-1)*pdelta*(pphi - 1); ...
    ppsi*(palpha - 1)/(1+pbeta * pgamma * (pphi - 1)), ...
    ptheta/(1+pbeta * pgamma * (pphi - 1))];
%[a,b] = eig(A);
ALI = A - eye(2)* plambda;
dALI = det(ALI);
sol = solve(dALI==0,plambda);
