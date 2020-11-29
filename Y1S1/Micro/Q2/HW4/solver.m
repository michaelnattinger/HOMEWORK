clear; close all; clc

syms p q a c
f1 = (1-p-q)*(1-a-c)+a*p == (p+q)*(1-a-c);
f2 = (1-p-q)*(a-c) == - (p+q)*c;
sol = solve([f1 f2],[a c]);