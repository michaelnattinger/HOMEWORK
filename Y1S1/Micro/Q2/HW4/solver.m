clear; close all; clc

syms p q a c
f1 = (1-p-q)*(-a-c)+a*p == (p+q)*(1-a-c)+(1-p-q);
f2 = (1-p-q)*(a-c) == - (p+q)*c;
sol = solve([f1 f2],[a c]);

q = exp(-10:-1);

p = 2*q;%0*q;
q = 0.5-q;

c = (p+q-1)./(8*p.^2 + 14*p.*q - 8*p + 6*q.^2 - 7*q + 2);
