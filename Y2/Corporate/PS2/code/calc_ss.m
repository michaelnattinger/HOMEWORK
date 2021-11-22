clear; close all; clc
ptheta = 0.7;
pr = 0.04;
pdelta = 0.15;
ppsi0 = 0.01;
prhoz = 0.7;
psigmae = 0.1;
peta = 0.5;
prhoc = 0.2;
pgamma = 2;

k = ((pr + pdelta)/ptheta)^(1/(ptheta - 1));
c = k^ptheta - pdelta*k;
I = pdelta*k;
M = 1/(1+pr);