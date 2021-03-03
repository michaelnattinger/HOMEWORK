clear; close all; clc
alpha = 0.0001:0.0001:1;
kappa = alpha';
nume = 2 + 8*alpha.*kappa.^2- alpha.^2.*kappa.^2 + 6 *alpha.^2 .*kappa.^4 - alpha.^3.*kappa.^3 + 2*alpha.*kappa + 8*alpha.^2.*kappa.^3 + 6 * alpha.^3.*kappa.^5;
nume = nume(:);
check = nume<0;
disp(sum(check));
disp(min(nume));