% Define model for symbolic toolbox
% Michael B Nattinger 9/9/2020
syms k k_p I I_p
syms pR pkst pdelta
% Define model equations
f1 = k_p == (1-pdelta)*k + I;
f2 = I_p == pR*(I + k-pkst)/(1-pdelta);
% Create mod with parameters, variables, functions
mod.p  = [pR pkst pdelta];
mod.y  = [k I];
mod.yp = [k_p I_p];
mod.f  = [f1 f2];
mod.pn = fieldnames(param);