% Define model for symbolic toolbox
% Michael B Nattinger 9/9/2020
syms k k_p c c_p
syms pbeta palpha pz pdelta
% Define model equations
f1 = k_p == pz*k^(palpha) + (1-pdelta)*k - c;
f2 = pbeta/c_p == (c)^(-1)*(1-pdelta + palpha*pz*k_p^(palpha - 1))^(-1);
% Create mod with parameters, variables, functions
mod.p  = [pbeta palpha pz pdelta];
mod.y  = [k c];
mod.yp = [k_p c_p];
mod.f  = [f1 f2];
mod.pn = fieldnames(param);