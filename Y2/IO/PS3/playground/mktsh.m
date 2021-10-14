function f = mktsh(mval, expmu)
% This function computes the market share for each product

% Written by Aviv Nevo, May 1998.

global ns 
f = sum((ind_sh(mval,expmu))')/ns;
f = f';