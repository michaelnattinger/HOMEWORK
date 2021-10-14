function f = var_cov(theta2)
% This function computes the VCov matrix of the estimates

% Written by Aviv Nevo, May 1998.

global invA IV

load mvalold
load ps2
load gmmresid

N = size(x1,1);
Z = size(IV,2);
temp = jacob(mvalold,theta2);
a = [x1 temp]'*IV;
IVres = IV.*(gmmresid*ones(1,Z));
b = IVres'*IVres;

%f = gmmresid'*gmmresid/N*inv(a*invA*a');
f = inv(a*invA*a')*a*invA*b*invA*a'*inv(a*invA*a');

