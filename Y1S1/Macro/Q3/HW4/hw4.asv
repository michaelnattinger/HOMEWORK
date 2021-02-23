clear; close all; clc
rng(999); % for reproducability
W = 1; % set wage to 1
eps = 1e-9;
prho = 1+eps; 
ptheta = 5; 
Nk = 20; 
K = 100000;
A = exp(normrnd(0,1,K,Nk)); % firm-level productivity draws
s0 = (1/Nk)+ 0*A; % initial guess for shares
diff = 999; 
tol = 1e-6; 
maxiter = 1e10;
iter = 1;
s = s0;
tune = 0.6; % how far are we moving s? increase as eps decreases
while (diff>tol) && (iter<maxiter) % repeat until converged
    if mod(iter,100) == 0; disp(['iter ' num2str(iter) ', diff = ' num2str(diff)]); end
    % calculate pik | s
    pik = (W./A).*(1 - (1./((1-ptheta) + s .*(ptheta - prho ))));
    % calculate pk | pik
    pk = (sum(pik.^(1-ptheta),2)).^(1/(1-ptheta));
    % calculate s |  pk,pik
    s0 = s;
    s = (pik./pk).^(1-ptheta);
    diff = sum(sum(abs(s-s0),2)); % check for convergence
    s = s0*tune + s*(1-tune); % update new guess for s: convex combination of old s0 and s|pik,pk
    iter = iter+1;
end
% Solved.
p = ((1/K)*sum(pk.^(1-prho))).^(1/(1-prho));
C = W/p; % c here is the real wage
pikspp = W./A; % best allocation
pkspp = (sum(pikspp.^(1-ptheta),2)).^(1/(1-ptheta));
pspp = ((1/K)*sum(pkspp.^(1-prho))).^(1/(1-prho));
Cspp = W/pspp;
disp(['C = ' num2str(C)])
disp(['Cspp = ' num2str(Cspp)])