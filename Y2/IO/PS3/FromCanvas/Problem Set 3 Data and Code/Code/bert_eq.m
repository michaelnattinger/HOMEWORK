function f=bert_eq(p,expall_i,alpha_i, mc_hat1, own_dum)
% This function computes the FOC of a Nash-Bertrand equilibrium.
% When equal to zero the FOC is satisfied.
% inputs:  	p--price(over which we will preform a search);
%					expall_i--all the other components of demand, which will not be allowed to vary;
%					alpha_i --individual price coeff;
%					mc_hat1 -- marginal costs;
%					own_dum -- ownership structure

global ns

eg=expall_i.*exp(p*ones(1,ns).*(kron(alpha_i,ones(24,1))));
shar_i=eg./(ones(size(eg,1),1)*(1+sum(eg)));
f=omega1(p, shar_i,alpha_i',own_dum)*(p-mc_hat1)+(mean(shar_i')');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%9.1  Calculates the value of:  ownership(p)*(p-mc)-share(p).  In equilibrium
%we have:
%	ownership(p)*(p-mc)-share(p)=0
%	