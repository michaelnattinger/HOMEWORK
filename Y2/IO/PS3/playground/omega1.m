function f=omega1(p,shar_i,alpha_i,own_dum)
%This function computes the matrix Omega defined in the solution to PS2
%inputs:		p--price(over which we will perform a search);
%				shar_i--individual purchase probabilities
%				alpha_i--individual price coeff
%				mc_hat1--marginal costs;
%				own_dum--ownership structure

global ns
n=24;
o1=(ones(n,1)*alpha_i').*shar_i*shar_i';
o2=diag(sum(((ones(n,1)*alpha_i').*shar_i)'));
omega=(o2-o1)/ns;
f=omega.*(own_dum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%9.1.1  Calculates the ownership share matrix x the (partial share)/(partial price)
%





