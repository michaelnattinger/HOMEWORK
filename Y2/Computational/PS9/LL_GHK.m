function GHK = LL_GHK(XX,YY,ZZ,params,draws_1,draws_2)
[T,nx] = size(XX);
GHK = 0;
sig = 1/(1-params.rho)^2;
sig = sqrt(sig);
for i_t=1:T % for each point in time
if YY(i_t,1) 
    prob = normcdf(-params.alpha0 - XX(i_t,:)*params.beta - ZZ(i_t,1)*params.gamma);
    GHK = GHK + log(prob);
elseif YY(i_t,2)
    prob1 = normcdf(params.alpha0 + XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma);
    eps0 = sig*norminv(draws_1*prob1); % drawn from truncated distribution
    prob2 = normcdf(-params.alpha1 - XX(i_t,:)*params.beta - ZZ(i_t,2)*params.gamma - params.rho*eps0);
    GHK = GHK + log(mean(prob1.*prob2));
elseif YY(i_t,3)
    prob1 = normcdf(params.alpha0 + XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma);
    eps0 = sig*norminv(draws_2(:,1)*prob1); % drawn from truncated distribution
    prob2 = normcdf(params.alpha1 + XX(i_t,:)*params.beta + ZZ(i_t,2)*params.gamma + params.rho*eps0);
    eta1 = norminv(draws_2(:,2).*prob2);
    prob3 = normcdf(-params.alpha2 - XX(i_t,:)*params.beta - params.rho^2*eps0 - params.rho*eta1 - ZZ(i_t,3)*params.gamma);
    GHK = GHK + log(mean(prob1.*prob2.*prob3));
else 
    prob1 = normcdf(params.alpha0 + XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma);
    eps0 = sig*norminv(draws_2(:,1)*prob1); % drawn from truncated distribution
    prob2 = normcdf(params.alpha1 + XX(i_t,:)*params.beta + ZZ(i_t,2)*params.gamma + params.rho*eps0);
    eta1 = norminv(draws_2(:,2).*prob2);
    prob3 = normcdf(params.alpha2 + XX(i_t,:)*params.beta + params.rho^2*eps0 + params.rho*eta1 + ZZ(i_t,3)*params.gamma);
    GHK = GHK + log(mean(prob1.*prob2.*prob3));
end
%if isinf(abs(GHK)); GHK = -10e10; end
end
end