function accrej = LL_accrej(XX,YY,ZZ,params,draws_1,draws_2,draws_3)
[T,nx] = size(XX);
accrej = 0;
sig = 1/(1-params.rho)^2;
sig = sqrt(sig);
eps01d = sig*norminv(draws_1); % I predrew these from a halton sequence
eps02d = sig*norminv(draws_2(:,1));
eps03d = sig*norminv(draws_3(:,1));

eta12d = norminv(draws_2(:,2));
eta13d = norminv(draws_3(:,2));
eta23d = norminv(draws_3(:,3));

for i_t=1:T % for each point in time
if YY(i_t,1) 
    eps0 = eps01d;
    ind = eps0<(-params.alpha0 - XX(i_t,:)*params.beta - ZZ(i_t,1)*params.gamma);
    accrej = accrej + log(mean(ind));
elseif YY(i_t,2)
    eps0 = eps02d;
    eta1 = eta12d;
    ind1 = eps0<(params.alpha0 + XX(i_t,:)*params.beta  + ZZ(i_t,1)*params.gamma);
    ind2 = eta1<(-params.alpha1 - XX(i_t,:)*params.beta - params.rho*eps0 - ZZ(i_t,2)*params.gamma);
    accrej = accrej+log(mean(ind1.*ind2));
elseif YY(i_t,3)
    eps0 = eps03d;
    eta1 = eta13d;
    eta2 = eta23d;
    ind1 = eps0<(params.alpha0 + XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma);
    ind2 = eta1<(params.alpha1 + XX(i_t,:)*params.beta + params.rho*eps0 + ZZ(i_t,2)*params.gamma);
    ind3 = eta2<(-params.alpha2 - XX(i_t,:)*params.beta - params.rho^2*eps0 - params.rho*eta1 - ZZ(i_t,3)*params.gamma);
    accrej = accrej+log(mean(ind1.*ind2.*ind3));
else 
    eps0 = eps03d;
    eta1 = eta13d;
    eta2 = eta23d;
    ind1 = eps0<(params.alpha0 + XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma);
    ind2 = eta1<(params.alpha1 + XX(i_t,:)*params.beta + params.rho*eps0 + ZZ(i_t,2)*params.gamma);
    ind3 = eta2>(-params.alpha2 - XX(i_t,:)*params.beta - params.rho^2*eps0 - params.rho*eta1 - ZZ(i_t,3)*params.gamma);
    accrej = accrej+log(mean(ind1.*ind2.*ind3));
end
if isinf(abs(accrej)); accrej = -10e10; end % handle negative infinities, which will happen
end
end