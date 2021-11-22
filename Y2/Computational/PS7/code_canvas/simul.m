function x = simul(P,e)
if nargin<2
    e = normrnd(0,1,P.T,P.H);
end
x(1,:) = P.px0 + e(1,:);
for tt=2:P.T
    x(tt,:) = P.prho * x(tt-1,:) + P.psigma*e(tt,:);
end
end