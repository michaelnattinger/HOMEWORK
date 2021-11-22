function LL = LL_mn(B,XX,YY)
LL = 0;
[T,~] = size(XX);
for tt= 1:T
    Pr = Lam(XX(tt,:)*B);
    LL = LL + log(Pr^YY(tt,1) * (1-Pr)^(1-YY(tt,1)));
end