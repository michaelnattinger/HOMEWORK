function scr = score_mn(B,XX,YY)
[T,nx] = size(XX);
scr = zeros(nx,1);
for tt=1:T
    scr = scr + (YY(tt) - Lam(XX(tt,:)*B)).*XX(tt,:)';
end