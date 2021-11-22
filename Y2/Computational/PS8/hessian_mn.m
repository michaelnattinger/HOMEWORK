function H = hessian_mn(B,XX)
[T,nx] = size(XX);
H = zeros(nx,nx);
for tt=1:T
    H = H + Lam(XX(tt,:)*B)*(1-Lam(XX(tt,:)*B))*XX(tt,:)'*XX(tt,:);
end
H = -1*H;