function num_f = num_f_mn(B,XX,YY,s)
LL0 = LL_mn(B,XX,YY);
num_f = zeros(size(XX,2),1);
for ii=1:size(XX,2)
    del_b = zeros(size(XX,2),1);
    del_b(ii) = s;
    LL = LL_mn(B+del_b,XX,YY)-LL0;
    num_f(ii) = LL/s;
end