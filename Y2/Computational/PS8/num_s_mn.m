function num_s = num_s_mn(B,XX,YY,s)
[~,N] = size(XX);
num_s = zeros(N);
for ii = 1:N
    delB = zeros(N,1);
    delB(ii) = s;
    num_s(:,ii) = num_f_mn(B,XX,YY,s);
    num_s(:,ii) = num_f_mn(B+delB,XX,YY,s) - num_s(:,ii);
    num_s(:,ii) = num_s(:,ii)./s;
end