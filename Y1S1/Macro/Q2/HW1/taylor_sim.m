function [tay_sim,y0] = taylor_sim(lin,ss,T,t0,y0)
% 1. solve for c_1,c_2
c_1 = 0; % Nonexplosive soln
c_2 = (y0(1) - ss.y.k)/(lin.sseigvecs(1,2)*lin.sseigvals(2,2)^(t0));
y0(2) = c_2*lin.sseigvecs(2,2)*lin.sseigvals(2,2)^(t0) + ss.y.I;
t = t0:T;
tay_sim = zeros(length(y0),length(t));
count = 1;
for tt=t
    tay_sim(:,count) = (lin.sseigvecs .* [c_1 c_2]) * diag(lin.sseigvals.^tt) + [ss.y.k ; ss.y.I];
    count = count+1;
end