function sad = calc_saddle(kg,kss,css,D,psigma,palpha,pbeta,pdelta, ns)
% Calculates the saddle path over a grid of capital levels using the
% shooting method: evaluates the ns step ahead value of K and chooses the
% point which results in capital levels as close 
ng = length(kg);
sad = 0*kg;
for i=1:ng
    C = 0:0.0001:8; % evaluate over a grid of points - dprop parallelizes
    K = kg(i)+0*C;  % the calculations so it is decently efficient
    for t=1:ns
        [K,C] = dprop(K,C,D,palpha,psigma,pbeta,pdelta);
    end
    [~,ind] = min(abs(K-kss)+abs(C-css));
    sad(i) = C(ind);
end
end