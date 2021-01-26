function sad = calc_saddle(kg,kss,D,psigma,palpha,pbeta,pdelta, ns)
ng = length(kg);
sad = 0*kg;
for i=1:ng
    C = 0:0.0001:10;
    K = kg(i)+0*C;
    for t=1:ns
        [K,C] = dprop(K,C,D,palpha,psigma,pbeta,pdelta);
    end
    [~,ind] = min(abs(K-kss));
    sad(i) = C(ind);
end
end