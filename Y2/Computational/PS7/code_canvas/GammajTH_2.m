function Gamm = GammajTH_2(j,x1,P)
m1 = zeros(2,P.H);
for ii=1:P.H
    mx = mean(x1(:,ii));
    m1(:,ii) = [((x1(:,ii) - mx)'*(x1(:,ii) - mx))/P.T ((x1(2:end,ii) - mx)'*(x1(1:end-1,ii) - mx))/(P.T-1)]';
end
MTH = mean(m1,2);
Gamm = zeros(2);
for hh = 1:P.H
    xbar = mean(x1(:,hh));
    for tt = j+2:P.T
        Gamm = Gamm + (1/((P.T-1)*P.H))*(([(x1(tt,hh) -xbar)^2; (x1(tt,hh)-xbar)*(x1(tt-1,hh)-xbar)] -MTH)*([ (x1(tt-j,hh) -xbar)^2; (x1(tt-j,hh) -xbar)*(x1(tt-1-j,hh) -xbar) ]-MTH)');
    end
end
end
