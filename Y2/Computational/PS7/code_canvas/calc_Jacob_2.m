function Jacob = calc_Jacob_2(x1,e,P,s)
% Calculates numerical approximation to Jacobian
m1 = zeros(2,P.H);
for ii=1:P.H
    mx = mean(x1(:,ii));
    m1(:,ii) = [((x1(:,ii) - mx)'*(x1(:,ii) - mx))/P.T  ((x1(2:end,ii) - mx)'*(x1(1:end-1,ii) - mx))/(P.T-1)]';
end
MTH = mean(m1,2);
P1 = P;
P1.prho = P1.prho - s;
x2 = simul(P1,e);
m1 = zeros(2,P1.H);
for ii=1:P1.H
    mx = mean(x2(:,ii));
    m1(:,ii) = [ ((x2(:,ii) - mx)'*(x2(:,ii) - mx))/P1.T ((x2(2:end,ii) - mx)'*(x2(1:end-1,ii) - mx))/(P1.T-1)]';
end
MTH2 = mean(m1,2);
part1 = (MTH - MTH2)/s;
P1 = P;
P1.psigma = P1.psigma - s;
x2 = simul(P1,e);
m1 = zeros(2,P1.H);
for ii=1:P1.H
    mx = mean(x2(:,ii));
m1(:,ii) = [ ((x2(:,ii) - mx)'*(x2(:,ii) - mx))/P1.T ((x2(2:end,ii) - mx)'*(x2(1:end-1,ii) - mx))/(P1.T-1)]';
end
MTH2 = mean(m1,2);
part2 = (MTH - MTH2)/s;
Jacob = [part1 part2];
end