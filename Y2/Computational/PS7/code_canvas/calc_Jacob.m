function Jacob = calc_Jacob(x1,e,P,s)
% Calculates numerical approximation to Jacobian
m1 = zeros(2,P.H);
for ii=1:P.H
m1(:,ii) = [mean(x1(:,ii)) ((x1(:,ii) - mean(x1(:,ii)))'*(x1(:,ii) - mean(x1(:,ii))))/P.T]';
end
MTH = mean(m1,2);
P1 = P;
P1.prho = P1.prho - s;
x2 = simul(P1,e);
m1 = zeros(2,P1.H);
for ii=1:P1.H
m1(:,ii) = [mean(x2(:,ii)) ((x2(:,ii) - mean(x2(:,ii)))'*(x2(:,ii) - mean(x2(:,ii))))/P1.T]';
end
MTH2 = mean(m1,2);
part1 = (MTH - MTH2)/s;
P1 = P;
P1.psigma = P1.psigma - s;
x2 = simul(P1,e);
m1 = zeros(2,P1.H);
for ii=1:P1.H
m1(:,ii) = [mean(x2(:,ii)) ((x2(:,ii) - mean(x2(:,ii)))'*(x2(:,ii) - mean(x2(:,ii))))/P1.T]';
end
MTH2 = mean(m1,2);
part2 = (MTH - MTH2)/s;
Jacob = [part1 part2];
end