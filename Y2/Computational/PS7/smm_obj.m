function obj = smm_obj(MT,W,e,P)
% Given parameter guesses in P,
% residuals e, weight matrx W,
% computes obj function
x = simul(P,e);
m1 = zeros(2,P.H);
for ii=1:P.H
m1(:,ii) = [mean(x(:,ii)) ((x(:,ii) - mean(x(:,ii)))'*(x(:,ii) - mean(x(:,ii))))/P.T]';
end
MTH = mean(m1,2);
obj = (MT-MTH)'*W*(MT-MTH);
end