function obj = smm_obj_2(MT,W,e,P)
% Given parameter guesses in P,
% residuals e, weight matrx W,
% computes obj function
x = simul(P,e);
m1 = zeros(2,P.H);
for ii=1:P.H
    mx = mean(x(:,ii));
    m1(:,ii) = [(((x(:,ii) - mx))'*(x(:,ii) - mx))/P.T ((x(2:end,ii) - mx)'*(x(1:end-1,ii) - mx))/(P.T-1)]';
end
MTH = mean(m1,2);
obj = (MT-MTH)'*W*(MT-MTH);
end