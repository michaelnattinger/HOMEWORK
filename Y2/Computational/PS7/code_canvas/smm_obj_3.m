function obj = smm_obj_3(MT,W,e,P)
% Given parameter guesses in P,
% residuals e, weight matrx W,
% computes obj function
x = simul(P,e);
m1 = zeros(3,P.H);
for ii=1:P.H
    mx = mean(x(:,ii));
    m1(:,ii) = [mx (((x(:,ii) - mx))'*(x(:,ii) - mx))/P.T ((x(2:end,ii) - mx)'*(x(1:end-1,ii) - mx))/(P.T-1)]';
end
MTH = mean(m1,2);
obj = (MT-MTH)'*W*(MT-MTH);
end