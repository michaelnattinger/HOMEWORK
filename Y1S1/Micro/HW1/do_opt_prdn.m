function [a,b,c,d] = do_opt_prdn(xs,ys)
% Does optimization to fit production function of y=Ax^b
% no intercept
A0 = 1; b0 = 1; C0 = 1; d0 = 1;
fun = @(x)diffs(x(1),x(2),x(3),x(4),ys,xs);
p = fminsearch(fun,[A0 b0 C0 d0]);
a = p(1); b=p(2); c=p(3); d=p(4);
end

function [objMin] = diffs(A,b,C,d,ys,xs)
objMin = sum((ys - (A*xs.^b + C*xs.^d)).^2);
end