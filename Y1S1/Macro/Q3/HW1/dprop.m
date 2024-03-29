function [Kp,Cp] = dprop(K,C,D,palpha,psigma,pbeta,pdelta)
% propogates the difference equations by one period given D, the
% parameterizations, and starting positions of K, C
% Built to handle parallelization via propogating K,C as vectors or
% matrices - this comes in handy for the shooting method
Kp = (1-pdelta)*K + K.^palpha - C - D;
Cp = C.*(pbeta * (palpha*K.^(palpha - 1) + 1 - pdelta)).^(1/psigma);
end