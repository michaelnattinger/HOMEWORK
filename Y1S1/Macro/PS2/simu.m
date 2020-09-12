function simul = simu(lin,T,y1)
%%% This function simulates y (in dev from ss) forward from an initial value
simul.ys = zeros(size(lin.ssjac,1),T);
simul.ys(:,1) = y1;
for tt=2:T
    simul.ys(:,tt) = lin.ssjac*simul.ys(:,tt-1);
end
