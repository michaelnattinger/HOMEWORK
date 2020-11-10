function simres = dosims(mod,ss,lin,mod2,ss2,lin2)
% Evaluates simulation computations required for pset 1
% This is specifically for the pset so it is less automated than I
% typically want my code to be.
% Michael B Nattinger 9/11/2020
t0 = 5;
T = 20;
simres.sim1 = zeros(2,T);

for i=2:t0-1
simres.sim1(:,i) = lin.ssjac*simres.sim1(:,i-1);
end
% need to do rest of problem here
for j=1:length(mod.yn) %At t0, z is shocked. K and C values at t0 were 
                       %already determined at t0-1 so at t0 we are at the
                       %old steady state, but now the system has a new
                       %steady state so our deviation from the new ss is
                       %the difference between the steady states.
simres.sim1(j,t0) = ss.y.(mod.yn{j})- ss2.y.(mod.yn{j});
end

for i=t0+1:T % naive way of calculating, using jac at ss. Will change after
             % section this afternoon.
simres.sim1(:,i) = lin2.ssjac*simres.sim1(:,i-1); 
end

for j=1:length(mod.yn)
simres.sim1(j,1:t0-1) = simres.sim1(j,1:t0-1)+ss.y.(mod.yn{j}); % adding SS back in
simres.sim1(j,t0:end) = simres.sim1(j,t0:end)+ss2.y.(mod.yn{j});
end
end