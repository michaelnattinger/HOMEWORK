%function lin = linearize(mod,ss)
%Linearize about the steady state
%Michael B Nattinger 9/9/2020
lin.yp = solve(mod.f,mod.yp); % solve for yp
lin.f = [];
ypn = fieldnames(lin.yp);
for i=1:length(ypn)
    lin.f = [lin.f; lin.yp.(ypn{i})]; %build vector for calculating jac
end
lin.jac = jacobian(lin.f,mod.y); % jacobian
lin.ssjac = lin.jac;
% eval w/ params
for i = 1:length(mod.pn)
eval(['lin.ssjac = subs(lin.ssjac, ' mod.pn{i} ',param.(mod.pn{i}));']); % plug in parameters
end
% eval w/ vars
for i=1:length(mod.yn)
lin.ssjac = subs(lin.ssjac,mod.y(i),ss.y.(mod.yn{i}));
end
% Convert to double
lin.ssjac = double(subs(lin.ssjac));
[lin.sseigvecs,lin.sseigvals] = eig(lin.ssjac);