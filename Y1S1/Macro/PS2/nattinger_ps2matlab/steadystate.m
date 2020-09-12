% Evaluate steady state values for model
% Michael B Nattinger, 9/9/2020
ss.f = subs(mod.f,mod.yp,mod.y); % steady state: yp = y = yss
for i = 1:length(mod.pn)
eval(['ss.f = subs(ss.f, ' mod.pn{i} ',param.(mod.pn{i}));']); % plug in parameters
end
ss.y = solve(ss.f,mod.y); % solve for y
mod.yn = fieldnames(ss.y);
for i = 1:length(mod.yn)
ss.y.(mod.yn{i}) = double(subs(ss.y.(mod.yn{i}))); % store ss value as a double
end