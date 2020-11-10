function [mod,ss,lin,tay_sim,ex_sim,sh_sim] = calcmodss(param,do_sims,y0)
% Calculates model and ss as defined in : model.m script and steadystate.m
% script - Michael B. Nattinger, 9/9/20
model
steadystate
linearize
if do_sims
T = 20; t0 = 5;
[tay_sim,y0] = taylor_sim(lin,ss,T,t0,y0);
T = 16;
exact_sim
ex_sim = simsim;
sh_method_sim
else
sh_sim=[]; sh_c=[]; ex_sim=[]; simsim=[];
end
end

