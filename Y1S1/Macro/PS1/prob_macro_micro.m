function [prob, iter] = prob_macro_micro(threshold,maxiter)
x=0;
iter=0;
while x<threshold && iter<=maxiter % maybe it's never cooler!
    iter=iter+1;
    x=rand;
end
prob=1/iter;

%% A function
% a) starts with function
% b) fill in the outputs in the square brackets
% c) give function a name which MUST be the same as the program file (in
%       this case prob_macro_micro
% d) define outputs in the parentheses
% e) make sure you define all the elements.
% f) NEVER change the inputs in the function itself
%       e.g. in the loop you could ad threshold=threshold-0.01, but never do
%       that
