clear; close all; clc
% Simulates store game
x = 0:1e-4:1;
c = 3;
left = @(x1,x2)(2*c*(x2^2 - x1^2)/9)+2*c*(x2 - x1)/9+c*(x2^2 - x1^2)^2/(18*(x2-x1));
right = @(x2,x1)(c*(x2^2 - x1^2)/9)+4*c*(x2 - x1)/9+c*(x2^2 - x1^2)^2/(18*(x2-x1));
nx = length(x);
maxiter = 1e3;
BR = 0*x; 
count = 0;
for t=1:nx %check all starting positions for the first-moving firm
x2 = t;
converge = 0;
x20 = 1;
x10 = 1;
iter = 1;
while (converge<1)&&(iter<maxiter)
    %x1profits
    x1pr = 0*x;
    for i=1:x2
        x1pr(i) = left(x(i),x(x2));
    end
    for i = x2:nx
        x1pr(i) = right(x(i),x(x2));
    end
    [~,x1] = max(x1pr); % firm 1 moves to max profits
    if iter==1; BR(t) = x(x1); end
    %x2profits
    x2pr = 0*x;
    for i=1:x1
        x2pr(i) = left(x(i),x(x1));
    end
    for i = x1:nx
        x2pr(i) = right(x(i),x(x1));
    end
    [~,x2] = max(x2pr);
    if (x20==x2)&&(x10==x1); converge = 1; end
    x20 = x2;
    x10 = x1;
    iter = iter+1;
end
if (x1<1&&x1>0)||(x2<1&&x2>0); count = count+1; end
end
if count>0
disp(['did not end up at corner solution ' num2str(count) ' times'])
else
    disp('always ended up at corner solution')
end

figure
plot(x,BR,'k')
hold on
plot(BR,x,'k')
hold off
set(gcf,'Color',[1 1 1])

