mkdir('pings')
clear; close all; clc

del = 0.01; % width of cgrid
cg = 0:del:1-del; % cgrid
rgrid = linspace(0,1,1e6)';
nc = length(cg);
r1 = 0*cg; % prealloc
r2 = r1;
p1 = r1;
p2 = r1;
for i=1:nc
    [p1(i),r1(i)] = min(pr1(rgrid,cg(i)));
    [p2(i),r2(i)] = min(pr2(rgrid,cg(i)));
end
% options = optimoptions('fmincon','Display','off','FunctionTolerance',1e-12,'StepTolerance',1e-12,'MaxFunctionEvaluations',2000);
% for i = 1:nc % loop through cgrid and optimize at each point
%     obj1 = @(r) pr1(r,cg(i));
%     [r1(i),p1(i)] = fmincon(obj1,0.5,[1;-1],[1;0],[],[],[],[],[],options);
%     obj2 = @(r) pr2(r,cg(i));
%     [r2(i),p2(i)] = fmincon(obj2,0.5,[1;-1],[1;0],[],[],[],[],[],options);
% end

% fmincon fails too much - switch to grid search

p1 = -1*p1; % convert from negative of profit to profit
p2 = -1*p2;
r1 = rgrid(r1);
r2 = rgrid(r2);
figure
plot(cg,r1,'b--')
hold on
plot(cg,r2,'r-.')
plot(cg,cg,'k-')
hold off
title('(Expected) Profit Maximizing r')
ylabel('Optimal r')
xlabel('c')
legend('First-price','Second-price','c','Location','SouthEast')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'argmax.png')
cd('..')
figure
plot(cg,p1,'b--')
hold on
plot(cg,p2,'r-.')
hold off
title('Expected Profit')
legend('First-price','Second-price','Location','NorthEast')
ylabel('Expected profit')
xlabel('c')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'profit.png')
cd('..')

function f = pr1(r,c) % negative profit so optimum is the minimum - to use with fmincon
f = -1*(( (1/3).*(1-r.^3) - c ).*(r.^2-r.^3) + ( (1/6).*(1-r).^2.*(2.*r.^3 + r.^2 + 2*r +1) - c ).*( (1/3)*(1-r).^2.*(2*r+1)) + (r.^2 -(4/3)*r.^3 +(1/3) - c ).*(r-r.^3) + ((1/12)*(1-r).^2.*(7*r.^2+2*r+3)  - c ).*(1/3).*( r.^3 - 3*r + 2 )) ;
end
function f = pr2(r,c) % negative profit so optimum is the minimum - to use with fmincon
f = -1*((r-c).*(r+r.^2-2*r.^3) + (1/3).*(((1/6)*(3*r.^4 - 4*r.^3 + 1) - c).*((r-1).^2.*(2*r+1)) + ((1/4)*(r.^2 - 1).^2 - c).*(r.^3 - 3*r+2))) ;
end
