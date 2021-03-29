mkdir('pings')
clear; close all; clc

% del = 0.01; % width of cgrid
% cg = 0:del:1-del; % cgrid
% rgrid = linspace(0,1,1e6)';
% nc = length(cg);
% r1 = 0*cg; % prealloc
% r2 = r1;
% p1 = r1;
% p2 = r1;
% for i=1:nc
%     %[p1(i),r1(i)] = min(pr1(rgrid,cg(i)));
%     [p2(i),r2(i)] = min(pr2(rgrid,cg(i)));
% end
% 
% 
% %p1 = -1*p1; % convert from negative of profit to profit
% p2 = -1*p2;
% %r1 = rgrid(r1);
% r2 = rgrid(r2);
% figure
% plot(cg,r1,'b--')
% hold on
% plot(cg,r2,'r-.')
% plot(cg,cg,'k-')
% hold off
% title('(Expected) Profit Maximizing r')
% ylabel('Optimal r')
% xlabel('c')
% legend('First-price','Second-price','c','Location','SouthEast')
% set(gcf,'Color',[1 1 1])
% cd('pings')
% saveas(gcf,'argmax.png')
% cd('..')
% figure
% plot(cg,p1,'b--')
% hold on
% plot(cg,p2,'r-.')
% hold off
% title('Expected Profit')
% legend('First-price','Second-price','Location','NorthEast')
% ylabel('Expected profit')
% xlabel('c')
% set(gcf,'Color',[1 1 1])
% cd('pings')
% saveas(gcf,'profit.png')
% cd('..')



% function f = pr1(r,c) % negative profit so optimum is the minimum - to use with fmincon
% f = -1*(( (1/3).*(1-r.^3) - c ).*(r.^2-r.^3) + ( (1/6).*(1-r).^2.*(2.*r.^3 + r.^2 + 2*r +1) - c ).*( (1/3)*(1-r).^2.*(2*r+1)) + (r.^2 -(4/3)*r.^3 +(1/3) - c ).*(r-r.^3) + ((1/12)*(1-r).^2.*(7*r.^2+2*r+3)  - c ).*(1/3).*( r.^3 - 3*r + 2 )) ;
% end

r=0.5;
grid = linspace(r,1,1000);
b10 = r+(grid-r)./2000;
b20 = b10;
tol = 1e-3;
maxiter = 1e6;
tune = 0.5;
[b1,b2] = calc_b(grid,r,b10,b20,tol,maxiter,tune);

function f = pr2(r,c) % negative profit so optimum is the minimum - to use with fmincon
f = -1*((r-c).*(r+r.^2-2*r.^3) + (1/3).*(((1/6)*(3*r.^4 - 4*r.^3 + 1) - c).*((r-1).^2.*(2*r+1)) + ((1/4)*(r.^2 - 1).^2 - c).*(r.^3 - 3*r+2))) ;
end
function [b1,b2] = calc_b(grid,r,b10,b20,tol,maxiter,tune)
b1 = 0*b10;
b2 = 0*b20;
iter = 1;
diff = 999;
ng = length(grid);
while (diff>tol)&&(iter<maxiter)
    b1(1) = r; % boundary conditions
    b2(1) = r;
    b10p = (b10(2:end) - b10(1:end-1))./(grid(2:end) - grid(1:end-1)); % derivative
    b20p = (b20(2:end) - b20(1:end-1))./(grid(2:end) - grid(1:end-1)); % derivative
    for i=2:ng
        [~,ind1] = min(b10-b10(i));
        [~,ind2] = min(b20-b20(i));
        if ind1==1; ind1=2; end
        if ind2==1; ind2=2; end
        b1(i) = grid(i) - (1/2)*b20p(ind1-1)*grid(ind1);
        b2(i) = grid(i) - b10p(ind2-1)*grid(ind2);
    end
    diff = sum(abs(b2-b20)+abs(b1-b10)); 
    iter = iter+1;
    b10 = tune*b1+(1-tune)*b10;
    b20 = tune*b2+(1-tune)*b20;
end
end