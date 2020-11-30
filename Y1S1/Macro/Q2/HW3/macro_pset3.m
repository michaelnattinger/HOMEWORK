mkdir('pings')
clear; close all; clc

recalc = 1;
if recalc
% Part 1: stationary labor distribution
Q = [0.85 0.15; 0.05 0.95];
P0 = [1 0];
tol=1e-6;
maxiter=1e6;
iter=1;
diff=999;
while ((iter<maxiter)&&(diff>tol))
    P = P0*Q; % distribution transition
    diff = norm(P-P0); % distribution change
    iter = iter+1;
    P0 = P;
end

% Part 2: numerical value function solution
pbeta = 0.95; % parameter values
pgamma = 3;
pr = 0.03;
pw = 1.1;
l = [0.7 1.1];
abar = 3; % upper bound for asset grid
agr = 5e-4;
assets = agr:agr:abar; % asset grid
na = length(assets);
V0 = ones(na,2); % value function initial guess
legal = true(na,na,2);
for nl = 1:2
   for nap = 1:na
       legal(nap,:,nl) = (pw*l(nl) + (1+pr)*assets(nap) - assets')>0;
   end
end

V=0*V0;
aprime = V;
aind = V;
diff=999;
iter=1;
while ((iter<maxiter)&&(diff>tol))
    V=0*V0; %don't think I need this one
    for nl = 1:2
    for ai = 1:na
        Val = (pw*l(nl) + (1+pr)*assets(ai) - assets(legal(ai,:,nl))' ).^(1-pgamma)./(1-pgamma) + pbeta*(V0(legal(ai,:,nl),1)*Q(nl,1) + V0(legal(ai,:,nl),2)*Q(nl,2));
        [V(ai,nl),ind] = max(Val);
        aprime(ai,nl) = assets(ind);
        aind(ai,nl) = ind;
    end
    end
    iter = iter+1;
    diff = sum(abs(V(:)-V0(:)));
    V0 = V;
end
save results
else
    load results
end

figure
plot(assets,V(:,1),'k')
hold on
plot(assets,V(:,2),'r--')
hold off
title('Value function')
xlabel('Assets')
ylabel('Value')
legend('Low labor','High labor','Location','SouthEast')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'valfunc.png')
cd('..')

figure
plot(assets,aprime(:,1),'k')
hold on
plot(assets,aprime(:,2),'r--') %plot(assets,aprime(aprime(:,1)<assets,1))
hold off
title('Asset holdings for next period')
xlabel('Assets')
ylabel('a prime')
legend('Low labor','High labor','Location','SouthEast')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'aprime.png')
cd('..')

% stationary distribution
dist0 = ones(size(V));
dist0 = dist0./(sum(dist0(:)));
dist = dist0;
diff=999;
iter=1;
while ((iter<maxiter)&&(diff>tol))
    dist = 0*dist0;
    for nl = 1:2
       for ai = 1:na
           if dist0(ai,nl)>1e-14 % any weight? Otherwise skip
           targ = aind(ai,nl);%assets==aprime(ai,nl);
           dist(targ,1) = dist(targ,1)+dist0(ai,nl)*Q(nl,1);
           dist(targ,2) = dist(targ,2)+dist0(ai,nl)*Q(nl,2);
           end
       end
    end
    iter = iter+1;
    diff = sum(abs((dist0-dist)));
    dist0 = dist;
end
amarg = sum(dist,2);
ave = sum(amarg.*assets');
figure
plot(assets,amarg,'k')
title('Marginal asset distribution')
xlabel('Assets')
ylabel('Probability mass')
legend(['Mean assets: ' num2str(ave)],'Location','NorthEast')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'marginal.png')
cd('..')