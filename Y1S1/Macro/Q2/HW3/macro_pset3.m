mkdir('pings')
clear; close all; clc

recalc = 0;
if recalc
% Part 1: stationary labor distribution
Q = [0.85 0.15; 0.05 0.95];
P0 = [1 0];
tol=1e-4;
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
agr = 3e-4;
assets = agr:agr:abar; % asset grid
na = length(assets);
%V0 = ones(na,2); % value function initial guess
V0 = [log(assets') log(assets')];
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
tic
while ((iter<maxiter)&&(diff>tol))
    toc
    disp(['beginning iteration ' num2str(iter) ', diff = ' num2str(diff)])
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
toc
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
xlim([-0.1 3])
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
xlim([-0.1 3])
cd('pings')
saveas(gcf,'aprime.png')
cd('..')

% abar
% check for a upper bar
change = [0 0];
for nl = 1:2
    check = 1;
    for ai = 1:na
        if check>0
            if aprime(ai,nl)<=assets(ai); check = 0; change(nl) = ai; end
        end
    end
end
figure
plot(assets,aprime(:,1),'k')
hold on
plot(assets,aprime(:,2),'k')
plot(assets(change(1):end),aprime(change(1):end,1),'r-')
plot(assets(change(2):end),aprime(change(2):end,2),'r-')
plot(assets(change(1)),aprime(change(1),1),'rx')
plot(assets(change(2)),aprime(change(2),2),'rx')
plot(assets,assets,'m--')
hold off
set(gcf,'Color',[1 1 1])
legend('a prime, l_l','a prime, l_h','a prime_h < a','a prime_l < a','cutoff','cutoff','a prime = a','Location', 'SouthEast')
title('abar for high and low l')
xlabel('Assets')
ylabel('a prime')
xlim([-0.1 3])
cd('pings')
saveas(gcf,'abar.png')
cd('..')

% stationary distribution
dist0 = ones(size(V));
dist0 = dist0./(sum(dist0(:)));
%dist = dist0;
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
    diff = sum(abs((dist0(:)-dist(:))));
    dist0 = dist;
end
amarg = sum(dist,2);
ave = sum(amarg.*assets');
% sum(dist(:,1)) = 0.25 can check the following
% sum(dist(:,2)) = 0.75
% sum(amarg) = 1
figure
plot(assets,amarg,'k')
title('Marginal asset distribution')
xlabel('Assets')
ylabel('Probability mass')
legend(['Mean assets: ' num2str(ave)],'Location','NorthEast')
set(gcf,'Color',[1 1 1])
xlim([-0.1 3])
cd('pings')
saveas(gcf,'marginal.png')
cd('..')

figure
plot(assets,dist(:,1),'k')
hold on
plot(assets,dist(:,2),'r')
hold off
title('Asset distribution')
xlabel('Assets')
ylabel('Probability mass')
legend('unemployed','employed','Location','NorthEast')
set(gcf,'Color',[1 1 1])
xlim([-0.1 3])
cd('pings')
saveas(gcf,'asdist.png')
cd('..')