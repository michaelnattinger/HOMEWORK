clear; close all; clc;

% we will define f(k)= log(k).

pbeta = 0.5;
pdelta = 0.3;
pgamma = 0.2;
peta = 0.2;
g0 = 1.4;
g1 = 1.9;

Kss0 = 1/(g0^(-peta*(1-pgamma))*(1/pbeta) - 1 + pdelta);
Css0 = log(Kss0) - pdelta*Kss0;

Kss1 = 1/(g1^(-peta*(1-pgamma))*(1/pbeta) - 1 + pdelta);
Css1 = log(Kss1) - pdelta*Kss1;

T = 50;
traj = zeros(2,T);
traj(:,1) = [Kss0; Css0];
for tt=2:50
    traj(1,tt) = (1-pdelta)*traj(1,tt-1) + log(traj(1,tt-1)) - traj(2,tt-1);
    traj(2,tt) = traj(2,tt-1)*(g1^(peta*(1-pgamma))*pbeta*(1-pdelta+(1/traj(1,tt))))^(1/pgamma);
end

figure
subplot(2,1,1)
plot(1:T,traj(1,:),'k')
legend('K')
subplot(2,1,2)
plot(1:T,traj(2,:),'k')
legend('C')
set(gcf,'Color',[1 1 1])
