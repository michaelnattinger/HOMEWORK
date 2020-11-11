mkdir('pings')
clear; close all; clc;

g0 = 1.02;
g1 = 1.5;
eta0 = 0.1;
eta1 = 0.9;
a = 0.5;
gamm = 0.5;
delt =0.3;
bet = 0.9;
% let f = log(x)
kbar0leta = 1/(g0^(-eta0*(1-gamm))*bet^(-1) - 1 + delt);
kbar1leta = 1/(g1^(-eta0*(1-gamm))*bet^(-1) - 1 + delt);
kbar0heta = 1/(g0^(-eta1*(1-gamm))*bet^(-1) - 1 + delt);
kbar1heta = 1/(g1^(-eta1*(1-gamm))*bet^(-1) - 1 + delt);

cbar0leta = log(kbar0leta) - delt*kbar0leta;
cbar1leta = log(kbar1leta) - delt*kbar1leta;
cbar0heta = log(kbar0heta) - delt*kbar0heta;
cbar1heta = log(kbar1heta) - delt*kbar1heta;

[cshmleta,xtrajleta,ytrajleta] = calc_shm_p4(kbar0leta,cbar1leta,kbar1leta,0,0.3,100000000,delt,bet,gamm,eta0,g1);
[cshmheta,xtrajheta,ytrajheta] = calc_shm_p4(kbar0heta,cbar1heta,kbar1heta,0,0.3,100000,delt,bet,gamm,eta1,g1);

% xtrajleta = [kbar0leta kbar0leta kbar1leta];
% ytrajleta = [cbar0leta cshmleta cbar1leta];
% xtrajheta = [kbar0heta kbar0heta kbar1heta];
% ytrajheta = [cbar0heta cshmheta cbar1heta];

figure
plot([kbar0leta xtrajleta],[cbar0leta ytrajleta],'k')
hold on
plot(kbar0leta,cbar0leta,'ro')
plot(kbar0leta,cshmleta,'r+')
plot(kbar1leta,cbar1leta,'rx')
hold off
set(gcf,'Color',[1 1 1])
legend('trajectory','initial position','jump','new SS','Location','SouthEast')
xlabel('k')
ylabel('C')
title('Transition dynamics: low \eta')
cd('pings')
saveas(gcf,'loweta.png')
cd('..')

figure
plot([kbar0heta xtrajheta],[cbar0heta ytrajheta],'k')
hold on
plot(kbar0heta,cbar0heta,'ro')
plot(kbar0heta,cshmheta,'r+')
plot(kbar1heta,cbar1heta,'rx')
hold off
set(gcf,'Color',[1 1 1])
legend('trajectory','initial position','jump','new SS','Location','SouthEast')
xlabel('k')
ylabel('C')
title('Transition dynamics: high \eta')
xlim([2 4.5])
cd('pings')
saveas(gcf,'higheta.png')
cd('..')