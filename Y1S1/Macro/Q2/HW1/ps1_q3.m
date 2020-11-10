mkdir('pings')
clear; close all; clc
param  = parameters();
[mod,ss,lin] = calcmodss(param,0,[]);
param2 = param;
param2.pR = param.pR+0.01;
tic
[mod2,ss2,lin2,tay_sim,ex_sim,sh_sim] = calcmodss(param2,1,[ss.y.k ss.y.I]');
toc

%% Plot results for pset solutions
figure
for j=1:2
subplot(2,1,j)
plot([ones(1,4)*ss.y.(mod.yn{j}) tay_sim(j,:)],'k')
hold on
plot(ones(20,1)*ss.y.(mod.yn{j}),'b:')
plot(ones(20,1)*ss2.y.(mod.yn{j}),'r:')
hold off
title(['Taylor approximation results: ' mod.yn{j}])
legend(mod.yn{j},'original ss','ss post prod. shock','Location','NorthWest')
end
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'taylor.png')
cd('..')
figure
for j=1:2
subplot(2,1,j)
plot([ones(1,4)*ss.y.(mod.yn{j}) ex_sim(j,:)],'k')
hold on
plot(ones(20,1)*ss.y.(mod.yn{j}),'b:')
plot(ones(20,1)*ss2.y.(mod.yn{j}),'r:')
hold off
title(['Exact results using taylor c_{t0}: ' mod.yn{j}])
legend(mod.yn{j},'original ss','ss post prod. shock','Location','NorthWest')
end
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'exact.png')
cd('..')
figure
for j=1:2
subplot(2,1,j)
plot([ones(1,4)*ss.y.(mod.yn{j}) sh_sim(j,:)],'k')
hold on
plot(ones(20,1)*ss.y.(mod.yn{j}),'b:')
plot(ones(20,1)*ss2.y.(mod.yn{j}),'r:')
hold off
xlim([1 20])
title(['Simulation results: ' mod.yn{j} ', shooting method'])
legend(mod.yn{j},'original ss','ss post prod. shock','Location','NorthWest')
end
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'sh.png')
cd('..')