mkdir('pings')
clear; close all; clc
%% Codefile for quarter 3 macro pset 1
% Michael Nattinger, 1/26/2021
%% Define parameters and capital grid
kg = 150.5696:1:200; % capital grid for phase diagram
D0 = 0;
D1 = 1;
T = 12;
ns = 200;
psigma = 1;
palpha = 1/3;
pbeta = 0.99^(1/12);
pdelta = 0.01;
F = @(K) K.^palpha;
Fp = @(K) palpha*K.^(palpha - 1);
U = @(C) ( C.^(1-psigma)- 1)./(1-psigma);
Up = @(C) (C.^(-psigma));

%% phase diagram for question 2
delk = F(kg) - pdelta*kg - D0;
% delta c = 0 implies k is at kss
kss = (((1/pbeta) - 1 + pdelta)/palpha)^(1/(palpha - 1));
css = F(kss) - pdelta*kss - D0;
% check SS values by propogating and seeing if it is unchanged
[kcheck,ccheck] = dprop(kss,css,D0,palpha,psigma,pbeta,pdelta);
if abs(kss- kcheck)>1e-10; warning('k changed'); end
if abs(css- ccheck)>1e-10; warning('c changed'); end
% no warnings means we are ready to rumble

recalc_saddle = 0;
if recalc_saddle
    sad = calc_saddle(kg,kss,css,D0,psigma,palpha,pbeta,pdelta, ns);
    save 'saddle_data' 'sad'
else
    load 'saddle_data' 
end

figure
plot(kg,delk,'r-')
hold on
plot([kss kss],[0 10],'b-')
plot(kg,sad,'k-')
hold off
set(gcf,'Color',[1 1 1])
title('Phase diagram - saddle path')
xlabel('K')
ylabel('C')
set(gca,'xtick',[])
set(gca,'ytick',[])
legend('\Delta K = 0','\Delta C = 0','Saddle path','Location', 'SouthEast')
annotation('arrow',[0.3 0.3],[0.2 0.25])
annotation('arrow',[0.3 0.35],[0.2 0.2])
annotation('arrow',[0.3 0.3],[0.8 0.85])
annotation('arrow',[0.3 0.25],[0.8 0.8])
annotation('arrow',[0.75 0.75],[0.35 0.3])
annotation('arrow',[0.75 0.8],[0.35 0.35])
annotation('arrow',[0.7 0.7],[0.8 0.75])
annotation('arrow',[0.7 0.65],[0.8 0.8])
ylim([3.75 3.88])
cd('pings')
saveas(gcf,'phasesad.png')
cd('..')

kg = 0:0.2:1100; % capital grid for phase diagram
delk = F(kg) - pdelta*kg - D0;
ind = delk>=0;
kg = kg(ind);
delk = delk(ind);
figure
plot(kg,delk,'r-')
hold on
plot([kss kss],[0 10],'b-')
plot(kss,css,'k*')
plot(0,0,'ko')
plot(kg(end),0,'ko')
hold off
set(gcf,'Color',[1 1 1])
title('Phase diagram - equilibria')
set(gca,'xtick',[])
set(gca,'ytick',[])
ylim([0 5])
xlabel('K')
ylabel('C')
legend('\Delta K = 0','\Delta C = 0','Steady state of interest','Other steady states','Location', 'NorthEast')
cd('pings')
saveas(gcf,'phaseeq.png')
cd('..')

%% question 3 transition
ns = 600; % 50 years in the future
D = D0+zeros(1,ns); 
D(T+1) = D1; % first index is 0 so earthquake occurs at index T+1
recalc_eq = 0;
if recalc_eq
[Ktraj,Ctraj,diff] = eq_traj(kss,css,D,palpha,psigma,pbeta,pdelta);
save 'eq' 'Ktraj' 'Ctraj'
else
load 'eq'
end
% tighter grid for zoomed in phase diagram
kg = linspace(169.6,170.8,100); % capital grid for phase diagram [169.6 170.8]
delk = F(kg) - pdelta*kg - D0;
% recalc_saddle2 = 2;
% if recalc_saddle2
%     sad = calc_saddle(kg,kss,D0,psigma,palpha,pbeta,pdelta, ns);
%     save 'saddle_data2' 'sad'
% else
%     load 'saddle_data2' 
% end
% plot trajectories
t = (-99:ns) - 1;
figure
subplot(2,1,1)
plot(t,[repmat(kss,1,100) Ktraj],'b-')
hold on
plot(t,kss+ t*0,'k-')
%plot([0 0],[169.6 170.8],'m-')
plot([0 0],[169.6 170.8],'m-')
plot([12 12],[169.6 170.8],'r-')
plot(t,[repmat(kss,1,100) Ktraj],'b-')
hold off
title('Capital trajectory over time')
legend('capital','capital steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Months since news of shock (shock at 12)')
ylim([169.6 170.8])

subplot(2,1,2)
plot(t,[repmat(css,1,100) Ctraj],'b-')
hold on
plot(t,css+ t*0,'k-')
plot([0 0],[3.82 3.85],'m-')
plot([12 12],[3.82 3.85],'r-')
plot(t,[repmat(css,1,100) Ctraj],'b-')
hold off
title('Consumption trajectory over time')
legend('consumption','consumption steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Months since news of shock (shock at 12)')
ylim([3.82 3.85])
set(gcf,'Color',[1 1 1])
suptitle('Long-run return to steady state')
cd('pings')
saveas(gcf,'longrun.png')
cd('..')

figure
subplot(2,1,1)
plot(t,[repmat(kss,1,100) Ktraj],'b-')
hold on
plot(t,kss+ t*0,'k-')
plot([0 0],[169.6 170.8],'m-')
plot([12 12],[169.6 170.8],'r-')
plot(t,[repmat(kss,1,100) Ktraj],'b-')
hold off
title('Capital trajectory over time')
legend('capital','capital steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Time since news of shock (shock at T)')
set(gca,'xtick',[])
set(gca,'ytick',[])
ylim([169.6 170.8])

subplot(2,1,2)
plot(t,[repmat(css,1,100) Ctraj],'b-')
hold on
plot(t,css+ t*0,'k-')
plot([0 0],[3.82 3.85],'m-')
plot([12 12],[3.82 3.85],'r-')
plot(t,[repmat(css,1,100) Ctraj],'b-')
hold off
title('Consumption trajectory over time')
legend('consumption','consumption steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Time since news of shock (shock at T)')
ylim([3.82 3.85])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gcf,'Color',[1 1 1])
suptitle('Long-run return to steady state')
cd('pings')
saveas(gcf,'longrunnolab.png')
cd('..')

t = (-20:120);
figure
subplot(2,1,1)
plot(t,[repmat(kss,1,20) Ktraj(:,1:121)],'b-')
hold on
plot(t,kss+ t*0,'k-')
plot([0 0],[169.6 170.8],'m-')
plot([12 12],[169.6 170.8],'r-')
plot(t,[repmat(kss,1,20) Ktraj(:,1:121)],'b-')
hold off
title('Capital trajectory over time')
legend('capital','capital steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Months since news of shock (shock at 12)')
ylim([169.6 170.8])

subplot(2,1,2)
plot(t,[repmat(css,1,20) Ctraj(:,1:121)],'b-')
hold on
plot(t,css+ t*0,'k-')
plot([0 0],[3.82 3.85],'m-')
plot([12 12],[3.82 3.85],'r-')
plot(t,[repmat(css,1,20) Ctraj(:,1:121)],'b-')
hold off
title('Consumption trajectory over time')
legend('consumption','consumption steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Months since news of shock (shock at 12)')
ylim([3.82 3.85])
set(gcf,'Color',[1 1 1])
suptitle('Short-run dynamics')
cd('pings')
saveas(gcf,'shortrun.png')
cd('..')

figure
subplot(2,1,1)
plot(t,[repmat(kss,1,20) Ktraj(:,1:121)],'b-')
hold on
plot(t,kss+ t*0,'k-')
plot([0 0],[169.6 170.8],'m-')
plot([12 12],[169.6 170.8],'r-')
plot(t,[repmat(kss,1,20) Ktraj(:,1:121)],'b-')
hold off
title('Capital trajectory over time')
legend('capital','capital steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Time since news of shock (shock at T)')
ylim([169.6 170.8])
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,1,2)
plot(t,[repmat(css,1,20) Ctraj(:,1:121)],'b-')
hold on
plot(t,css+ t*0,'k-')
plot([0 0],[3.82 3.85],'m-')
plot([12 12],[3.82 3.85],'r-')
plot(t,[repmat(css,1,20) Ctraj(:,1:121)],'b-')
hold off
title('Consumption trajectory over time')
legend('consumption','consumption steady state','shock is anticipated','shock hits','Location','SouthEast')
xlabel('Time since news of shock (shock at T)')
ylim([3.82 3.85])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gcf,'Color',[1 1 1])
suptitle('Short-run dynamics')
cd('pings')
saveas(gcf,'shortrunnolab.png')
cd('..')

figure
plot(kg,delk,'r-')
hold on
plot([kss kss],[0 5],'b-')
plot(Ktraj(1),Ctraj(1),'mo')
plot(Ktraj(T),Ctraj(T),'mx')
plot(Ktraj(T+1),Ctraj(T+1),'m+')
plot(kss,css,'m*')
plot([kss Ktraj],[css Ctraj],'k-')
plot(Ktraj(1),Ctraj(1),'mo')
plot(Ktraj(T),Ctraj(T),'mx')
plot(Ktraj(T+1),Ctraj(T+1),'m+')
plot(kss,css,'m*')
hold off
set(gcf,'Color',[1 1 1])
title('Phase diagram - transition')
xlabel('K')
ylabel('C')
legend('\Delta K = 0','\Delta C = 0','initial jump','position before shock','position at shock','final position - steady state','transition path','Location', 'NorthWest')
xlim([169.6 170.8])
ylim([3.82 3.85])
cd('pings')
saveas(gcf,'phtr.png')
cd('..')

figure
plot(kg,delk,'r-')
hold on
plot([kss kss],[0 5],'b-')
plot(Ktraj(1),Ctraj(1),'mo')
plot(Ktraj(T),Ctraj(T),'mx')
plot(Ktraj(T+1),Ctraj(T+1),'m+')
plot(kss,css,'m*')
plot([kss Ktraj],[css Ctraj],'k-')
plot(Ktraj(1),Ctraj(1),'mo')
plot(Ktraj(T),Ctraj(T),'mx')
plot(Ktraj(T+1),Ctraj(T+1),'m+')
plot(kss,css,'m*')
hold off
set(gcf,'Color',[1 1 1])
title('Phase diagram - transition')
xlabel('K')
ylabel('C')
legend('\Delta K = 0','\Delta C = 0','initial jump','position before shock','position at shock','final position - steady state','transition path','Location', 'NorthWest')
xlim([169.6 170.8])
ylim([3.82 3.85])
set(gca,'xtick',[])
set(gca,'ytick',[])
cd('pings')
saveas(gcf,'phtrnolab.png')
cd('..')

close all