% Solutions for Michael B Nattinger's HW1
% This program generates and plots the price dynamics given the first-order 
% difference equation discussed in the problem set given some initial price 
%
% Prepared by Fu Tan; heavily modified by Michael B Nattinger 9/2020
mkdir('pings')
clear; clc; close all;
%% Question 3
counter=1;
for p0 = [100 90 110]
%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
r = 0.01; % interest rate
d = 1; % constant dividend
% p0 = 100; % initial price
dim = 99; % terminal period t = 99

%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
pvector = zeros(dim+1,1); % creating a vector of price from t=0 to t=99
tvector = linspace(0,dim,dim+0)'; % creating a vector for time from 0 to 99
pvector(1) = p0; % giving value to the first element of the price vector
                 % with the initial price

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%

for n = 2:dim+1 % starting from t = 1 to t = 100 
    pvector(n) = (1+r)*pvector(n-1)-d; % updating the price in the next period 
                                       % with the first-order difference
                                       % equation
end

switch counter %save pvectors - will use them later in p vs t plot for (2)
    case 1
        pst1 =pvector;
    case 2
        pst2 =pvector;
    case 3
        pst3 =pvector;
    otherwise
        caution('Case did not match?')
end
counter = counter+1;
%%%%%%%%
% PLOTS
%%%%%%%%
figure();
plot(tvector(:),pvector(:),'b-');
hold on
plot(tvector(:),(d/r)+ 0*pvector(:),'k:')
hold off
title('Price Dynamics');
xlabel('Time t'); ylabel('Price p_t');
legend(['p_0 = ' num2str(p0)],'steady state p^{*}','Location','Northeast');
set(gcf,'Color',[1 1 1])
ylim([50 150])
%axis([0 dim 0 150])
cd('pings')
saveas(gcf,['dynamics_' num2str(p0) '.png'])
cd('..')
end

%% Question 4
r = 0.01; % interest rate
d = 1; % constant dividend
p0 = 100; % initial price
dim = 99; % terminal period t = 99
announ = 20;
swtch = 50; % point in time at which fed changes interest rates
rs = 0.02; % new interest rate
pvector = zeros(dim+1,1); % creating a vector of price from t=0 to t=99
tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99
pvector(1) = p0; % giving value to the first element of the price vector
                 % with the initial price
pvector(swtch) = 50;
for n = 2:announ % starting from t = 1 to t = 100 
    pvector(n) = (1+r)*pvector(n-1)-d; % updating the price in the next period 
                                       % with the first-order difference
                                       % equation
end

for n = swtch+1:dim+1
    pvector(n) = (1+rs)*pvector(n-1)-d;
end
for n = swtch-1:-1:announ+1
   pvector(n)=(pvector(n+1)+d)/(1+r); %solve backwards
end
figure();
plot(tvector(:),pvector(:),'b-');
title('Price Dynamics: Fed Funds Rate Hike Announcement Exercise ');
xlabel('Time t'); ylabel('Price P_t');
legend('Price over time','Location','Northeast');
set(gcf,'Color',[1 1 1])
ylim([40 110])
cd('pings')
saveas(gcf,'dynamics_announcement.png')
cd('..')

%% Question 2 (Phase diagram and pt over time)
r_pd = .1; % Slope given parameters from 3 is too visually similar to
                % 1 for phase diagram to look good. So I am just using a
                % different value for r so it looks better on the phase
                % diagram
x = 0:20; %just a grid for the phase diagram
xp = x*(1+r_pd)-d; %pt+1 = (1+r)pt - d
xss = 10+0*x; %steady state value


figure
plot(x,x,'k-')
hold on
plot(x,xp,'b-.')
plot(x,xss,'r--')
xline(xss(1),'r--')
hold off
legend('p_{t+1} = p_{t}','p_{t+1} =(1+r) p_t - d','p^{*}','Location','SouthEast')
%set(gca,'xtick',[])
set(gca,'xticklabel',[])
%set(gca,'ytick',[])
set(gca,'yticklabel',[])
% ylim([99 101]) % the slope given the correct parameterization is visually 
%                % too close to 1 for this graph to look good - so zoom in 
%                % for the y dimension!
% xlim([0 200])
title('Phase diagram')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'phase.png')
cd('..')

figure
plot(tvector(:),pst1,'k')
hold on
plot(tvector(:),pst2,'r')
plot(tvector(:),pst3,'b')
hold off
legend('p_0 = p^{*}','p_0 < p^{*}','p_0 > p^{*}','Location','NorthWest')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
ylim([50 150])
title('p_t over time')
xlabel('time')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'pttime.png')
cd('..')