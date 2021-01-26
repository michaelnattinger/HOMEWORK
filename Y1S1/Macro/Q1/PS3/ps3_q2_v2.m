mkdir('pings')
clear; close all; clc

grid = -1000:0.01:1000;
prat = exp(grid); % get grid of price ratios
c1st = 10./(8 + 2*(prat).^2);
c2st = 2 - prat.*c1st;
U = 10*c1st- 4*c1st.^2 + 4 - (prat.*c1st).^2;

symb = {'-' '--' '-.' ':'};
pnt = {'*' '+' 'x' 'o'};
psel = [0.25 0.5 1 2]*2;
figure
hold on
for i=1:4
    [~,i_sel] = min(abs(psel(i) - prat));
   %calc U curve
   Ux = linspace(1.25 - ((41/4 - U(i_sel))/4)^0.5,1.25 + ((41/4 - U(i_sel))/4)^0.5,1000);
   Up = 2+((41/4 - U(i_sel)) - 4*(Ux - 1.25).^2).^0.5;
   Ul = 2-((41/4 - U(i_sel)) - 4*(Ux - 1.25).^2).^0.5;
   %calc budg cons
   c1g = 0:0.01:3;
   c2g = 2 - prat(i_sel).*c1g;
   plot(Ux,Up,['k' symb{i}])
   plot(Ux,Ul,['k' symb{i}])
   plot(c1g,c2g,['b' symb{i}])
   plot(c1st(i_sel),c2st(i_sel),['r' pnt{i}])
end
plot(c1st,c2st,'r-')
hold off
title('Indifference curves and budget constraints')
ylim([0 5])
set(gcf,'Color',[1 1 1])
ylabel('c_2')
xlabel('c_1')
cd('pings')
saveas(gcf,'indiff1.png')
cd('..')

figure
plot(c1st,c2st -2,'r-')
hold on
plot([-10 10],[0 0],'k')
plot([0 0],[-10 10],'k')
hold off
xlim([-1 2])
ylim([-2 1])
set(gcf,'Color',[1 1 1])
ylabel('x_2')
xlabel('x_1')
title('Offer Curve')
cd('pings')
saveas(gcf,'off1.png')
cd('..')
