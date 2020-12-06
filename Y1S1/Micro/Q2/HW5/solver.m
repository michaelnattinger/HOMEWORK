clear; close all; clc
syms x y c
f1 = 0==(4*c*x/9) + (2*c/9) + (3*x-y)*(x+y);
f2 = 0==(4*c/9) + (2*c*y/9) -(x - 3*y)*(x+y);

soln = solve([f1 f2],[x y],'MaxDegree',3);

cgrid = exp(-10:0.1:10);
x1 = 0*cgrid; x2=x1;x3=x1;
for i=1:length(cgrid)
x1(i) = double(subs(soln.x(1),c,cgrid(i)));
x2(i) = double(subs(soln.x(2),c,cgrid(i)));
x3(i) = double(subs(soln.x(3),c,cgrid(i)));
end

figure
plot(cgrid,x1,'k')
hold on
plot(cgrid,x2,'b')
plot(cgrid,x3,'r')
hold off
ylim([0 1])
set(gcf,'Color',[1 1 1])