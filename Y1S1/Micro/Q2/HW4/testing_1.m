clear; close all;clc
% first: p to 1
p = 1:10;
p = 1-10.^(-p);
q = (1-p)./3;
v = 0:0.01:1;
a = (1-p-2*q)./p;
c = (1-p-q).*((1-p-2*q)./(p))./((1-p-q) - ((1-p-2*q)./(p)).*(p+q));
b = 0*a'*v;
for i=1:length(p)
b(i,:) = a(i)*v +c(i);  
end
figure
for i=1:length(p)
%plot(v,b(i,:)); hold on
end
hold off

p = 1-p;
q = 0.5+p;
a = (1-p-2*q)./p;
c = (1-p-q).*((1-p-2*q)./(p))./((1-p-q) - ((1-p-2*q)./(p)).*(p+q));
b = 0*a'*v;
for i=1:length(p)
b(i,:) = a(i)*v +c(i);  
end
figure
for i=1:length(p)
plot(v,b(i,:)); hold on
end
hold off