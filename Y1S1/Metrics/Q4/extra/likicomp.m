clear; close all; clc

n = 1000;
pbeta = 10;
x = rand(n,1);
s = 2;
eps = normrnd(0,s^2,n,1);
yst = exp(x*pbeta + eps);
y = 0*yst;
y(yst>=100000) = y(yst>=100000)+1;
y(yst>=50000) = y(yst>=50000)+1;
y(yst>=20000) = y(yst>=20000)+1;
y(yst>=10000) = y(yst>=10000)+1;
y(yst<10000) = y(yst<10000)+1;
liki_bh = logliki_bh(pbeta,s,x,y);
liki_mn = logliki_mn(pbeta,s,x,y);

disp(['BH log likelihood: ' num2str(liki_bh)])
disp(['MN log likelihood: ' num2str(liki_mn)])

function liki =  logliki_bh(pbeta,s,x,y)
% log likelihood calculated according to formula in BH's solutions
p1 = normcdf((log(10000) - x*pbeta)/s);
p2 = normcdf((log(20000) - x*pbeta)/s) - normcdf((log(10000) - x*pbeta)/s);
p3 = normcdf((log(50000) - x*pbeta)/s) - normcdf((log(20000) - x*pbeta)/s);
p4 = normcdf((log(100000) - x*pbeta)/s) - normcdf((log(50000) - x*pbeta)/s);
p5 = 1 - normcdf((log(100000) - x*pbeta)/s);
liki = sum(log(p1).*(y==1) + log(p2).*(y==2)+ log(p3).*(y==3)+ log(p4).*(y==4)+ log(p5).*(y==5));
end
function liki = logliki_mn(pbeta,s,x,y)
% log likelihood calculated according to formula in my solutions
p1 = logncdf(10000./exp(x*pbeta),0,s);
p2 = logncdf(20000./exp(x*pbeta),0,s) - logncdf(10000./exp(x*pbeta),0,s);
p3 = logncdf(50000./exp(x*pbeta),0,s) - logncdf(20000./exp(x*pbeta),0,s);
p4 = logncdf(100000./exp(x*pbeta),0,s) - logncdf(50000./exp(x*pbeta),0,s);
p5  = 1-logncdf(100000./exp(x*pbeta),0,s);
liki = sum(log(p1.*(y==1) + p2.*(y==2) + p3.*(y==3) + p4.*(y==4) + p5.*(y==5)));
end