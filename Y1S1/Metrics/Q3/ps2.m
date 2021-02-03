clear; clc; close all

syms k c kp cp i ip p A palpha pdelta kbar ibar cbar
f1 = i == A*(kbar/ibar)*k - (cbar/ibar)*c;
f2 = ip == A*(kbar/ibar)*kp - (cbar/ibar)*cp;
f3 = kp == (1-pdelta)*kp + pdelta * i;
f4 = kp*(A*palpha * pdelta * (palpha - 1)/kbar - (1-pdelta)* ibar/kbar)/(A * palpha * pdelta * kbar^(palpha - 1) + (1-pdelta)* ibar/kbar) ...
    + ip*((1-pdelta)* ibar/kbar)/(A * palpha * pdelta * kbar^(palpha - 1) + (1-pdelta)* ibar/kbar) ...
    + c + (1-pdelta)*k + (pdelta - 1)*i;
sol = solve([f1 f2 f3 f4],[kp cp i ip]);
solcp = simplify(sol.cp);