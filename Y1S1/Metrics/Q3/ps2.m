clear; close all; clc
reload = 0;
if reload
dta = readtable('AE80.csv');
save 'AE80' 'dta'
else
load 'AE80'
end

x1 = dta.morekids;