use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS1\AB1991.dta", clear

xtabond k, lags(1) vce(robust)

xtdpdsys k, lags(1) vce(robust)

