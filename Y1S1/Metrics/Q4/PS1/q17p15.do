use "AB1991.dta", clear

xtabond k, lags(1) vce(robust)

xtdpdsys k, lags(1) vce(robust)

