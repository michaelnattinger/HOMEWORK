pause on
use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS3\RR2010.dta", clear
gen glag = L.gdp
gen dlag = L.debt

gen dlag60 = 0
replace dlag60 = dlag-60 if dlag>=60
gen dlag40 = 0
replace dlag40 = dlag-40 if dlag>=40
gen dlag80 = 0
replace dlag80 = dlag-80 if dlag>=80

* no spline
reg gdp glag dlag
predict yhat, xb
gen y0 = yhat
estat ic
drop yhat

reg gdp glag dlag dlag60
predict yhat, xb
predict stdy, stdp
gen y_ub = yhat + 1.96*stdy
gen y_lb = yhat - 1.96*stdy
gen y1 = yhat
estat ic
drop yhat

reg gdp glag dlag dlag40 dlag80
predict yhat, xb
estat ic
gen y2 = yhat


line y0 y1 y2 dlag, sort xtitle("lagged debt") ytitle("growth") title("three point estimates") 
pause
twoway rarea y_ub y_lb dlag, sort xtitle("lagged debt") ytitle("growth") title("one knot confidence interval") legend(off)|| line y1 dlag, sort legend(off)

