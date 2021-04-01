use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS3\cps09mar.dta", clear

gen lwage = log(earnings/(hours*week))
egen educmin = min(education)
egen educmax = max(education)
gen educ = (education - educmin)/(educmax - educmin)
gen educ2 = educ^2
gen educ3 = educ^3
gen educ4 = educ^4
gen educ5 = educ^5
gen educ6 = educ^6
reg lwage educ educ2-educ6, r
predict yhat, xb
predict stdy, stdp
gen y_ub = yhat + 1.96*stdy
gen y_lb = yhat - 1.96*stdy

twoway rarea y_ub y_lb educ, sort xtitle("rescale education") ytitle("log(wage)") title("6th polynomials regression") legend(off)|| line yhat educ, sort legend(off)


