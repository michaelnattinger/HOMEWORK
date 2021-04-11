use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS4\PSS2017.dta", clear
gen lny = log(EG_total)
rename EC_c_alt x1
rename EC_d_alt x2
drop if missing(lny, x1, x2)
nl (lny = {beta} + ({nu}/{rho}))*log({alpha}*x1^{rho} + (1-{alpha})*x2^{rho}), initial(alpha 0.5 beta 1.5 nu 1 rho 0.5)