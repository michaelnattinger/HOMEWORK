use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS1\Card1995.dta", clear
gen lwage = lwage76
gen edu = ed76
gen exp = age76 - edu - 6
gen exp2per = exp^2/100
gen south = reg76r
gen urban = smsa76r
gen public = nearc4a
gen private = nearc4b
gen pubage = nearc4a*age76
gen pubage2 = nearc4a*age76^2/100

ivregress 2sls lwage exp exp2per south black urban (edu = public private), r
ivregress gmm lwage exp exp2per south black urban (edu = public private), r

estat overid

ivregress 2sls lwage exp exp2per south black urban (edu = public private pubage pubage2), r
ivregress gmm lwage exp exp2per south black urban (edu = public private pubage pubage2), r

estat overid



