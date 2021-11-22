// code for hw4

// clear workspace
clear

// import data
infile cnty log_cnty_pop log_cnty_rtl perc_pop midwest log_dist_B south kmart walmart num_stores dist_kmart dist_walmart opt1 opt2 opt3 using "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y2\IO\PS4\XMat.out"

////////1. probit regressions of entry////////
// Walmart entry (all variables included)
probit walmart kmart cnty log_cnty_pop log_cnty_rtl num_stores midwest south perc_pop log_dist_B dist_kmart dist_walmart

// Walmart entry (specifications that best fit the data)
probit walmart kmart log_cnty_pop log_cnty_rtl num_stores south perc_pop log_dist_B dist_kmart

// KMart entry (all variables included)
probit kmart walmart cnty log_cnty_pop log_cnty_rtl num_stores midwest south perc_pop log_dist_B dist_kmart dist_walmart

// KMart entry (specifications that best fit the data)
probit kmart walmart log_cnty_pop log_cnty_rtl midwest south perc_pop log_dist_B

////////2. probit regressions with instrumenting strategy////////
// Walmart entry
ivprobit walmart log_cnty_pop log_cnty_rtl num_stores south perc_pop log_dist_B dist_kmart (kmart = cnty midwest dist_walmart opt1 opt2 opt3)

// KMart entry
ivprobit kmart log_cnty_pop log_cnty_rtl midwest south perc_pop log_dist_B (walmart = cnty num_stores dist_kmart dist_walmart opt1 opt2 opt3)

////////3. Bresnahan and Reiss analysis of industry////////
gen large_stores = kmart+walmart
gen total_stores = num_stores+large_stores

// dependent variable = number of large players
oprobit large_stores log_cnty_pop log_cnty_rtl num_stores south perc_pop log_dist_B dist_kmart

// dependent variable = total number of small and large players
oprobit total_stores log_cnty_pop log_cnty_rtl num_stores south perc_pop log_dist_B dist_kmart

////////4. two-step method from Bajari et al. (2012)////////
// Walmart first stage, Kmart second stage
// first stage: find estimates using probit regression from question 1
probit walmart kmart log_cnty_pop log_cnty_rtl num_stores south perc_pop log_dist_B dist_kmart
predict walmart_hat

// second stage: use predicted value of Walmart entry in regression of Kmart entry
probit kmart walmart_hat log_cnty_pop log_cnty_rtl midwest south perc_pop log_dist_B

// Kmart first stage, Walmart second stage
// first stage: find estimates using probit regression from question 1
probit kmart walmart log_cnty_pop log_cnty_rtl midwest south perc_pop log_dist_B
predict kmart_hat

// second stage: use predicted value of Walmart entry in regression of Kmart entry
probit walmart kmart_hat log_cnty_pop log_cnty_rtl num_stores south perc_pop log_dist_B dist_kmart
