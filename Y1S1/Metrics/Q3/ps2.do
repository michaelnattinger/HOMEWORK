cd "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q3"
use "AE80.dta"
pause on
// Michael Nattinger, with help from Sarah Bass
local cntrl agem agefstm boy1st boy2nd blackm hispm othracem
local ys workedm weeksm1 hourswm
local ym workedd weeksd1 hourswd 

// reduced form
reg morekids samesex `cntrl'

//pause
foreach i in `ys'{
    reg `i' morekids `cntrl' // col 1
	matrix row=r(table)
	local beta1_`i'=row[1,1]
	local se1_`i'=row[2,1]
	ivregress 2sls `i' (morekids=samesex) `cntrl' // col 2
	matrix row=r(table)
	local beta2_`i'=row[1,1]
	local se2_`i'=row[2,1]
	ivregress 2sls `i' (morekids=samesex) `cntrl' if msample==1 // col 5
	matrix row=r(table)
	local beta3_`i'=row[1,1]
	local se3_`i'=row[2,1]
	local beta1_`i' = round(`beta1_`i'', .001)
	local beta2_`i' = round(`beta2_`i'', .001)
	local beta3_`i' = round(`beta3_`i'', .001)
	local se1_`i' = round(`se1_`i'', .001)
	local se2_`i' = round(`se2_`i'', .001)
	local se3_`i' = round(`se3_`i'', .001)
}

foreach i in `ym'{
    reg `i' morekids `cntrl' if msample==1 // col 7 
	matrix row=r(table)
	local beta4_`i'=row[1,1]
	local se4_`i'=row[2,1]
	ivregress 2sls `i' (morekids=samesex) `cntrl' if msample ==1 // col 8
	matrix row=r(table)
	local beta5_`i'=row[1,1]
	local se5_`i'=row[2,1]
	local beta4_`i' = round(`beta4_`i'', .001)
	local beta5_`i' = round(`beta5_`i'', .001)
	local se4_`i' = round(`se4_`i'', .001)
	local se5_`i' = round(`se5_`i'', .001)
}


file open resultsfile using "ps2_results.tex", write replace
    file write resultsfile                                                          ///
        "\begin{tabular}{c | c c c c c}"                                               _newline    ///
		            "\hline"             _newline    ///
            "& All W., OLS & All W., 2SLS & M. W., 2SLS & H. of M.W., OLS & H. of M.W., 2SLS \\"             _newline    ///
		            "\hline"             _newline    ///
            "Worked for pay & `beta1_workedm' & `beta2_workedm' & `beta3_workedm' & `beta4_workedd' & `beta5_workedd' \\"             _newline    ///
						  " & (`se1_workedm') & (`se2_workedm') & (`se3_workedm') & (`se4_workedd') & (`se5_workedd')\\"       _newline    ///
            "Weeks worked & `beta1_weeksm1' & `beta2_weeksm1' & `beta3_weeksm1' & `beta4_weeksd1' & `beta5_weeksd1' \\"             _newline    ///
						  " & (`se1_weeksm1') & (`se2_weeksm1') & (`se3_weeksm1') & (`se4_weeksd1') & (`se5_weeksd1') \\"       _newline    ///
            "Hours per week & `beta1_hourswm' & `beta2_hourswm' & `beta3_hourswm' & `beta4_hourswd' & `beta5_hourswd' \\"             _newline    ///
						  " & (`se1_hourswm') & (`se2_hourswm') & (`se3_hourswm') & (`se4_hourswd') & (`se5_hourswd')"       _newline    ///
        "\end{tabular}"
    file close resultsfile