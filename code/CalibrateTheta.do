clear all
do Globals

/*================================================================
Backs out the capital share theta_l directly from the capital first-order
condition.

1. delta_l is El-Shagi and Yamarik's own state-level depreciation rate, time-
   averaged over 1994-2021 (this project's analysis window). Their series is reported at 
   the state-YEAR level

2. theta_l is the state-level time average of the year-by-year ratio
   (r+delta_l)*K_lt/Y_lt, i.e. mean_t[(r+delta_l)*K_lt/Y_lt], not
   (r+delta_l)*mean_t[K_lt]/mean_t[Y_lt]. One state (a resource-heavy outlier,
   K/Y up to ~10) pulls theta_l up to ~0.69 -- left as-is for now, not
   winsorized or otherwise adjusted.

Input: data/CapByState/state_capital_yesdata21.dta (2009$ vintage -- NOT the
newer state_capital_yesdata26.dta, which is in 2017$ and not yet merged into
this project) and data/StateAnalysisPreTfp.dta (GDP, CapStock).
Output: data/ThetaDelta.dta, one row per state (statefip, delta_l, theta_l).
================================================================*/

loc beta = 0.96                    // matches Parameters()'s default (Solve_Baseline_Functions.jl)
loc r    = 1/`beta' - 1

/******************* STATE-LEVEL DEPRECIATION, TIME-AVERAGED *********************/
frame create Dep
frame Dep {

    use "${CapStock}/state_capital_yesdata21.dta", clear

    tostring fips, gen(statefip)
    drop fips
    replace statefip = "0" + statefip if strlen(statefip) == 1

    keep if inrange(year, 1994, 2021)
    keep statefip year deprate
    count if missing(deprate)

    collapse (mean) delta_l = deprate, by(statefip)
    la var delta_l "El-Shagi-Yamarik combined depreciation rate, 1994-2021 state mean"

    tempfile dep
    save `dep'
}

/******************* CAPITAL SHARE FROM THE CAPITAL FOC ***************************/
use "${Data}/StateAnalysisPreTfp.dta", clear
keep statefip year GDP CapStock

merge m:1 statefip using `dep'
assert _merge == 3
drop _merge

gen double theta_lt = (`r' + delta_l) * CapStock / GDP

collapse (mean) theta_l = theta_lt (first) delta_l = delta_l, by(statefip)

la var theta_l "Capital share, mean_t[(r+delta_l)*CapStock_lt/GDP_lt], 1994-2021"
la var delta_l "El-Shagi-Yamarik combined depreciation rate, 1994-2021 state mean"

summ theta_l delta_l, detail

save "${Data}/ThetaDelta.dta", replace
di as text "Wrote ${Data}/ThetaDelta.dta"
