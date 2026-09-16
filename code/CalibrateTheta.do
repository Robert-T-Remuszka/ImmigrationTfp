clear all
do Globals

/*================================================================
Backs out the capital share theta_l directly from the capital first-order
condition.

1. delta is a single national depreciation rate: El-Shagi and Yamarik's
   state-level rates (reported at the state-year level, time-averaged over
   1994-2021) collapsed to one scalar via a state-employment-weighted
   average (1996 weights), to keep the paper's exposition to one
   depreciation rate rather than 51 state-specific ones.

2. theta is likewise a single national capital share: (r+delta) applied to
   the ratio of NATIONAL (summed across states) capital stock to GDP each
   year, i.e. mean_t[(r+delta)*(sum_l K_lt)/(sum_l Y_lt)], then time-averaged
   over 1994-2021 -- a ratio of national aggregates, not an average of state
   ratios, so the one resource-heavy outlier state (K/Y up to ~10) only
   enters via its own (small) share of the national totals rather than
   distorting a simple state-by-state average.

Input: data/CapByState/state_capital_yesdata21.dta (2009$ vintage -- NOT the
newer state_capital_yesdata26.dta, which is in 2017$ and not yet merged into
this project), data/StateAnalysisPreTfp.dta (GDP, CapStock), and data/PiMat.dta
(Domestic_1996/Foreign_1996, for the depreciation-rate employment weights).
Output: data/ThetaDelta.dta, one row per state (statefip, delta_l, theta_l) --
both delta_l and theta_l are identical across every row by construction.
================================================================*/

loc beta = 0.96                    // matches Parameters()'s default (Solve_Baseline_Functions.jl)
loc r    = 1/`beta' - 1

/******************* NATIONAL DEPRECIATION RATE, EMPLOYMENT-WEIGHTED *********************/
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

    rename statefip Origin
    merge 1:1 Origin using "${Data}/PiMat.dta", keepusing(Domestic_1996 Foreign_1996)
    // Origin==2 unmatched "from using" is expected: PiMat.dta's 52nd row is
    // Rest-of-World, which has no counterpart in the (51-state) depreciation panel.
    assert _merge != 1
    keep if _merge == 3
    drop _merge
    gen double emp1996 = Domestic_1996 + Foreign_1996

    summ delta_l [aw = emp1996], meanonly
    replace delta_l = r(mean)
    la var delta_l "El-Shagi-Yamarik rate, national average weighted by 1996 state employment"

    keep Origin delta_l
    rename Origin statefip

    tempfile dep
    save `dep'
}

/******************* CAPITAL SHARE FROM THE CAPITAL FOC ***************************/
use "${Data}/StateAnalysisPreTfp.dta", clear
keep statefip year GDP CapStock

merge m:1 statefip using `dep'
assert _merge == 3
drop _merge

preserve
    collapse (sum) GDP CapStock (mean) delta_l, by(year)
    gen double theta_t = (`r' + delta_l) * CapStock / GDP
    summ theta_t, meanonly
    loc theta = r(mean)
restore

collapse (first) delta_l, by(statefip)
gen double theta_l = `theta'

la var theta_l "Capital share, national aggregate mean_t[(r+delta)*CapStock_t/GDP_t], 1994-2021"
la var delta_l "El-Shagi-Yamarik rate, national average weighted by 1996 state employment"

summ theta_l delta_l, detail

save "${Data}/ThetaDelta.dta", replace
di as text "Wrote ${Data}/ThetaDelta.dta"
