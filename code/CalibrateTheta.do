clear all
do Globals

/*================================================================
Backs out the capital share theta_l directly from the capital first-order
condition, using 2015 only (the same calibration year as the bilateral
migration costs).

1. delta is a single national depreciation rate: El-Shagi and Yamarik's 2015
   state-level rates collapsed to one scalar via a 2015-state-employment-
   weighted average, to keep the paper's exposition to one depreciation rate
   rather than 51 state-specific ones.

2. theta is likewise a single national capital share: (r+delta) applied to
   the ratio of NATIONAL (summed across states) 2015 capital stock to 2015
   GDP -- a ratio of national aggregates, not an average of state ratios, so
   the one resource-heavy outlier state (K/Y up to ~10) only enters via its
   own (small) share of the national totals rather than distorting a simple
   state-by-state average.

Input: data/CapByState/state_capital_yesdata21.dta (2009$ vintage -- NOT the
newer state_capital_yesdata26.dta, which is in 2017$ and not yet merged into
this project) and data/StateAnalysisPreTfp.dta (GDP, CapStock, Supply_Domestic,
Supply_Foreign).
Output: data/ThetaDelta.dta, one row per state (statefip, delta_l, theta_l) --
both delta_l and theta_l are identical across every row by construction.
================================================================*/

loc beta = 0.96                    // matches Parameters()'s default (SsSolve_Functions.jl)
loc r    = 1/`beta' - 1

/******************* NATIONAL DEPRECIATION RATE, 2015-EMPLOYMENT-WEIGHTED ****************/
frame create Dep
frame Dep {

    use "${Data}/StateAnalysisPreTfp.dta", clear
    keep if year == 2015
    keep statefip Supply_Domestic Supply_Foreign
    tempfile supply2015
    save `supply2015'

    use "${CapStock}/state_capital_yesdata21.dta", clear

    keep if year == 2015
    tostring fips, gen(statefip)
    drop fips
    replace statefip = "0" + statefip if strlen(statefip) == 1
    keep statefip deprate

    merge 1:1 statefip using `supply2015'
    assert _merge == 3
    drop _merge
    gen double emp2015 = Supply_Domestic + Supply_Foreign

    summ deprate [aw = emp2015], meanonly
    gen double delta_l = r(mean)
    la var delta_l "El-Shagi-Yamarik rate, 2015 national average weighted by 2015 state employment"

    keep statefip delta_l

    tempfile dep
    save `dep'
}

/******************* CAPITAL SHARE FROM THE CAPITAL FOC ***************************/
use "${Data}/StateAnalysisPreTfp.dta", clear
keep if year == 2015
keep statefip GDP CapStock

merge m:1 statefip using `dep'
assert _merge == 3
drop _merge

qui summ GDP, meanonly
loc GDP_agg = r(sum)
qui summ CapStock, meanonly
loc CapStock_agg = r(sum)
qui summ delta_l, meanonly
loc delta = r(mean)
loc theta = (`r' + `delta') * `CapStock_agg' / `GDP_agg'

gen double theta_l = `theta'

la var theta_l "Capital share, 2015 national aggregate (r+delta)*CapStock/GDP"
la var delta_l "El-Shagi-Yamarik rate, 2015 national average weighted by 2015 state employment"

keep statefip theta_l delta_l
summ theta_l delta_l, detail

save "${Data}/ThetaDelta.dta", replace
di as text "Wrote ${Data}/ThetaDelta.dta"
