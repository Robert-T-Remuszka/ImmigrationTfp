clear all
do Globals

/*================================================================
ACS-based build of eq. (27):
β=0.96 calibrated, single endogenous regressor = wage ratio at t+1,
instrumented with the t-dated wage ratio and the t-1-dated flow-ratio
term, route-specic FE absorbed. Input panel is ACS's AcsPiPanel.dta (built by
MakeAcsPi.do)

DECISION (2026-08-21): The
endog() test below does NOT reject OLS==IV on either nativity. We're 
using the OLS estimate as the baseline value going forward (ν^D ≈ 2.6, ν^F ≈ 2.8), 
reporting the DWH test alongside it.
================================================================*/

loc bta = 0.96   // calibrated, matches Parameters()'s β default (annual equivalent of CDP's quarterly 0.99)

use "${Data}/AcsPiPanel.dta", clear
egen pairid = group(Origin Destination)

/******************* STAY (DIAGONAL) FLOW LOOKUP, BUILT ONCE ********************/
preserve
    keep if Origin == Destination
    keep Origin t Flow_Domestic Flow_Foreign
    ren Origin State
    ren Flow_Domestic Flow_Domestic_stay
    ren Flow_Foreign  Flow_Foreign_stay
    tempfile StayFlow
    save `StayFlow'
restore

/******************* M1: ORIGIN'S OWN STAY-FLOW AT t (LHS denominator) **********/
preserve
    use `StayFlow', clear
    ren State Origin
    tempfile M1
    save `M1'
restore
merge m:1 Origin t using `M1', keep(1 3) nogen
ren Flow_Domestic_stay Flow_Domestic_ll_t
ren Flow_Foreign_stay  Flow_Foreign_ll_t

/******************* M3: DESTINATION'S OWN STAY-FLOW AT t+1 (flow term's own denom) */
preserve
    use `StayFlow', clear
    ren State Destination
    replace t = t - 1
    ren Flow_Domestic_stay Flow_Domestic_l2l2_lead
    ren Flow_Foreign_stay  Flow_Foreign_l2l2_lead
    tempfile M3
    save `M3'
restore
merge m:1 Destination t using `M3', keep(1 3) nogen

/******************* M4: SAME (Origin,Destination) PAIR'S FLOW AT t+1 (flow term's own num) */
preserve
    keep Origin Destination t Flow_Domestic Flow_Foreign
    ren Flow_Domestic Flow_Domestic_lead
    ren Flow_Foreign  Flow_Foreign_lead
    replace t = t - 1
    tempfile M4
    save `M4'
restore
merge m:1 Origin Destination t using `M4', keep(1 3) nogen

/******************* M5: SAME PAIR'S WAGES AT t (lag instrument for the wage ratio) */
preserve
    keep Origin Destination t Wage_Domestic_Origin Wage_Foreign_Origin ///
         Wage_Domestic_Dest Wage_Foreign_Dest
    ren Wage_Domestic_Origin Wage_Domestic_Origin_lag
    ren Wage_Foreign_Origin  Wage_Foreign_Origin_lag
    ren Wage_Domestic_Dest   Wage_Domestic_Dest_lag
    ren Wage_Foreign_Dest    Wage_Foreign_Dest_lag
    replace t = t + 1
    tempfile M5
    save `M5'
restore
merge m:1 Origin Destination t using `M5', keep(1 3) nogen

/******************* M6: SAME PAIR'S FLOW AT t-1 (new instrument's own numerator) */
preserve
    keep Origin Destination t Flow_Domestic Flow_Foreign
    ren Flow_Domestic Flow_Domestic_lag2
    ren Flow_Foreign  Flow_Foreign_lag2
    replace t = t + 1
    tempfile M6
    save `M6'
restore
merge m:1 Origin Destination t using `M6', keep(1 3) nogen

/******************* M7: DESTINATION'S OWN STAY-FLOW AT t-1 (new instrument's own denom) */
preserve
    use `StayFlow', clear
    ren State Destination
    replace t = t + 1
    ren Flow_Domestic_stay Flow_Domestic_l2l2_lag2
    ren Flow_Foreign_stay  Flow_Foreign_l2l2_lag2
    tempfile M7
    save `M7'
restore
merge m:1 Destination t using `M7', keep(1 3) nogen

/******************* ONE CLEAN (State, Year) -> SUPPLY LOOKUP ********************/
preserve
    keep Origin t Supply_Domestic_Origin Supply_Foreign_Origin
    ren Origin State
    ren t Year_stock
    ren Supply_Domestic_Origin Supply_Domestic
    ren Supply_Foreign_Origin  Supply_Foreign
    tempfile stockA
    save `stockA'
restore
preserve
    keep Origin Year Supply_Domestic_Origin_lead Supply_Foreign_Origin_lead
    ren Origin State
    ren Year Year_stock
    ren Supply_Domestic_Origin_lead Supply_Domestic
    ren Supply_Foreign_Origin_lead  Supply_Foreign
    tempfile stockB
    save `stockB'
restore
preserve
    keep Destination Year Supply_Domestic_Dest Supply_Foreign_Dest
    ren Destination State
    ren Year Year_stock
    ren Supply_Domestic_Dest Supply_Domestic
    ren Supply_Foreign_Dest  Supply_Foreign
    tempfile stockC
    save `stockC'
restore
preserve
    use `stockA', clear
    append using `stockB'
    append using `stockC'
    collapse (mean) Supply_Domestic Supply_Foreign, by(State Year_stock)
    tempfile StockLookup
    save `StockLookup'
restore

/* Origin's stock at t-1. StockLookup is keyed on the stock's OWN calendar year,
   so the key must be shifted UP by one. */
preserve
    use `StockLookup', clear
    ren State Origin
    replace Year_stock = Year_stock + 1
    ren Year_stock t
    ren Supply_Domestic Supply_Domestic_Origin_lag2
    ren Supply_Foreign  Supply_Foreign_Origin_lag2
    tempfile OriginStockLag2
    save `OriginStockLag2'
restore
merge m:1 Origin t using `OriginStockLag2', keep(1 3) nogen

preserve
    use `StockLookup', clear
    ren State Destination
    replace Year_stock = Year_stock + 1
    ren Year_stock t
    ren Supply_Domestic Supply_Domestic_Dest_lag2
    ren Supply_Foreign  Supply_Foreign_Dest_lag2
    tempfile DestStockLag2
    save `DestStockLag2'
restore
merge m:1 Destination t using `DestStockLag2', keep(1 3) nogen

/******************* BUILD THE REGRESSION VARIABLES *************************/
gen ln_flow_ratio_D = ln(Flow_Domestic / Flow_Domestic_ll_t)
gen ln_flow_ratio_F = ln(Flow_Foreign  / Flow_Foreign_ll_t)

gen ln_pi_lead_ratio_D = ln(Flow_Domestic_lead / Flow_Domestic_l2l2_lead) ///
                       + ln(Supply_Domestic_Dest / Supply_Domestic_Origin_lead)
gen ln_pi_lead_ratio_F = ln(Flow_Foreign_lead  / Flow_Foreign_l2l2_lead) ///
                       + ln(Supply_Foreign_Dest  / Supply_Foreign_Origin_lead)

gen lhs_D = ln_flow_ratio_D - `bta' * ln_pi_lead_ratio_D
gen lhs_F = ln_flow_ratio_F - `bta' * ln_pi_lead_ratio_F

gen ln_wage_ratio_D = ln(Wage_Domestic_Dest / Wage_Domestic_Origin)
gen ln_wage_ratio_F = ln(Wage_Foreign_Dest  / Wage_Foreign_Origin)

gen ln_wage_ratio_D_lag = ln(Wage_Domestic_Dest_lag / Wage_Domestic_Origin_lag)
gen ln_wage_ratio_F_lag = ln(Wage_Foreign_Dest_lag   / Wage_Foreign_Origin_lag)

gen ln_pi_ratio_D_lag2 = ln(Flow_Domestic_lag2 / Flow_Domestic_l2l2_lag2) ///
                       + ln(Supply_Domestic_Dest_lag2 / Supply_Domestic_Origin_lag2)
gen ln_pi_ratio_F_lag2 = ln(Flow_Foreign_lag2  / Flow_Foreign_l2l2_lag2) ///
                       + ln(Supply_Foreign_Dest_lag2  / Supply_Foreign_Origin_lag2)

/******************* SAMPLE RESTRICTION: WELL-OBSERVED CORRIDORS ONLY ***********/
loc restrict_D "Flow_Domestic_lag2>0 & Flow_Domestic>0 & Flow_Domestic_lead>0"
loc restrict_F "Flow_Foreign_lag2>0  & Flow_Foreign>0  & Flow_Foreign_lead>0"
count if `restrict_D'
count if `restrict_F'

/******************* POWER DIAGNOSTIC: HOW MUCH VARIATION IS LEFT? *************/
foreach n in D F {
    qui reghdfe lhs_`n' if `restrict_`n'', absorb(pairid) resid
    predict double e_lhs_`n' if e(sample), resid
    qui reghdfe ln_wage_ratio_`n' if `restrict_`n'', absorb(pairid) resid
    predict double e_x_`n' if e(sample), resid
}
di as text "--- pair-demeaned SDs (the variation the FE regression actually uses) ---"
summ e_lhs_D e_x_D e_lhs_F e_x_F

/******************* ESTIMATE: DOMESTIC (n = D) ******************************/
/* endog() prints the Durbin-Wu-Hausman test of OLS vs. IV. Added 2026-08-21:
   The OLS/IV gap on BOTH equations turns out to be statistically indistinguishable
   from sampling noise, so this test is reported for both, not just Foreign. */
ivreghdfe lhs_D (ln_wage_ratio_D = ln_wage_ratio_D_lag ln_pi_ratio_D_lag2) ///
    if `restrict_D', absorb(pairid) cluster(pairid) first endog(ln_wage_ratio_D)

nlcom (nu_D: `bta' / _b[ln_wage_ratio_D])
loc betaOverNu_D = _b[ln_wage_ratio_D]
loc seBoN_D      = _se[ln_wage_ratio_D]
loc nu_D = `bta' / `betaOverNu_D'

reghdfe lhs_D ln_wage_ratio_D if `restrict_D', absorb(pairid) cluster(pairid)
loc betaOverNu_D_ols = _b[ln_wage_ratio_D]
loc nu_D_ols = `bta' / `betaOverNu_D_ols'
reghdfe lhs_D ln_wage_ratio_D_lag ln_pi_ratio_D_lag2 if `restrict_D', ///
    absorb(pairid) cluster(pairid)

/******************* ESTIMATE: FOREIGN (n = F) *******************************/
ivreghdfe lhs_F (ln_wage_ratio_F = ln_wage_ratio_F_lag ln_pi_ratio_F_lag2) ///
    if `restrict_F', absorb(pairid) cluster(pairid) first endog(ln_wage_ratio_F)

nlcom (nu_F: `bta' / _b[ln_wage_ratio_F])
loc betaOverNu_F = _b[ln_wage_ratio_F]
loc seBoN_F      = _se[ln_wage_ratio_F]
loc nu_F = `bta' / `betaOverNu_F'

reghdfe lhs_F ln_wage_ratio_F if `restrict_F', absorb(pairid) cluster(pairid)
loc betaOverNu_F_ols = _b[ln_wage_ratio_F]
loc nu_F_ols = `bta' / `betaOverNu_F_ols'

di as text "β (calibrated)         = " as result `bta'
di as text "(β/ν^D)-hat  IV        = " as result `betaOverNu_D' as text "  (se " as result `seBoN_D' as text ")"
di as text "ν̂^D          IV        = " as result `nu_D'
di as text "(β/ν^D)-hat  OLS       = " as result `betaOverNu_D_ols' as text "   <- USED"
di as text "ν̂^D          OLS       = " as result `nu_D_ols' as text "   <- USED"
di as text "(β/ν^F)-hat  IV        = " as result `betaOverNu_F' as text "  (se " as result `seBoN_F' as text ")"
di as text "ν̂^F          IV        = " as result `nu_F'
di as text "(β/ν^F)-hat  OLS       = " as result `betaOverNu_F_ols' as text "   <- USED"
di as text "ν̂^F          OLS       = " as result `nu_F_ols' as text "   <- USED"

/******************* SAVE A SMALL RESULTS TABLE ******************************/
clear
set obs 2
gen str1 nativity = "D" in 1
replace nativity  = "F" in 2
gen beta_calibrated  = `bta'         in 1
replace beta_calibrated = `bta'      in 2
gen betaOverNu_hat   = `betaOverNu_D' in 1
replace betaOverNu_hat = `betaOverNu_F' in 2
gen betaOverNu_se    = `seBoN_D'      in 1
replace betaOverNu_se  = `seBoN_F'    in 2
gen nu_hat           = `nu_D'        in 1
replace nu_hat        = `nu_F'        in 2
gen betaOverNu_ols   = `betaOverNu_D_ols' in 1
replace betaOverNu_ols = `betaOverNu_F_ols' in 2
gen nu_ols           = `nu_D_ols'    in 1
replace nu_ols        = `nu_F_ols'   in 2

la var betaOverNu_hat "IV estimate of β/ν^n (pair FE, 2 lagged instruments), ACS 2001-2022"
la var betaOverNu_ols "OLS estimate of β/ν^n, same sample and pair FE -- THE ONE WE USE (see header: DWH doesn't reject OLS==IV on either nativity)"

save "${Data}/NuBetaEstimatesAcs.dta", replace
