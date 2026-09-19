clear all
do Globals

/*================================================================
Estimate bilateral migration costs (f_ll') for
the 51x51 domestic and foreign submatrices via the symmetric-cost
ratio-of-products :

    f_ll' = -(nu/2) * ln[(pi_ll'/pi_ll)*(pi_l'l/pi_l'l')]


Before doing that inversion, first need to pick WHICH YEAR (or window) of
ACS choice probabilities to treat as the steady-state calibration
snapshot. This file's first task is exactly that: check year-over-year
stability of the domestic and foreign Pi matrices, both (1) the bilateral
PATTERN (off-diagonal correlation) and (2) the aggregate migration RATE,
to find a defensible calibration period.
================================================================*/

use "${Data}/AcsPiPanel.dta", clear
gen pi_Domestic = Flow_Domestic / Supply_Domestic_Origin
gen pi_Foreign  = Flow_Foreign  / Supply_Foreign_Origin

/******************* (1) YEAR-OVER-YEAR CORRELATION, OFF-DIAGONAL PATTERN *******/
preserve
    keep Origin Destination t pi_Domestic pi_Foreign
    drop if Origin == Destination
    reshape wide pi_Domestic pi_Foreign, i(Origin Destination) j(t)
    tempfile wide
    save `wide'
restore

use `wide', clear
qui ds pi_Domestic*
loc yearvars `r(varlist)'
loc years
foreach v of loc yearvars {
    loc y = substr("`v'", strlen("pi_Domestic") + 1, .)
    loc years `years' `y'
}
loc years : list sort years
loc n : word count `years'

matrix stab = J(`n' - 1, 3, .)
forval i = 2/`n' {
    loc y0 : word `=`i' - 1' of `years'
    loc y1 : word `i' of `years'
    qui corr pi_Domestic`y0' pi_Domestic`y1'
    matrix stab[`i' - 1, 1] = `y1'
    matrix stab[`i' - 1, 2] = r(rho)
    qui corr pi_Foreign`y0' pi_Foreign`y1'
    matrix stab[`i' - 1, 3] = r(rho)
}
matrix colnames stab = year corr_domestic corr_foreign

clear
svmat stab, names(col)
la var year "Year (correlation vs. year-1)"
la var corr_domestic "corr(Pi^D_t, Pi^D_(t-1)), off-diagonal pairs"
la var corr_foreign  "corr(Pi^F_t, Pi^F_(t-1)), off-diagonal pairs"
save "${Data}/CPStabilityByYear.dta", replace

di as text "--- Year-over-year off-diagonal correlation ---"
list, sep(0)

twoway (line corr_domestic year, lcolor(navy) lwidth(medthick)) ///
       (line corr_foreign  year, lcolor(maroon) lwidth(medthick)), ///
    ytitle("corr(Pi_t, Pi_(t-1)), off-diagonal") xtitle("Year") ///
    ylabel(0(0.2)1, angle(0)) yscale(range(0 1)) ///
    legend(order(1 "Domestic" 2 "Foreign") pos(6) rows(1)) ///
    title("Stability of the bilateral choice-probability pattern") ///
    name(stability, replace)
graph export "${Graphs}/CPStabilityByYear.pdf", replace as(pdf)

/******************* (2) AGGREGATE (STOCK-WEIGHTED) MIGRATION RATE OVER TIME ****/
use "${Data}/AcsPiPanel.dta", clear
gen pi_Domestic = Flow_Domestic / Supply_Domestic_Origin
gen pi_Foreign  = Flow_Foreign  / Supply_Foreign_Origin

preserve
    keep if Origin == Destination
    gen stay_D = pi_Domestic * Supply_Domestic_Origin
    gen stay_F = pi_Foreign  * Supply_Foreign_Origin
    collapse (sum) stay_D stay_F Supply_Domestic_Origin Supply_Foreign_Origin, by(t)
    gen agg_rate_domestic = 1 - stay_D / Supply_Domestic_Origin
    gen agg_rate_foreign  = 1 - stay_F / Supply_Foreign_Origin
    keep t agg_rate_domestic agg_rate_foreign
    ren t year
    la var year "Year"
    la var agg_rate_domestic "1 - stock-weighted stay probability, domestic"
    la var agg_rate_foreign  "1 - stock-weighted stay probability, foreign"
    save "${Data}/AggMigrationRateByYear.dta", replace
restore

use "${Data}/AggMigrationRateByYear.dta", clear
di as text "--- Aggregate (stock-weighted) migration rate by year ---"
list, sep(0)

twoway (line agg_rate_domestic year, lcolor(ebblue%85) lwidth(medthick)) ///
       (line agg_rate_foreign  year, lcolor(orange%85) lwidth(medthick)) if year >= 2001, ///
    ytitle("Aggregate interstate migration rate") xtitle("") ///
    ylabel(, nogrid) xlabel(2001(4)2021, nogrid) xscale(range(2001 2021)) ///
    legend(order(1 "Domestic" 2 "Foreign") pos(6) rows(1)) ///
    name(aggrate, replace)
graph export "${Graphs}/AggMigrationRateByYear.pdf", replace as(pdf)

di as text "Wrote ${Graphs}/CPStabilityByYear.pdf and ${Graphs}/AggMigrationRateByYear.pdf"

/*================================================================
(3) INVERT F_LL' FROM THE 2015 CHOICE PROBABILITIES
================================================================*/
use "${Data}/AcsPiPanel.dta", clear
keep if t == 2015
keep Origin Destination Flow_Domestic Flow_Foreign

* Floor zero flows at half the smallest nonzero flow observed (per
* nativity). At t=2015 this affects 553/2550 (22%) domestic and
* 1597/2550 (63%) foreign off-diagonal pairs
qui summ Flow_Domestic if Flow_Domestic > 0
loc floor_D = r(min) / 2
qui summ Flow_Foreign if Flow_Foreign > 0
loc floor_F = r(min) / 2
replace Flow_Domestic = `floor_D' if Flow_Domestic == 0
replace Flow_Foreign  = `floor_F' if Flow_Foreign  == 0

* Own-stay flow (Flow_ii), one row per state
preserve
    keep if Origin == Destination
    keep Origin Flow_Domestic Flow_Foreign
    ren Flow_Domestic Flow_Domestic_ii
    ren Flow_Foreign  Flow_Foreign_ii
    tempfile stay
    save `stay'
restore
merge m:1 Origin using `stay', keep(1 3) nogen

* Destination's own-stay flow (Flow_jj)
preserve
    use `stay', clear
    ren Origin Destination
    ren Flow_Domestic_ii Flow_Domestic_jj
    ren Flow_Foreign_ii  Flow_Foreign_jj
    tempfile stay2
    save `stay2'
restore
merge m:1 Destination using `stay2', keep(1 3) nogen

* Reverse flow (Flow_ji): swap Origin/Destination labels, then merge back
* on (Origin,Destination) -- for the ORIGINAL row (i,j), this picks up
* the flow that was originally recorded at (j,i).
preserve
    keep Origin Destination Flow_Domestic Flow_Foreign
    ren Origin Dtmp
    ren Destination Origin
    ren Dtmp Destination
    ren Flow_Domestic Flow_Domestic_ji
    ren Flow_Foreign  Flow_Foreign_ji
    tempfile reverse
    save `reverse'
restore
merge 1:1 Origin Destination using `reverse', keep(1 3) nogen

* nu^D, nu^F from the already-estimated OLS values (EstimateScaleBetaAcs.do)
preserve
    use "${Data}/NuBetaEstimatesAcs.dta", clear
    qui summ nu_ols if nativity == "D"
    loc nu_D = r(mean)
    qui summ nu_ols if nativity == "F"
    loc nu_F = r(mean)
restore
di as text "nu^D = " as result `nu_D' as text "   nu^F = " as result `nu_F'

gen ratio_D = (Flow_Domestic / Flow_Domestic_ii) * (Flow_Domestic_ji / Flow_Domestic_jj)
gen ratio_F = (Flow_Foreign  / Flow_Foreign_ii)  * (Flow_Foreign_ji  / Flow_Foreign_jj)

gen f_Domestic = -(`nu_D' / 2) * ln(ratio_D)
gen f_Foreign  = -(`nu_F' / 2) * ln(ratio_F)

keep Origin Destination f_Domestic f_Foreign
la var Origin      "Origin state (l), FIPS"
la var Destination "Destination state (l'), FIPS"
la var f_Domestic  "Bilateral migration cost f_ll' (eq. 3.3), domestic-born, 2015"
la var f_Foreign   "Bilateral migration cost f_ll' (eq. 3.3), foreign-born, 2015"
sort Origin Destination

save "${Data}/BilateralCosts2015.dta", replace
export delimited "${Data}/BilateralCosts2015.csv", replace

di as text "--- f_Domestic summary ---"
summ f_Domestic, detail
di as text "--- f_Foreign summary ---"
summ f_Foreign, detail

di as text "Wrote ${Data}/BilateralCosts2015.dta and ${Data}/BilateralCosts2015.csv"
