clear all
do Globals

loc samp STATEFIP != "11"

/*================================================================
Analysis file for the LpivRefactor pipeline: reads data/StateAnalysisPanel.dta
(built by MakeLpivPanel.do) and runs the empirical local-projection IV used to
discipline the model's structural IRFs. Unlike the legacy MakeIRF.do, this
file does not call into Functions.do -- all regression logic lives here.

Migration shock (fg) and its instrument (BartikNewLoo_1990) are both built in
MakeLpivPanel.do, along with everything else in data/StateAnalysisPanel.dta --
fg = Inflow_Foreign / L.StockForeignTot, the foreign-born gross in-migration
flow (ACS MIGPLAC1, realized between t-1 and t) scaled by the state's own
foreign-born stock as of t-1; BartikNewLoo_1990 is the 1990 share of each of 11
world regions in the state's foreign-born workforce interacted with that
region's national growth rate from arrivals from abroad only, leaving out the
state's own arrivals and stock, summed across regions. This file only does
analysis: first-stage checks, controls, and the local-projection IV itself. No
population or employment weighting is used anywhere (see MakeLpivPanel.do's
header).
================================================================*/

use "${Data}/StateAnalysisPanel.dta", clear

encode STATEFIP, gen(state)
xtset state Year
gen emp = Supply_Domestic + Supply_Foreign

/*****************************
    First Stage
*****************************/
corr fg BartikNewLoo_1990 BartikNewLoo_2000
summ fg BartikNewLoo_1990 BartikNewLoo_2000

* Year FE, cluster-robust by state - within state serial correlation
* DC excluded throughout (decided earlier: it's a small-denominator outlier --
* fg sits at the ~99th percentile every single year -- dropping it raises F
* despite a smaller point estimate, since it was adding noise, not signal).
di "*************************** BARTIKNEWLOO_1990 FIRST STAGE (EX-DC) ***************"
reghdfe fg BartikNewLoo_1990 if `samp', absorb(Year) vce(cluster state)

di "*************************** BARTIKNEWLOO_2000 FIRST STAGE (EX-DC) ***************"
reghdfe fg BartikNewLoo_2000 if `samp', absorb(Year) vce(cluster state)

/*****************************
    Visual IV: first-stage scatter, within-year variation only
*****************************/
reghdfe fg, absorb(Year) residuals(fg_tilde)
reghdfe BartikNewLoo_1990, absorb(Year) residuals(bartik_tilde)

twoway (scatter fg_tilde bartik_tilde, mcolor(%40) msize(small)) ///
       (lfit fg_tilde bartik_tilde, lcolor(red) lwidth(medthick)), ///
       legend(off) ytitle("Migration rate") xtitle("Bartik new-arrivals leave-one-out 1990") ///
       xlab(,nogrid) ylab(,nogrid) note("Residualized on time fixed effects")

/*****************************
    Lag-exogeneity test (Stock-Watson 2018, Condition LP-IV(iii)): the
    instrument should be unforecastable by lags of the endogenous variable.
    Year FE only, no other controls.
*****************************/
di "*************************** LAG-EXOGENEITY TEST (BARTIK_1990) *******************"
reghdfe BartikNewLoo_1990 L(1/5).fg, absorb(Year) vce(cluster state)
test L.fg L2.fg L3.fg L4.fg L5.fg

* Table: the same test for each version of the instrument (dependent variable =
* the instrument in each column), DC excluded, with the joint F-test on the lags.
loc lagmodels ""
foreach z in Bartik_1990 BartikNew_1990 BartikNewLoo_1990 BartikNewLoo_2000 {
    eststo lag_`z': qui reghdfe `z' L(1/5).fg if `samp', absorb(Year) vce(cluster state)
    qui test L.fg L2.fg L3.fg L4.fg L5.fg
    estadd scalar Fjoint = r(F)
    estadd scalar pjoint = r(p)
    loc lagmodels "`lagmodels' lag_`z'"
}

esttab `lagmodels' using "${Tables}/Lag_exog.tex", replace booktabs se label nonum ///
    mtitles("All arrivals" "Abroad only" "\shortstack{Abroad,\\leave-one-out}" "\shortstack{Abroad, leave-one-out,\\2000 shares}") ///
    keep(L.fg L2.fg L3.fg L4.fg L5.fg) ///
    varlabels(L.fg "Lag 1 of migration flow" L2.fg "Lag 2 of migration flow" L3.fg "Lag 3 of migration flow" L4.fg "Lag 4 of migration flow" L5.fg "Lag 5 of migration flow") ///
    stats(Fjoint pjoint N r2_a, fmt(%9.2f %9.4f %6.0fc %9.3f) labels("Joint F, five lags" "Joint p-value" "Observations" "Adj. \$R^2\$")) ///
    star(* 0.1 ** 0.05 *** 0.01) nonotes ///
    addnotes("Standard errors clustered by state, in parentheses." ///
    "All regressions include year fixed effects; the sample excludes DC." ///
    "The dependent variable is the instrument named in each column." ///
    "\sym{*} \(p<0.1\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)")

/*****************************
    Local Projection IV, no lags: ln(y_{t+h}/y_{t-1}) = year FE + fg + e
    fg instrumented by BartikNewLoo_1990 or BartikNewLoo_2000. Horizons 0-10, cluster by
    state. Results frames are named by the pre-period year the instrument's
    shares are drawn from.
*****************************/
loc horizon 10
loc ytitles `""Labor Productivity" "Labor Aggregate" "Wage Domestic" "Wage Foreign""'

foreach baseyear in 1990 {

    frame create LpIv`baseyear'_Results
    frame LpIv`baseyear'_Results {
        gen variable = ""
        gen h = .
        gen beta = .
        gen se = .
        gen first_stage_F = .
    }

    foreach v of varlist Z L Wage_Domestic_Unresid Wage_Foreign_Unresid {
        
        forval h = 0/`horizon' {

            cap drop D_h_y
            gen double D_h_y = ln(F`h'.`v' / L.`v')

            qui ivreg2 D_h_y (fg = BartikNewLoo_`baseyear') i.Year if `samp', cluster(state)

            frame LpIv`baseyear'_Results: insobs 1
            frame LpIv`baseyear'_Results {
                replace variable = "`v'" if _n == _N
                replace h = `h' if _n == _N
                replace beta = _b[fg] if _n == _N
                replace se   = _se[fg] if _n == _N
                replace first_stage_F = e(widstat) if _n == _N
            }

            drop D_h_y
        }
    }

    /*****************************
        Plot: beta(h) +/- 1.645*se (90% CI), one panel per outcome variable
    *****************************/
    frame LpIv`baseyear'_Results {
        gen double upper = beta + 1.645*se
        gen double lower = beta - 1.645*se

        loc plots ""
        loc i = 1
        foreach v in Z L Wage_Domestic_Unresid Wage_Foreign_Unresid {

            loc ytitle: word `i' of `ytitles'

            twoway (rarea lower upper h if variable == "`v'", color(ebblue%20) lwidth(none)) ///
                   (line beta h if variable == "`v'", lcolor(ebblue) lwidth(thick)) ///
                   (function y=0, range(0 `horizon') lcolor(black)), ///
                   legend(off) xtitle("Horizon (h)") ytitle("") ///
                   xlab(0(1)`horizon', nogrid labsize(small)) ylab(, nogrid) ytitle("`ytitle'") ///
                   name(irf_`baseyear'_`v', replace) nodraw

            loc plots "`plots' irf_`baseyear'_`v'"
            loc ++i
        }

        graph combine `plots', rows(2) cols(2) name(Baseline_`baseyear')
        graph export "${Graphs}/LpivBaseline_`baseyear'.pdf", name(Baseline_`baseyear') replace
    }
}
