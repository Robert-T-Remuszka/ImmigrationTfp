clear all
do Globals

loc samp STATEFIP != "11"

/*================================================================
Reads data/StateAnalysisPanel.dta (built by MakeLpivPanel.do) and runs the
empirical local-projection IV used to discipline the model's structural IRFs.
All regression logic lives in this file.

Migration shock (fg) and the instruments (BartikNewLoo_1990, BartikNewLoo_2000)
are built in MakeLpivPanel.do, along with everything else in
data/StateAnalysisPanel.dta --
fg = Inflow_Foreign / L.StockForeignTot, the foreign-born gross in-migration
flow (ACS MIGPLAC1, realized between t-1 and t) scaled by the state's own
foreign-born stock as of t-1; BartikNewLoo_<year> is the <year> share of each of 11
world regions in the state's foreign-born workforce interacted with that
region's national growth rate from arrivals from abroad only, leaving out the
state's own arrivals and stock, summed across regions. The local projection
uses BartikNewLoo_1990 with lags of fg as controls. This file only does
analysis: first-stage checks, controls, and the local-projection IV itself. No
population or employment weighting is used anywhere (see MakeLpivPanel.do's
header).

The last section repeats every table and figure with the migration flow
defined as the percent change in the foreign-born stock,
dlnLF = ln(Supply_Foreign_t) - ln(Supply_Foreign_{t-1}), so the local-projection
coefficients are elasticities with respect to the foreign-born stock.
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
reghdfe fg if `samp', absorb(Year) residuals(fg_tilde)
reghdfe BartikNewLoo_1990 if `samp', absorb(Year) residuals(bartik_tilde)

twoway (scatter fg_tilde bartik_tilde, mcolor(ebblue%40) msize(small)) ///
       (lfit fg_tilde bartik_tilde, lcolor(black) lwidth(medthick)), ///
       legend(off) ytitle("Migration flow") xtitle("Shift-share instrument, 1990 enclave") ///
       xlab(,nogrid) ylab(,nogrid) name(VisualIv_1990, replace)
graph export "${Graphs}/VisualIv_1990.pdf", name(VisualIv_1990) replace

/*****************************
    Lag-exogeneity test (Stock-Watson 2018, Condition LP-IV(iii)): the
    instrument should be unforecastable by lags of the endogenous variable.
    Table, three columns, all with year FE and DC excluded:
      (1) the instrument on six lags of the migration flow
      (2) the instrument on five lags of the migration flow
      (3) the first stage: the migration flow on the instrument and five lags
    Column (3) reports the first-stage F-statistic on the instrument.
*****************************/
eststo lagexog6:   reghdfe BartikNewLoo_1990 L(1/6).fg if `samp', absorb(Year) vce(cluster state)
eststo lagexog5:   reghdfe BartikNewLoo_1990 L(1/5).fg if `samp', absorb(Year) vce(cluster state)
eststo firststage: reghdfe fg BartikNewLoo_1990 L(1/5).fg if `samp', absorb(Year) vce(cluster state)
test BartikNewLoo_1990
estadd scalar FirstStageF = r(F)

esttab lagexog6 lagexog5 firststage using "${Tables}/Lag_exog.tex", replace booktabs se label nonum ///
    mtitles("\shortstack{1990 Enclave\\Instrument}" "\shortstack{1990 Enclave\\Instrument}" "Migration Flow") ///
    keep(BartikNewLoo_1990 L.fg L2.fg L3.fg L4.fg L5.fg L6.fg) ///
    order(BartikNewLoo_1990 L.fg L2.fg L3.fg L4.fg L5.fg L6.fg) ///
    varlabels(BartikNewLoo_1990 "1990 enclave instrument" L.fg "Lag 1 of migration flow" L2.fg "Lag 2 of migration flow" L3.fg "Lag 3 of migration flow" L4.fg "Lag 4 of migration flow" L5.fg "Lag 5 of migration flow" L6.fg "Lag 6 of migration flow") ///
    stats(N r2_a FirstStageF, fmt(%6.0fc %9.3f %9.2f) labels("Observations" "Adj. \$R^2\$" "First-Stage F")) ///
    star(* 0.1 ** 0.05 *** 0.01) nonotes ///
    addnotes("Standard errors clustered by state, in parentheses." ///
    "All regressions include year fixed effects." ///
    "\sym{*} \(p<0.1\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)")

/*****************************
    Local Projection IV: ln(y_{t+h}/y_{t-1}) = year FE + fg + lags of fg + e
    fg instrumented by BartikNewLoo_1990, controlling for the first `nlags'
    lags of fg (the lags that forecast the instrument in the test above).
    Horizons 0-10, cluster by state. Results frames are named by the
    pre-period year the instrument's shares are drawn from.
*****************************/
loc nlags 5

di "*************************** BARTIKNEWLOO_1990 FIRST STAGE, LAGS OF FG (EX-DC) ***"
reghdfe fg BartikNewLoo_1990 L(1/`nlags').fg if `samp', absorb(Year) vce(cluster state)

loc horizon 9
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

    foreach v of varlist Z L Wage_Domestic_Resid Wage_Foreign_Resid {
        
        forval h = 0/`horizon' {

            cap drop D_h_y
            gen double D_h_y = ln(F`h'.`v' / L.`v')

            qui ivreg2 D_h_y (fg = BartikNewLoo_`baseyear') L(1/`nlags').fg i.Year if `samp', cluster(state)

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
        foreach v in Z L Wage_Domestic_Resid Wage_Foreign_Resid {

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

/*================================================================
    PERCENT-CHANGE FLOW: the same first stage, lag-exogeneity table and
    local projections with the migration flow measured as the percent change
    in the foreign-born stock, dlnLF. The instrument is unchanged. Outputs
    carry the suffix PctChg.
================================================================*/
gen double dlnLF = ln(Supply_Foreign) - ln(L.Supply_Foreign)
la var dlnLF "Migration flow: ln(Supply_Foreign_t) - ln(Supply_Foreign_{t-1})"

/*****************************
    First Stage
*****************************/
di "*************************** BARTIKNEWLOO_1990 FIRST STAGE, PERCENT-CHANGE FLOW (EX-DC) ***"
reghdfe dlnLF BartikNewLoo_1990 if `samp', absorb(Year) vce(cluster state)

/*****************************
    Visual IV: binned first-stage scatter. The flow and the instrument are
    residualized on year FE and the first `nlags' lags of the flow, so the
    fitted line has the slope of the first stage in Lag_exog_PctChg.tex.
    Points are means within `nbins' equal-sized bins of the residualized
    instrument.
*****************************/
loc nbins 20

reghdfe dlnLF L(1/`nlags').dlnLF if `samp', absorb(Year) residuals(dlnLF_tilde)
reghdfe BartikNewLoo_1990 L(1/`nlags').dlnLF if `samp', absorb(Year) residuals(bartik_tilde_pc)

qui reg dlnLF_tilde bartik_tilde_pc
loc fsslope = _b[bartik_tilde_pc]
qui summ bartik_tilde_pc
loc xlo = r(min)
loc xhi = r(max)

preserve
    keep if !mi(dlnLF_tilde, bartik_tilde_pc)
    xtile bin = bartik_tilde_pc, nq(`nbins')
    collapse (mean) dlnLF_tilde bartik_tilde_pc, by(bin)

    twoway (function y = `fsslope'*x, range(`xlo' `xhi') lcolor(black) lwidth(medthick)) ///
           (scatter dlnLF_tilde bartik_tilde_pc, mcolor(ebblue) msize(medium)), ///
           legend(off) ytitle("Migration flow (percent change)") xtitle("Shift-share instrument, 1990 enclave") ///
           xlab(,nogrid) ylab(,nogrid) name(VisualIvPctChg_1990, replace)
    graph export "${Graphs}/VisualIvPctChg_1990.pdf", name(VisualIvPctChg_1990) replace
restore

/*****************************
    Lag-exogeneity test and first stage, same three columns as Lag_exog.tex
*****************************/
eststo pc_lagexog6:   reghdfe BartikNewLoo_1990 L(1/6).dlnLF if `samp', absorb(Year) vce(cluster state)
eststo pc_lagexog5:   reghdfe BartikNewLoo_1990 L(1/5).dlnLF if `samp', absorb(Year) vce(cluster state)
eststo pc_firststage: reghdfe dlnLF BartikNewLoo_1990 L(1/5).dlnLF if `samp', absorb(Year) vce(cluster state)
test BartikNewLoo_1990
estadd scalar FirstStageF = r(F)

esttab pc_lagexog6 pc_lagexog5 pc_firststage using "${Tables}/Lag_exog_PctChg.tex", replace booktabs se label nonum ///
    mtitles("\shortstack{1990 Enclave\\Instrument}" "\shortstack{1990 Enclave\\Instrument}" "Migration Flow") ///
    keep(BartikNewLoo_1990 L.dlnLF L2.dlnLF L3.dlnLF L4.dlnLF L5.dlnLF L6.dlnLF) ///
    order(BartikNewLoo_1990 L.dlnLF L2.dlnLF L3.dlnLF L4.dlnLF L5.dlnLF L6.dlnLF) ///
    varlabels(BartikNewLoo_1990 "1990 enclave instrument" L.dlnLF "Lag 1 of migration flow" L2.dlnLF "Lag 2 of migration flow" L3.dlnLF "Lag 3 of migration flow" L4.dlnLF "Lag 4 of migration flow" L5.dlnLF "Lag 5 of migration flow" L6.dlnLF "Lag 6 of migration flow") ///
    stats(N r2_a FirstStageF, fmt(%6.0fc %9.3f %9.2f) labels("Observations" "Adj. \$R^2\$" "First-Stage F")) ///
    star(* 0.1 ** 0.05 *** 0.01) nonotes ///
    addnotes("Standard errors clustered by state, in parentheses." ///
    "All regressions include year fixed effects." ///
    "\sym{*} \(p<0.1\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)")

/*****************************
    Local Projection IV: ln(y_{t+h}/y_{t-1}) = year FE + dlnLF + lags of dlnLF + e
    dlnLF instrumented by BartikNewLoo_1990. The relative wage
    (domestic over foreign) and the foreign-born and domestic-born stocks
    are added as outcomes: the foreign-born response is
    one at h = 0 by construction, so every coefficient is per one percent
    increase in the foreign-born stock on impact.
*****************************/
gen double Wage_Relative = Wage_Domestic_Resid / Wage_Foreign_Resid

loc ytitles `""Labor Productivity" "Labor Aggregate" "Wage Domestic" "Wage Foreign" "Relative Wage" "Foreign-Born Labor" "Domestic-Born Labor""'
loc outcomes "Z L Wage_Domestic_Resid Wage_Foreign_Resid Wage_Relative Supply_Foreign Supply_Domestic"

frame create LpIvPctChg_Results
frame LpIvPctChg_Results {
    gen variable = ""
    gen h = .
    gen beta = .
    gen se = .
    gen first_stage_F = .
}

foreach v of local outcomes {

    forval h = 0/`horizon' {

        cap drop D_h_y
        gen double D_h_y = ln(F`h'.`v' / L.`v')

        qui ivreg2 D_h_y (dlnLF = BartikNewLoo_1990) L(1/`nlags').dlnLF i.Year if `samp', cluster(state)

        frame LpIvPctChg_Results: insobs 1
        frame LpIvPctChg_Results {
            replace variable = "`v'" if _n == _N
            replace h = `h' if _n == _N
            replace beta = _b[dlnLF] if _n == _N
            replace se   = _se[dlnLF] if _n == _N
            replace first_stage_F = e(widstat) if _n == _N
        }

        drop D_h_y
    }
}

frame LpIvPctChg_Results {
    gen double upper = beta + 1.645*se
    gen double lower = beta - 1.645*se

    loc plots ""
    loc i = 1
    foreach v of local outcomes {

        loc ytitle: word `i' of `ytitles'

        twoway (rarea lower upper h if variable == "`v'", color(ebblue%20) lwidth(none)) ///
               (line beta h if variable == "`v'", lcolor(ebblue) lwidth(thick)) ///
               (function y=0, range(0 `horizon') lcolor(black)), ///
               legend(off) xtitle("Horizon (h)") ytitle("") ///
               xlab(0(1)`horizon', nogrid labsize(small)) ylab(, nogrid) ytitle("`ytitle'") ///
               name(irfpc_`v', replace) nodraw

        loc plots "`plots' irfpc_`v'"
        loc ++i
    }

    graph combine `plots', rows(2) cols(4) xsize(12) ysize(5) name(BaselinePctChg_1990)
    graph export "${Graphs}/LpivBaselinePctChg_1990.pdf", name(BaselinePctChg_1990) replace

    list variable h beta se first_stage_F if inlist(h, 0, 1, 2, 5, 10), noobs sepby(variable)
}
