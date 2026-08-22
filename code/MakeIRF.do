clear all
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing

loc depvarlags 3   // outcome's own lags, matching LagSelect.tex's BIC-preferred lag order
loc fglags     `depvarlags'   // fg's own lags, also BIC-preferred per LagSelect.tex
loc ivlags     `depvarlags'   // Bartik's own lags; significant at all three per LagExog.tex

* Standard errors clustered by state (Peri 2012's convention for this design), since
* residuals are persistent within state over time.
loc se_spec cluster(state)

* Generate macro containing first-differenced variables
loc vars "Z Wage_Domestic Wage_Foreign L"
foreach v in `vars' {
    
    gen D0`v' = ln(`v' / L.`v')

}

/*****************************
    Estimate Responses
*****************************/
* Start with the standard Bartik
EstimateIRF Z , endogenous(fg) instruments(Bartik_1990) depvarlags(L(1/`depvarlags').D0Z L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(Z_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990)

EstimateIRF Wage_Foreign , endogenous(fg) instruments(Bartik_1990) depvarlags(L(1/`depvarlags').D0Wage_Foreign L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(Wage_Foreign_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990)

EstimateIRF Wage_Domestic , endogenous(fg) instruments(Bartik_1990) depvarlags(L(1/`depvarlags').D0Wage_Domestic L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(Wage_Domestic_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990)

EstimateIRF L , endogenous(fg) instruments(Bartik_1990) depvarlags(L(1/`depvarlags').D0L L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(L_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990)

/*****************************
    SAVE IRF ESTIMATES (Iv1990 baseline) FOR JULIA'S INDIRECT INFERENCE
*****************************/
/* Stacks copies of the four Iv1990 frames (Z, Wage_Foreign, Wage_Domestic, L) into one
   long (outcome, h) panel and saves it -- operates on COPIES, not the originals, since
   the PLOTS section below still needs Beta_Iv1990/Se_Iv1990/F_Iv1990 by name to build
   the response-function and first-stage-F figures for both Iv1990 and Iv1990_LOO.
   Only the Iv1990 baseline is saved; Iv1990_LOO stays a plot-only robustness check. */
tempfile IRFStack
loc n = 0
foreach v in `vars' {
    cap frame drop IRFtmp
    frame copy `v'_Iv1990 IRFtmp
    frame IRFtmp {
        gen str14 outcome = "`v'"
        ren Beta_Iv1990 beta
        ren Se_Iv1990   se
        ren F_Iv1990    Fstat
        keep outcome h beta se Fstat
        order outcome h beta se Fstat
        loc ++n
        if `n' == 1 save `IRFStack', replace
        else {
            append using `IRFStack'
            save `IRFStack', replace
        }
    }
    frame drop IRFtmp
}

cap frame drop IRFEstimates
frame create IRFEstimates
frame IRFEstimates {
    use `IRFStack', clear
    sort outcome h
    la var outcome "Outcome variable: Z, L, Wage_Domestic or Wage_Foreign"
    la var h       "Horizon (years since the migration shock, 0-9)"
    la var beta    "LPIV coefficient, Iv1990 (Bartik) instrument"
    la var se      "State-clustered SE, Iv1990"
    la var Fstat   "First-stage Kleibergen-Paap/Wald F-stat, Iv1990"
    save "${Data}/IRFEstimates.dta", replace
}
frame drop IRFEstimates

* Look at the LOO Bartik
EstimateIRF Z , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(L(1/`depvarlags').D0Z L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(Z_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990_LOO)

EstimateIRF Wage_Foreign , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(L(1/`depvarlags').D0Wage_Foreign L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(Wage_Foreign_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990_LOO)

EstimateIRF Wage_Domestic , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(L(1/`depvarlags').D0Wage_Domestic L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(Wage_Domestic_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990_LOO)

EstimateIRF L , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(L(1/`depvarlags').D0L L(1/`fglags').fg) absorb(year) wt(emp) impulse(fg) ///
framename(L_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(`se_spec') exogenous(L(1/`ivlags').Bartik_1990_LOO)

/*****************************
            PLOTS
*****************************/
set graphics off

loc suffixes "Iv1990 Iv1990_LOO"
loc depvars "Z Wage_Foreign Wage_Domestic L"
loc ylabs "Z" "w{sup:F}" "w{sup:D}" "L"
loc counter = 1
foreach v in `depvars' {

    foreach suffix in `suffixes' {

        
        frame `v'_`suffix' {

            qui summ h
            loc hmin `r(min)'
            loc hmax `r(max)'
            gen Beta_upper = Beta_`suffix' + 1.645 * Se_`suffix'
            gen Beta_lower = Beta_`suffix' - 1.645 * Se_`suffix'
            
            loc ylab: word `counter' of "`ylabs'"
            loc yline "yline(0, lc(black%50) lp(solid))"

            * Impulse response
            tw line Beta_`suffix' h, lc("0 147 245") lw(thick) || rarea Beta_upper Beta_lower h, fcolor(ebblue%30) lwidth(none) ///
            xlab(`hmin'(1)`hmax', nogrid) ytitle("{&Delta}{sup:h}ln(`ylab')", size(large)) ylab(, nogrid) xtitle("h") legend(off) ///
            `yline' name(`v'_`suffix') 
            
            * Save
            graph export "${Graphs}/`v'_`suffix'.pdf", replace name(`v'_`suffix')

            * First stage diagnostics
            if inlist("`suffix'", "Iv1990", "Iv1990_LOO") {
                
                graph bar F_Iv1990 if h > -1, over(h) bar(1, color("0 147 245") fcolor("0 147 245")) ylab(, nogrid labsize(small)) ///
                legend(off) b1title("Horizon") ytitle("First Stage F Stat (`ylab')") name(`v'_`suffix'_F)

                * Save
                graph export "${Graphs}/`v'_`suffix'_F.pdf", replace name(`v'_`suffix'_F)

            }

        }
    
    }

    loc ++counter

}

***** COMBINED GRAPHS
set graphics on

* Responses
graph combine Z_Iv1990 L_Iv1990 Wage_Foreign_Iv1990 Wage_Domestic_Iv1990, ///
rows(2) cols(2) name(Responses_Iv1990)

graph export "${Graphs}/Responses_Iv1990.pdf", replace name(Responses_Iv1990)

graph combine Z_Iv1990_LOO L_Iv1990_LOO Wage_Foreign_Iv1990_LOO Wage_Domestic_Iv1990_LOO, ///
rows(2) cols(2) name(Responses_Iv1990_LOO)

graph export "${Graphs}/Responses_Iv1990_LOO.pdf", replace name(Responses_Iv1990_LOO)

* First Stages
graph combine Z_Iv1990_F L_Iv1990_F Wage_Foreign_Iv1990_F Wage_Domestic_Iv1990_F, ///
rows(2) cols(2) name(Responses_Iv1990_F)

graph export "${Graphs}/Responses_Iv1990_F.pdf", replace name(Responses_Iv1990_F)

graph combine Z_Iv1990_LOO_F L_Iv1990_LOO_F Wage_Foreign_Iv1990_LOO_F Wage_Domestic_Iv1990_LOO_F, ///
rows(2) cols(2) name(Responses_Iv1990_LOO_F)

graph export "${Graphs}/Responses_Iv1990_LOO_F.pdf", replace name(Responses_Iv1990_LOO_F)



