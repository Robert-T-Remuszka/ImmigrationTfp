clear all
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing

/*****************************
    Estimate Responses
*****************************/
EstimateIRF Zg , endogenous(fg) instruments(Bartik_1990) exogenous(l.fg) depvarlags(1) absorb(year statefip) wt(emp) impulse(fg) ///
errtype(cluster year) framename(Z_Iv1990) suffix(Iv1990) samp(`samp')

EstimateIRF Wage_Foreign , endogenous(fg) instruments(Bartik_1990) exogenous(l.fg) depvarlags(1) absorb(year statefip) wt(emp) impulse(fg) ///
errtype(cluster year) framename(Wage_Foreign_Iv1990) suffix(Iv1990) samp(`samp')

EstimateIRF Wage_Domestic , endogenous(fg) instruments(Bartik_1990) exogenous(l.fg) depvarlags(1) absorb(year statefip) wt(emp) impulse(fg) ///
errtype(cluster year) framename(Wage_Domestic_Iv1990) suffix(Iv1990) samp(`samp')

EstimateIRF Lg , endogenous(fg) instruments(Bartik_1990) exogenous(l.fg) depvarlags(1) absorb(year statefip) wt(emp) impulse(fg) ///
errtype(cluster year) framename(L_Iv1990) suffix(Iv1990) samp(`samp')

/*****************************
            PLOTS
*****************************/
set graphics off

loc suffixes "Iv1990"
loc depvars "Z Wage_Foreign Wage_Domestic L"
loc ylabs "Z" "w{superscript:F}" "w{superscript:D}" "L"
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

            * Impulse response
            tw line Beta_`suffix' h, lc("0 147 245") lw(thick) || rarea Beta_upper Beta_lower h, fcolor(ebblue%30) lwidth(none) ///
            xlab(`hmin'(1)`hmax', nogrid) ylab(, nogrid) ytitle("Response of `ylab' (%)") xtitle("Horizon") legend(off) ///
            yline(0, lc(black%50) lp(solid)) name(`v'_Response_`suffix')

            * Save
            graph export "${Graphs}/`v'_Response_`suffix'.pdf", replace name(`v'_Response_`suffix')

            * First stage diagnostics
            if inlist("`suffix'", "Iv1990") {
                
                graph bar F_Iv1990 if h > -1, over(h) bar(1, color("0 147 245") fcolor("0 147 245")) ylab(, nogrid labsize(small)) ///
                legend(off) b1title("Horizon") ytitle("First Stage F Stat (`ylab')") name(`v'_Response_`suffix'_F)

                * Save
                graph export "${Graphs}/`v'_Response_`suffix'_F.pdf", replace name(`v'_Response_`suffix'_F)

            }

        }

        loc ++counter
    }

}

***** COMBINED GRAPHS
set graphics on

* Responses
graph combine Z_Response_Iv1990 L_Response_Iv1990 Wage_Foreign_Response_Iv1990 Wage_Domestic_Response_Iv1990, ///
rows(2) cols(2) name(Responses_Iv1990)

graph export "${Graphs}/Responses_Iv1990.pdf", replace name(Responses_Iv1990)

* First Stages
graph combine Z_Response_Iv1990_F L_Response_Iv1990_F Wage_Foreign_Response_Iv1990_F Wage_Domestic_Response_Iv1990_F, ///
rows(2) cols(2) name(Responses_Iv1990_F)

graph export "${Graphs}/Responses_Iv1990_F.pdf", replace name(Responses_Iv1990_F)

