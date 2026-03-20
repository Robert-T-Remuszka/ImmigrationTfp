clear all
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing

loc dkraayband = 9
loc depvarlags = 3
loc migration_flow_lags "L(1/4).fg"

/*****************************
    Estimate Responses
*****************************/
* Start with the standard Bartik
EstimateIRF Z , endogenous(fg) instruments(Bartik_1990) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(Z_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990)

EstimateIRF Wage_Foreign , endogenous(fg) instruments(Bartik_1990) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(Wage_Foreign_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990)

EstimateIRF Wage_Domestic , endogenous(fg) instruments(Bartik_1990) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(Wage_Domestic_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990)

EstimateIRF L , endogenous(fg) instruments(Bartik_1990) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(L_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990)

EstimateIRF CapStock , endogenous(fg) instruments(Bartik_1990) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(CapStock_Iv1990) suffix(Iv1990) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990)

* Look at the LOO Bartik
EstimateIRF Z , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(Z_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990_LOO)

EstimateIRF Wage_Foreign , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(Wage_Foreign_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990_LOO)

EstimateIRF Wage_Domestic , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(Wage_Domestic_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990_LOO)

EstimateIRF L , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(L_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990_LOO)

EstimateIRF CapStock , endogenous(fg) instruments(Bartik_1990_LOO) depvarlags(`depvarlags') absorb(year state) wt(emp) impulse(fg) ///
framename(CapStock_Iv1990_LOO) suffix(Iv1990_LOO) samp(`samp') horizon(9) se_spec(dkraay(`dkraayband') partial(i.year i.state)) exogenous(`migration_flow_lags' L.Bartik_1990_LOO)

/*****************************
            PLOTS
*****************************/
set graphics off

loc suffixes "Iv1990 Iv1990_LOO"
loc depvars "Z Wage_Foreign Wage_Domestic L CapStock"
loc ylabs "Z" "w{sup:F}" "w{sup:D}" "L" "K"
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
graph combine Z_Iv1990 L_Iv1990 Wage_Foreign_Iv1990 Wage_Domestic_Iv1990 CapStock_Iv1990, ///
rows(3) cols(2) name(Responses_Iv1990)

graph export "${Graphs}/Responses_Iv1990.pdf", replace name(Responses_Iv1990)

graph combine Z_Iv1990_LOO L_Iv1990_LOO Wage_Foreign_Iv1990_LOO Wage_Domestic_Iv1990_LOO CapStock_Iv1990_LOO, ///
rows(3) cols(2) name(Responses_Iv1990_LOO)

graph export "${Graphs}/Responses_Iv1990_LOO.pdf", replace name(Responses_Iv1990_LOO)

* First Stages
graph combine Z_Iv1990_F L_Iv1990_F Wage_Foreign_Iv1990_F Wage_Domestic_Iv1990_F CapStock_Iv1990_F, ///
rows(3) cols(2) name(Responses_Iv1990_F)

graph export "${Graphs}/Responses_Iv1990_F.pdf", replace name(Responses_Iv1990_F)

graph combine Z_Iv1990_LOO_F L_Iv1990_LOO_F Wage_Foreign_Iv1990_LOO_F Wage_Domestic_Iv1990_LOO_F CapStock_Iv1990_LOO_F, ///
rows(3) cols(2) name(Responses_Iv1990_LOO_F)

graph export "${Graphs}/Responses_Iv1990_LOO_F.pdf", replace name(Responses_Iv1990_LOO_F)



