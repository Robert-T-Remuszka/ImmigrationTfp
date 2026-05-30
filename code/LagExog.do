clear all
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing
sort state year
loc dkraayband = 9
loc instrument Bartik_1990

/*******
Test Lag Exogeneity
*******/

* Regress the Bartik_1990 on lags of migration flows
loc models ""
loc vlabs ""

* Vars in the reduced form system (see LagSelect.do) - leave fg out here though since we are
* going to sequentially test it in this file 
loc vars "Z Wage_Domestic Wage_Foreign L CapStock"
loc vars_fd ""
foreach v in `vars' {
    
    if "`v'" != "fg" {
        gen D0`v' = ln(`v' / L.`v')
        loc vars_fd "`vars_fd' D0`v' "
    }
    else {
        loc vars_fd "`vars_fd' fg "
    }

}

* Regress the instrument on lags of endogenous migration flows
loc maxlagtest 4
loc maxlag_other 3 // this comes from LagSelect.do
foreach lag of numlist 1/`maxlagtest' {
    
    eststo m`lag': qui ivreg2 `instrument' l(1/`maxlag_other').(`vars_fd') l(1/`lag').fg i.year i.state [pw = emp] if `samp', dkraay(`dkraayband') partial(i.year i.state)
    
    loc models "`models' m`lag'"
    if `lag' > 1 loc vlabs `vlabs' L`lag'.fg  "Lag `lag' of migration flow"
    else loc vlabs `vlabs' L.fg "Lag 1 of migration flow"

}

esttab `models' using "${Tables}/Lag_exog.tex", replace booktabs varlabels(`vlabs') se label ///
stats(N r2_a, fmt(%6.0fc %9.3f %9.3f)) ///
subs("Standard errors in parentheses" "Driscoll-Kraay standard errors with bandwidth set to `dkraayband'. All regressions include state and year fixed effects." ///
"N" "Observations" "r2_a" "Adj. \$R^2$") star(* 0.1 ** 0.05 *** 0.01) keep(L.fg L2.fg L3.fg L4.fg)

