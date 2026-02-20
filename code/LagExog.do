clear all
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing
sort state year
/*******
Lag Exogeneity
*******/

* Regress the Bartik_1990 on lags of migration flows
loc models ""
loc vlabs ""
loc maxlag 5
foreach lag of numlist 1/`maxlag' {
    
    eststo m`lag': qui ivreg2 Bartik_1990 l(1/`lag').fg i.year i.state [pw = emp] if `samp', dkraay(9) partial(i.year i.state)
    
    loc models "`models' m`lag'"
    if `lag' > 1 loc vlabs `vlabs' L`lag'.fg  "Lag `lag' of migration flow"
    else loc vlabs `vlabs' L.fg "Lag 1 of migration flow"

}

esttab `models' using "${Tables}/Lag_exog.tex", replace booktabs varlabels(`vlabs') se label ///
stats(N r2_a, fmt(%6.0fc %9.3f %9.3f)) ///
subs("Standard errors in parentheses" "Driscoll-Kraay standard errors with bandwidth set to nine. All regressions include state and year fixed effects." ///
"N" "Observations" "r2_a" "Adj. \$R^2$") star(* 0.1 ** 0.05 *** 0.01)

