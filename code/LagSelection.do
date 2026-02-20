clear all
set more off
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing
sort state year

la var Z "\$\Delta^0 Z_t$"

********* Make Sequential F-test tables
loc numlags 6

* A local to store the lags that we are looking to test
loc vlabs   ""
loc models  ""
loc testset ""
foreach lag of numlist 1/`numlags' {
    
    eststo L`lag': qui ivreg2 fg Bartik_1990 L(1/`lag').D.Z i.year i.state [pw = emp] if `samp', dkraay(9) partial(i.year i.state)

    loc testset "`testset' L`lag'.D.Z"

    loc models "`models' L`lag'"

}    


esttab `models' using "${Tables}/Lag_order_Z.tex", replace booktabs label se nocons ///
stats(N r2_a, fmt(%6.0fc %9.3f)) ///
subs("Standard errors in parentheses" ///
"\makecell[l]{Driscroll-Kraay standard errors with bandwidth set to nine. All regressions are employment weighted and include \\state and year fixed effects.}" ///
"N" "Observations" "r2_a" "Adj. \$R^2$" ///
"LD." "First Lag " "L2D." "Second Lag " "L3D." "Third Lag " "L4D." "Fourth Lag " "L5D." "Fifth Lag " "L6D." "Sixth Lag ") nomtitles ///
star(* 0.10 ** 0.05 *** 0.01)
