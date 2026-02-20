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

********* Make Sequential F-test tables
loc numlags 4

* A local to store the lags that we are looking to test
loc vlabs   ""
loc models  ""
loc depvars "Z Wage_Domestic"
loc varlabs \$\Delta^0\ln{Z}_t$ \$\Delta^0\ln{w}^D_t$
loc i = 1

foreach y in `depvars' {
    

    gen ln`y' = log(`y')
    loc varlab: word `i' of `varlabs'
    la var ln`y' "`varlab'"
    loc ++i

    foreach lag of numlist 1/`numlags' {
    
        eststo `y'_L`lag': qui ivreg2 fg Bartik_1990 L(1/`lag').D.ln`y' i.year i.state L(1/4).fg [pw = emp] if `samp', dkraay(9) partial(i.year i.state)

        loc models "`models' `y'_L`lag'"

    }    

    
}

esttab `models' using "${Tables}/Lag_order_Z.tex", replace booktabs label se nocons drop(*.fg Bartik_1990) ///
stats(N r2_a, fmt(%6.0fc %9.3f)) ///
subs("Standard errors in parentheses" ///
"\makecell[l]{Driscroll-Kraay standard errors with bandwidth set to nine. All regressions are employment \\weighted and include state and year fixed effects.}" ///
"N" "Observations" "r2_a" "Adj. \$R^2$" ///
"LD." "First Lag " "L2D." "Second Lag " "L3D." "Third Lag " "L4D." "Fourth Lag " "L5D." "Fifth Lag ") nomtitles ///
star(* 0.10 ** 0.05 *** 0.01) mgroups("\$Z$ Regressions" "\$w^D$ Regressions", pattern(1 0 0 0 1 0 0 0) span prefix(\multicolumn{@span}{c}{) suffix(}) erepeat(\cmidrule(lr){@span}))

/*
esttab `models' using "${Tables}/Lag_order_Z.tex", replace booktabs label se nocons drop(*.fg Bartik_1990) ///
stats(N r2_a, fmt(%6.0fc %9.3f)) ///
subs("Standard errors in parentheses" ///
"\makecell[l]{Driscroll-Kraay standard errors with bandwidth set to nine. All regressions are employment \\weighted and include state and year fixed effects.}" ///
"N" "Observations" "r2_a" "Adj. \$R^2$" ///
"LD." "First Lag " "L2D." "Second Lag " "L3D." "Third Lag " "L4D." "Fourth Lag " "L5D." "Fifth Lag ") nomtitles ///
star(* 0.10 ** 0.05 *** 0.01) mgroups("\$Z$ Regressions" "\$w^D$ Regressions", pattern(1 0 0 0 1 0 0 0) span prefix(\multicolumn{@span}{c}{) suffix(}))
*/