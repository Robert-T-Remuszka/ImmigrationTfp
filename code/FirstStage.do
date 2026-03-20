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
loc maxlag     = 4

la var fg "\$\Delta^0L^F_t"

/*******
First Stage
*******/

* Build variable labels: fg lags first, then Bartik regressors
loc vlabs ""
forvalues lag = 1/`maxlag' {
    if `lag' == 1 loc vlabs `vlabs' L.fg     "Lag 1 of \$\Delta^0 \ln L^F_t$"
    else          loc vlabs `vlabs' L`lag'.fg "Lag `lag' of \$\Delta^0 \ln L^F_t$"
}
loc vlabs `vlabs' Bartik_1990   "Bartik IV"
loc vlabs `vlabs' L.Bartik_1990 "Lag 1 of Bartik IV"

loc models ""

* Columns 1-4: one lag of Bartik, varying lags of fg
forvalues lag = 1/`maxlag' {

    eststo m`lag': qui ivreg2 fg Bartik_1990 l(1/`lag').fg i.year i.state ///
        [pw = emp] if `samp', dkraay(`dkraayband') partial(i.year i.state)

    loc models "`models' m`lag'"

}

* Columns 5-8: two lags of Bartik, varying lags of fg
forvalues lag = 1/`maxlag' {

    eststo m`=`lag'+`maxlag'': qui ivreg2 fg Bartik_1990 L.Bartik_1990 l(1/`lag').fg i.year i.state ///
        [pw = emp] if `samp', dkraay(`dkraayband') partial(i.year i.state)

    loc models "`models' m`=`lag'+`maxlag''"

}

esttab `models' using "${Tables}/FirstStage.tex", replace booktabs ///
    keep(L.fg L2.fg L3.fg L4.fg Bartik_1990 L.Bartik_1990) ///
    order(L.fg L2.fg L3.fg L4.fg Bartik_1990 L.Bartik_1990) ///
    varlabels(`vlabs') se label nonumbers ///
    mgroups("Without Bartik Lag" "With Bartik Lag", pattern(1 0 0 0 1 0 0 0) span prefix(\multicolumn{@span}{c}{) suffix(}) erepeat(\cmidrule(lr){@span})) ///
    stats(N r2_a, fmt(%6.0fc %9.3f)) ///
    subs("Standard errors in parentheses" "Driscoll-Kraay standard errors with bandwidth set to `dkraayband'. All regressions include state and year fixed effects." ///
    "N" "Observations" "r2_a" "Adj. \$R^2$") star(* 0.1 ** 0.05 *** 0.01)

/**********************************************************************
Lag Selection for Dependent Variables
**********************************************************************/

loc models ""
loc depvars "Z Wage_Domestic Wage_Foreign"
loc varlabs "\$\Delta^0\text{ln}Z_t$ \$\Delta^0\text{ln}w^D_t$ \$\Delta^0\text{ln}w^F_t$"
loc i = 1

foreach y in `depvars' {

    gen ln`y' = log(`y')
    loc varlab: word `i' of `varlabs'
    la var ln`y' "`varlab'"
    loc ++i

    foreach lag of numlist 2/4 {

        eststo `y'_L`lag': qui ivreg2 fg Bartik_1990 L.Bartik_1990 L(1/`lag').D.ln`y' L(1/4).fg i.year i.state ///
            [pw = emp] if `samp', dkraay(`dkraayband') partial(i.year i.state)

        loc models "`models' `y'_L`lag'"

    }

}

esttab `models' using "${Tables}/DepvarLags.tex", replace booktabs label se nocons nonumbers ///
    drop(*.fg Bartik_1990 L.Bartik_1990) ///
    stats(N r2_a, fmt(%6.0fc %9.3f)) ///
    subs("Standard errors in parentheses" "Driscoll-Kraay standard errors with bandwidth set to `dkraayband'. All regressions include state and year fixed effects." ///
    "N" "Observations" "r2_a" "Adj. \$R^2$" ///
    "LD." "First Lag " "L2D." "Second Lag " "L3D." "Third Lag " "L4D." "Fourth Lag ") ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    mgroups("\$\Delta^0\ln Z_t$" "\$\Delta^0\ln w^D_t$" "\$\Delta^0\ln w^F_t$", pattern(1 0 0 1 0 0 1 0 0) span ///
    prefix(\multicolumn{@span}{c}{) suffix(}) erepeat(\cmidrule(lr){@span}))
