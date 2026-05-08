clear all
do Globals
do Functions

/*
Remark: This file requires downloading Paul Goldsmith-Pinkham's bartik_weight.ado file. 
See the README for how to do that.
*/
use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing

* Generate interactions and separated aggregate shifts
forvalues t = 1995/2021 {

    * Generate time fixed effects
    if `t' > 1995 gen tfe_`t' = year == `t'

    * Generate share x tfe interactions
    qui ds s_**_1990
    foreach var of varlist `r(varlist)' {
        gen t_`t'_`var' = (year == `t') * `var'
    }

    * Generate separated shift
    qui ds fg_agg_*
    foreach var of varlist `r(varlist)' {

        if substr("`var'", 8, 3) != "LOO" {
            gen t_`t'_`var'b = `var' if year == `t'
            egen t_`t'_`var' = max(t_`t'_`var'b), by(state)
            drop t_`t'_`var'b
        }
    }

}

* Create state fixed effects
forvalues s = 2/51 {
    gen sfe_`s' = state == `s'
}

/************************************************************
Rotemberg weights for baseline regressions
************************************************************/
xtset state year
sort state year
loc horizon = 9

* Create lhs - long differences
forvalues h = 0/`horizon' {

    loc operator "Delta`h'"
    gen `operator'Z = log(F`h'.Z/L.Z)
}

gen logZ = log(Z)
loc controls ""
forvalues lag = 1/4 { // aggregate shifts lags
    gen L`lag'_fg = L`lag'.fg
    loc controls "`controls' L`lag'_fg"
}
forvalues lag = 1/3 { // depvar lags
    gen L`lag'_D_logZ = L`lag'.D.logZ
    loc controls "`controls' L`lag'_D_logZ"
}


* Save the output to a matrix then a data frame
preserve
    bartik_weight, y(Delta0Z) x(fg) z(t_**_s_**_1990) weightstub(t_**_fg_*) controls(`controls' tfe_* sfe_*)
restore

mat beta = r(beta)
mat alpha = r(alpha)
mat gamma = r(gam)
mat pi = r(pi)
mat G = r(G)

* Going to loop through each instrument
qui desc t**_s_**_1990, varlist
local varlist = r(varlist)

frame create Rotemberg 
frame Rotemberg {
    
    loc K = rowsof(beta)
    set obs `K'
    svmat beta
    svmat alpha
    svmat gamma
    svmat pi
    svmat G

    * Label migrant goup and year
    gen group = ""
    gen year = ""
    local t = 1
    foreach var in `varlist' {
        
        replace year = substr("`var'", 3, 4) if _n == `t'
        replace group = subinstr(substr("`var'", 10, .), "_1990", "", 1) if _n == `t'
        loc ++t
    }

}
frame change Rotemberg