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

* Create lhs - long differences
loc horizon 9
forvalues h = 0/`horizon' {

    gen Delta`h'Z = log(F`h'.Z/L.Z)
    gen Delta`h'Wage_Domestic = log(F`h'.Wage_Domestic / L.Wage_Domestic)
    gen Delta`h'Wage_Foreign = log(F`h'.Wage_Foreign / L.Wage_Foreign)

}

* Generate lags
loc migflowlags 4
loc BartikLags  4

loc controls ""
forv lag = 1/`migflowlags' { // lags of migration rates
    gen L`lag'_fg = L`lag'.fg
    loc controls "`controls' L`lag'_fg"
}
forv lag = 1/`BartikLags' { // lags of the instrument

    gen L`lag'_Bartik1990 = L`lag'.Bartik_1990
    loc controls "`controls' L`lag'_Bartik1990"

}


* Create state fixed effects
forvalues s = 2/51 {
    gen sfe_`s' = state == `s'
}


/************************************************************
                COMPUTE ROTEMBERG WEIGHTS
Loop through each horizon and each variable. Store the results
in dataframe Rotemberg
************************************************************/
loc depvarlags  4
xtset state year
sort state year
frame create Rotemberg_appended
frame Rotemberg_appended {
    gen Variable = ""
    gen Horizon  = .
}
loc vars "Z Wage_Domestic Wage_Foreign"
loc sampstart = max(`depvarlags', `migflowlags') + 1995

foreach v in `vars' {

    gen log`v' = log(`v')
    * Generate lags of the dependent var
    loc laggeddep ""
    forv lag = 1/`depvarlags' { // depvar lags
        gen L`lag'_D_`v' = L`lag'.D.log`v'
        loc laggeddep "`laggeddep' L`lag'_D_`v'"
    }
    
    forvalues h = 0 / `horizon' {

        di "********************************************" _n "CALCULATING VARIABLE `v', HORIZON `h' Weights" _n "********************************************"

        * Generate interactions and separated aggregate shifts
        loc sampend   = 2021 - `h'
        forvalues t = `sampstart'/`sampend' {

            * Generate time fixed effects
            if `t' > `sampstart' gen tfe_`t' = year == `t'

            * Generate share x tfe interactions
            qui ds s_**_1990
            foreach var of varlist `r(varlist)' {
                gen t_`t'_`var' = (year == `t') * `var'
            }

            * Generate separated shift
            qui ds fg_agg_*
            foreach var of varlist `r(varlist)' {

                if substr("`var'", 8, 3) != "LOO" {
                    qui gen t_`t'_`var'b = `var' if year == `t'
                    egen t_`t'_`var' = max(t_`t'_`var'b), by(state)
                    drop t_`t'_`var'b
                }
            }

        }

        * Capture instrument names before bartik_weight overwrites r(varlist)
        qui ds t_*_s_*_1990
        local varlist `r(varlist)'

        * Save the output to a matrix then a data frame
        preserve
            keep if inrange(year,`sampstart', `sampend')
            qui bartik_weight, y(Delta`h'`v') x(fg) z(t_**_s_**_1990) weightstub(t_**_fg_*) controls(`controls' `laggeddep' tfe_* sfe_*) weight_var(emp)
        restore

        mat beta = r(beta)
        mat alpha = r(alpha)
        mat gamma = r(gam)
        mat pi = r(pi)
        mat G = r(G)

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

        frame Rotemberg_appended {
            
            xframeappend Rotemberg, drop
            replace Horizon = `h' if mi(Horizon)
            replace Variable = "`v'" if mi(Variable)
        }
        matrix drop beta alpha gamma pi G
        drop tfe_* t_*

    }
}

/**************************
        SUMMARY TABLE
***************************/

frame Rotemberg_appended {

    destring year, replace
    gen beta2 = alpha1 * beta1
    collapse (sum) alpha1, by(Variable group Horizon)
    
}

frame change Rotemberg_appended

/**************************
    ROTEMBERG WEIGHTS TABLE
***************************/

capture program drop post_alpha
program define post_alpha, eclass
    args b V
    ereturn post `b' `V', obs(1)
    ereturn local cmd "tabstat"
end

* Order groups by decreasing alpha (Wage_Domestic, h=0 as reference)
preserve
    keep if Variable == "Wage_Domestic" & Horizon == 0
    gsort -alpha1
    local groups ""
    forvalues i = 1/`=_N' {
        local groups "`groups' `=group[`i']'"
    }
restore
local ngroups : word count `groups'

eststo clear

foreach v in Wage_Domestic Wage_Foreign Z {
    foreach h in 0 9 {

        matrix b = J(1, `ngroups', .)
        matrix V = J(`ngroups', `ngroups', 0)
        matrix colnames b = `groups'
        matrix rownames V = `groups'
        matrix colnames V = `groups'

        local i = 1
        foreach g in `groups' {
            qui sum alpha1 if Variable == "`v'" & Horizon == `h' & group == "`g'"
            matrix b[1, `i'] = r(mean)
            local ++i
        }

        post_alpha b V
        eststo m_`v'_`h'
    }
}

esttab m_Wage_Domestic_0 m_Wage_Domestic_9 ///
       m_Wage_Foreign_0  m_Wage_Foreign_9  ///
       m_Z_0             m_Z_9             ///
    using "${Tables}/Rotemberg_alphas.tex", replace ///
    booktabs nostar nose not noobs nonumber ///
    mgroups("Domestic Wage" "Foreign Wage" "Productivity", pattern(1 0 1 0 1 0) ///
        prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
    mlabels("\$h=0\$" "\$h=9\$" "\$h=0\$" "\$h=9\$" "\$h=0\$" "\$h=9\$") ///
    coeflabels(LA "Latin America" WestEu "Western EU" AsiaOther "Asia Other" ///
               CaAuNz "Can., Aus. and Nzl." EastEu "Eastern EU") ///
    b(%9.3f)
