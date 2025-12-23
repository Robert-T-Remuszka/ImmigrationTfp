/***********************************************************************************************************************************************************
Pre-processing the state analysis file to prepare for estimation of
IRFs using LP.
************************************************************************************************************************************************************/
program PreRegProcessing

    * Total employment
    gen emp = Supply_Total

    * Generate foreign-born labor share and its growth rate
    gen f = Supply_Foreign
    bys statefip (year): gen fg = (f - f[_n-1]) / emp

    * Calculate the shares
    qui ds Supply_*
    foreach v in `r(varlist)' {

        loc region = subinstr("`v'", "Supply_", "", 1)

        if !inlist("`region'","Total", "Foreign", "Domestic", "US") {

            gen s_`region' = `v' / emp
            bys statefip (year): gen s_`region'_L1 = s_`region'[_n - 1]
            bys statefip (year): gen s_`region'_L2 = s_`region'[_n - 2]

        }

    }

    * Create fixed 1990 shares
    qui ds *1990
    loc vars1990 "`r(varlist)'"
    egen emp1990 = rowtotal(`vars1990')
    foreach v in `vars1990' {
        
        loc region = subinstr("`v'", "1990", "", 1)

        if "`region'" != "US" gen s_`region'_1990 = `v' / emp1990

    }

    * Create aggregate shifts
    egen emp_agg = total(emp), by(year)
    qui ds Supply_*
    foreach v in `r(varlist)' {

        loc region = subinstr("`v'", "Supply_", "", 1)

        if !inlist("`region'", "Total", "Foreign", "Domestic", "US") {
            
            replace `v' = 0 if mi(`v')
            egen Supply_Agg_`region' = total(`v'), by(year)
            bys statefip (year): gen fg_agg_`region' = log(Supply_Agg_`region') - log(Supply_Agg_`region'[_n-1])

        }
    }

    * Create the Bartik instruments
    qui ds fg_agg_*
    foreach v in `r(varlist)' {

        loc region = subinstr("`v'", "fg_agg_", "", 1)
        gen Bartik_1990_`region' = s_`region'_1990 * `v'
        gen Bartik_L1_`region' = s_`region'_L1 * `v'
        gen Bartik_L2_`region' = s_`region'_L2 * `v'

    }

    egen Bartik_1990 = rowtotal(Bartik_1990_*), missing     // Pre-period shares
    egen Bartik_L1   = rowtotal(Bartik_L1_*),   missing     // Lagged shares
    egen Bartik_L2   = rowtotal(Bartik_L2_*),   missing
    drop Bartik_1990_* Bartik_L1_* Bartik_L2_*

    * Define left hand side variables
    forvalues h = 0/9 { // LHS variables
        loc name = "F" + string(abs(`h'))
        bys statefip (year): gen Lg_`name' = log(L[_n + `h'] / L[_n - 1])
        bys statefip (year): gen Zg_`name' = log(Z[_n + `h'] / Z[_n - 1])
        bys statefip (year): gen Wage_Foreign_`name'  = log(Wage_Foreign[_n + `h'] / Wage_Foreign[_n - 1])
        bys statefip (year): gen Wage_Domestic_`name' = log(Wage_Domestic[_n + `h'] / Wage_Domestic[_n - 1])
        bys statefip (year): gen CapStock_`name'      = log(CapStock[_n + `h'] / CapStock[_n - 1])
    }

    ren state statename
    encode statefip, gen(state)
    xtset state year

end

/***********************************************************************************************************************************************************
Estimate the reponse of y to the impulse given in the option impulse.
    - namelist should be the stub of the dependent variable e.g. Zg, Lg etc
    - endogenous variables are the variables which are to be treated as endogenous e.g. fg. You should not specify any instruments
      if you do not specity any endogenous variables.
    - instruments gives an instrument for each endogenous variable (or multiple for overidentification). The order of instruments should correspond to the
      order which you specify the endogenous variables
    - exogenous gives any exogenous controls we would like to include excluding lags of the dependent variable
    - lagorderdepvar is used to include lags of the dependent variable
    - framename indicates the name of the frame created to stroe the results
    - colname tells the routine what to suffix the point estimate and se columns with
*************************************************************************************************************************************************************/
program EstimateIRF

    syntax namelist(max = 1) [, endogenous(varlist ts) instruments(varlist ts) exogenous(varlist ts) depvarlags(integer 0) ///
    horizon(integer 9) absorb(varlist) errtype(string) wt(string) framename(string) suffix(string) samp(string)] impulse(varname)

    cap frame drop `framename'
    frame create `framename'
    frame `framename' {
        gen h = .
        gen Beta_`suffix' = .
        gen Se_`suffix'   = .
        if "`endogenous'" != "" gen F_`suffix'    = .
    }
    

    * Loop through and estimate the LP at each horizon
    forvalues h = 0/`horizon' {

        * Need `horizon' local to call the correct dependent variables in the regression
        loc horizon = "F" + string(abs(`h'))

        * Construct a local to store lags of the dependent variable
        loc yvarlags ""
        if `depvarlags' > 0 {   
            forvalues l = 1/`depvarlags' { // Following Stock Watson (2018), Piger and Stockwell (2025)
                loc yvarlags "`yvarlags' l`l'.d.`1'_`horizon'"
            }
        }

        * Let the user know what is going on
        di "***********************************************************************************************************"
            di "Dep Var     : `1'_`horizon'"
            di "Exog.  RHS  : `depvarlags' `exogenous' "
            di "Endog. RHS  : `endogenous'"
            di "Absorbed    : `absorb'"
        di "***********************************************************************************************************"
        

        * Run the regression and save results in the provided frame
        qui ivreghdfe `1'_`horizon' (`endogenous' = `instruments') `exogenous' `yvarlags' [pw = `wt'] if `samp', vce(`errtype') absorb(`absorb')

        * Record the results in the Estimates frame
        frame `framename' {

            insobs 1
            replace h = `h' if _n == _N
            replace Beta_`suffix' = _b[`impulse']       if _n == _N
            replace Se_`suffix'   = _se[`impulse']      if _n == _N
            if "`endogenous'" != "" replace F_`suffix'    = `e(widstat)' if _n == _N

        }

        
        
    }

end