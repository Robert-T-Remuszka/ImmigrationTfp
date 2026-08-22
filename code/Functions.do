/***********************************************************************************************************************************************************
Pre-processing the state analysis file to prepare for estimation of
IRFs using LP.
************************************************************************************************************************************************************/
program PreRegProcessing

    * Local total employment, and its national total (used below for the Bartik LOO shifts)
    gen emp = Supply_Total
    egen emp_agg = total(emp), by(year)

    * Migration-flow regressor: the level change in the foreign-born stock, scaled by the
    * local labor force in the prior period, N_{l,t-1} (eq. 4.4). This matches the
    * normalization in Card (2001) and Peri (2012), and is required for (4.5)/(4.6)'s
    * share-weighted decomposition to be an exact identity (the shares below use this same
    * denominator).
    gen f = Supply_Foreign
    bys statefip (year): gen fg = (f - f[_n-1]) / emp[_n-1]

    * Calculate the shares
    qui ds Supply_*
    foreach v in `r(varlist)' {

        loc region = subinstr("`v'", "Supply_", "", 1)

        if !inlist("`region'","Total", "Foreign", "Domestic", "US") {

            gen s_`region' = `v' / emp

        }

    }

    * Fixed 1990 shares: s_{m,l,1990} = L^F_{m,l,1990}/N_{l,1990}, state l's own local 1990
    * labor force (eq. 4.5/4.6), matching Card (2001) and Peri (2012). rowtotal(`vars1990')
    * sums every *1990 column (all origin-region stocks plus US1990, the domestic count)
    * within a state, giving that state's own 1990 total labor force -- the same concept
    * `emp` measures contemporaneously, and the same denominator fg uses above.
    qui ds *1990
    loc vars1990 "`r(varlist)'"
    egen emp1990 = rowtotal(`vars1990')

    foreach v in `vars1990' {

        loc region = subinstr("`v'", "1990", "", 1)

        if "`region'" != "US" {
            gen s_`region'_1990 = `v' / emp1990
            replace s_`region'_1990 = 0 if mi(s_`region'_1990)
        }

    }

    * Create aggregate shifts
    gen  emp_agg_LOO = emp_agg - emp
    ds Supply_*
    foreach v in `r(varlist)' {

        loc region = subinstr("`v'", "Supply_", "", 1)

        if !inlist("`region'", "Total", "Foreign", "Domestic", "US") {
            
            replace `v' = 0 if mi(`v')
            egen Supply_Agg_`region' = total(`v'), by(year)
            * National growth rate for group `region', lagged denominator, matching the
            * instrument definition in eq. (4.5)-(4.6)
            bys statefip (year): gen fg_agg_`region' = (Supply_Agg_`region' - Supply_Agg_`region'[_n-1]) / Supply_Agg_`region'[_n-1]

            * LOO
            gen Supply_Agg_`region'_LOO = Supply_Agg_`region' - `v'
            bys statefip (year): gen fg_agg_LOO_`region' = (Supply_Agg_`region'_LOO - Supply_Agg_`region'_LOO[_n-1]) / Supply_Agg_`region'_LOO[_n-1]
            drop Supply_Agg_`region'_LOO

        }
    }

    * Create the Bartik instruments
    loc rowtotal_LOO ""
    loc rowtotal     ""
    ds fg_agg_*
    foreach v in `r(varlist)' {
        
        if substr("`v'", 8, 3) == "LOO" {
            loc region = subinstr("`v'", "fg_agg_LOO_", "", 1)
            gen Bartik_1990_`region'_LOO = s_`region'_1990 * `v'
            loc rowtotal_LOO "`rowtotal_LOO' Bartik_1990_`region'_LOO"
        }
        
        else {
            loc region = subinstr("`v'", "fg_agg_", "", 1)
            gen Bartik_1990_`region' = s_`region'_1990 * `v'
            loc rowtotal "`rowtotal' Bartik_1990_`region'"
        }
    
    }

    egen Bartik_1990_LOO   = rowtotal(`rowtotal_LOO'), missing
    egen Bartik_1990       = rowtotal(`rowtotal'), missing
    la var Bartik_1990 "Bartik IV, 1990 Shares"
    la var Bartik_1990_LOO "Bartik IV, 1990 Shares (LOO)"
    drop `rowtotal_LOO' `rowtotal'
    
    ren state statename
    encode statefip, gen(state)
    xtset state year

    * Save this file for use in Julia when computing Rotemberg weights
    save "${Data}/StateAnalysisRegReady.dta", replace
    
end

/***********************************************************************************************************************************************************
Estimate the reponse of y to the impulse given in the option impulse.
    - namelist should be the stub of the dependent variable e.g. Zg, Lg etc
    - endogenous variables are the variables which are to be treated as endogenous e.g. fg. You should not specify any instruments
      if you do not specity any endogenous variables.
    - instruments gives an instrument for each endogenous variable (or multiple for overidentification). The order of instruments should correspond to the
      order which you specify the endogenous variables
    - exogenous gives any exogenous controls we would like to include, excluding lags of the dependent variable
    - lagorderdepvar is used to include lags of the dependent variable
    - framename indicates the name of the frame created to stroe the results
    - colname tells the routine what to suffix the point estimate and se columns with
*************************************************************************************************************************************************************/
program EstimateIRF

    syntax namelist(max = 1) [, endogenous(varlist ts) instruments(varlist ts) exogenous(varlist ts) depvarlags(varlist ts) ///
    horizon(integer 9) absorb(varlist) wt(string) framename(string) suffix(string) samp(string) se_spec(string)] impulse(varname)

    * Create a place to store the results
    cap frame drop `framename'
    frame create `framename'
    frame `framename' {
        gen h = .
        gen Beta_`suffix' = .
        gen Se_`suffix'   = .
        if "`endogenous'" != "" gen F_`suffix' = .
    }

    * Create long differences of the dependent variable
    forvalues h = 0/`horizon' {
        bys statefip (year): gen Delta`h'`1' = log(`1'[_n + `h'] / `1'[_n - 1])
    }

    * Absorb fixed effects
    loc fes ""
    foreach v in `absorb' {
        loc fes "`fes' i.`v'"
    }

    * Loop through and estimate the LP at each horizon
    forvalues h = 0/`horizon' {

        * Need `horizon' local to call the correct dependent variables in the regression
        loc horizon = "F" + string(abs(`h'))

        * Let the user know what is going on
        di "***********************************************************************************************************"
            di "Dep Var     : `1'_`horizon'"
            di "Exog.  RHS  : `depvarlags' `exogenous' "
            di "Endog. RHS  : `endogenous'"
            di "Absorbed    : `absorb'"
        di "***********************************************************************************************************"

        * Run the regression and save results in the provided frame
        sort state year
        qui ivreg2 Delta`h'`1' (`endogenous' = `instruments') `exogenous' `depvarlags' `fes' [pw = `wt'] if `samp', `se_spec'

        * Record the results in the Estimates frame
        frame `framename' {

            insobs 1
            replace h = `h' if _n == _N
            replace Beta_`suffix' = _b[`impulse']       if _n == _N
            replace Se_`suffix'   = _se[`impulse']      if _n == _N
            if "`endogenous'" != "" replace F_`suffix'    = `e(widstat)' if _n == _N

        }

    }

    drop Delta*

end