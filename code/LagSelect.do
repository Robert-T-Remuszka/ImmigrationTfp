clear all
do Globals
do Functions


use "${Data}/StateAnalysis.dta", clear
qui PreRegProcessing
sort state year
loc highlight 1

/*********************************************************
    Estimate reduced form VAR equation by equation
**********************************************************/
loc maxlag 4
loc sampstart = 1995 + `maxlag'
loc vars "Z Wage_Domestic Wage_Foreign L CapStock fg"
loc K : word count `vars'
loc dkraayband 9

matrix BIC = J(`maxlag' + 1, `K', .)
loc rownames ""
forv p = 1/`maxlag' {
    loc rownames `"`rownames' "Lag `p'""'
}
matrix rownames BIC = `rownames' "N"
matrix colnames BIC = "\$\Delta^0\ln Z_{l,t}\$" "\$\Delta^0\ln w^D_{l,t}\$" ///
"\$\Delta^0\ln w^F_{l,t}\$" "\$\Delta^0\ln L_{l,t}\$" "\$\Delta^0\ln K_{l,t}\$" ///
"\$\Delta^0\ln L^F_{l,t}\$"

/***** 
Compute variables in the reduced form VAR representation - add them to a local
*****/
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

/*******************
Calculate BIC of a reduced form ARX model for each outcome variable.
*****************/
forv k = 1/`K' {

    
    loc v: word `k' of `vars_fd'
    
    forv p = 1/`maxlag' {

        qui ivreg2 `v' L(1/`p').(`vars_fd') i.year i.state [pw = emp] if year >= `sampstart', ///
        partial(i.year i.state)
        
        * Compute AIC for this model
        scalar bic = log(e(rss)/e(N)) + (log(e(N))/e(N)) * (`K' * `p' + 1)

        * Add aic and number of observations to the matrix
        matrix BIC[`p', `k'] = bic
        matrix BIC[ `maxlag' + 1 , `k'] = e(N)


    }

}

/******
Build Table
*******/
loc obs = e(N)

* Find minimum BIC per column, build substitute pairs
loc highlight_subs ""
if `highlight' {
    forv k = 1/`K' {
        scalar minval = BIC[1, `k']
        forv p = 2/`maxlag' {
            if BIC[`p', `k'] < minval scalar minval = BIC[`p', `k']
        }
        local mink = strtrim("`: display %9.3f minval'")
        local highlight_subs `"`highlight_subs' "`mink'" "\textcolor{mycolor2}{`mink'}""'
    }
}

esttab matrix(BIC, fmt(%9.3f)) using "${Tables}/LagSelect.tex", booktabs replace  ///
    nomtitles substitute("N " "\midrule N " "`obs'.000" "`obs'" `highlight_subs') ///
    addnotes("All regressions contain the indicated lags of all six variables shown in the column of this table.")

