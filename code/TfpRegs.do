clear all
do Globals

use "${Data}/StateAnalysis.dta", clear

loc samp inrange(year, 1994, 2021)
loc sampno21 inrange(year, 1994, 2020)

loc graphs "ZResponse_Iv1990 ZResponse_Iv1990_F"
/*************************************************
Some Cleaning
*************************************************/

* Total employment
gen emp = Supply_Total

* Generate foreign-born labor share and its growth rate
gen f = Supply_Foreign
bys statefip (year): gen fg = (f - f[_n - 1]) / emp

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

        egen Supply_Agg_`region' = total(`v'), by(year)
        bys statefip (year): gen fg_agg_`region' = (Supply_Agg_`region' - Supply_Agg_`region'[_n-1]) / emp_agg

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

* Calculate TFP Growth rates
forvalues h = -9/9 { // LHS variables
    if `h' < 0 loc name = "L" + string(abs(`h'))
    if `h' >= 0 loc name = "F" + string(abs(`h'))
    bys statefip (year): gen Zg_`name' = Z[_n + `h']/Z[_n - 1] - 1
}

ren state statename
encode statefip, gen(state)
xtset state year

/*************************************************
REGRESSIONS
*************************************************/

* Fixed-shares design
frame create Estimates
frame Estimates {

    gen h = .
    gen BetaIv1990 = .
    gen SeIv1990 = .
    gen FIv1990 = .

}

forvalues h = -6/9 {

    if `h' != -1 {
        di "***********************************************************************************************************"
        di "COMPUTING REGRESSION AT HORIZON `h'"
        di "***********************************************************************************************************"
        if `h' < 0 loc horizon = "L" + string(abs(`h'))
        if `h' >= 0 loc horizon = "F" + string(abs(`h'))

        qui ivreghdfe Zg_`horizon' (fg = Bartik_1990) l.Zg_`horizon' [pw = emp] if `samp', absorb(state year) vce(robust)

        * Record the results in the Estimates frame
        frame Estimates {
            insobs 1
            replace h = `h' if _n == _N
            replace BetaIv1990 = _b[fg] if _n == _N
            replace SeIv1990 = _se[fg] if _n == _N
            replace FIv1990 = `e(widstat)' if _n == _N
        }
    
    }

    * Fill in the horizon -1 results
    else {
        frame Estimates {
            insobs 1
            replace h = `h' if _n == _N
            replace BetaIv1990 = 0 if _n == _N
            replace SeIv1990 = . if _n == _N
            replace FIv1990 = . if _n == _N
        }
    }
    
}


/*************************************************
GRAPHS
*************************************************/

* IRF of pre-period share
frame Estimates {

    * Confidence intervals - 95%
    gen BetaIv1990_upper = BetaIv1990 + 1.96 * SeIv1990
    gen BetaIv1990_lower = BetaIv1990 - 1.96 * SeIv1990

    tw connected BetaIv1990 h if inrange(h,-6, 9), ms(oh) mc("0 147 245") xlab(-6(1)9, nogrid) sort || rcap BetaIv1990_upper BetaIv1990_lower h, lcolor("0 147 245") ylab(, nogrid) ///
    ytitle("{&eta}{subscript:Z}") xtitle("Horizon") legend(off) yline(0, lc(black%50) lp(solid)) name(ZResponse_Iv1990)
    rtr

}

* F Stats of IRF for pre-period
frame Estimates {

    graph bar FIv1990 if h > -1, over(h) bar(1, color("0 147 245") fcolor("0 147 245")) ylab(0(10)40, nogrid labsize(small)) ///
    yline(10,lc(black%70) lp(dash)) legend(off) b1title("Horizon") ytitle("First Stage F") name(ZResponse_Iv1990_F)
    
}

* Export
foreach g in `graphs' {

    graph export "${Graphs}/`g'.pdf", name(`g') replace

}
