clear all
do Globals
loc graphs "IvOlsTfp IvLooOlsTfp FirstStageF"
loc samp1 !inlist(statefip, "11")
loc samp2 year >= 1995

use "$Data/StateAnalysisFileTfp.dta", clear

/*********************************************************************************************
Cleaning
**********************************************************************************************/
encode StateName, gen(State)
gen wt = BodiesSupplied1 + BodiesSupplied0           // Total Employment
order State year
xtset State year
gen f = (f1.BodiesSupplied1 - BodiesSupplied1)/wt    // Create (endoegenous) migration flow
order State year
forvalues h = -2/10 { // LHS variables
    if `h' < 0 loc name = "Back" + string(abs(`h'))
    if `h' >= 0 loc name = "Forward" + string(abs(`h'))
    bysort State (year): gen z`name' = Z[_n + `h']/Z[_n] - 1
    bysort State (year): gen y`name' = Y[_n + `h']/Y[_n] - 1
    bysort State (year): gen k`name' = K[_n + `h']/K[_n] - 1
}

* Past settlement instrument shares - I will also calculate more recent shares in another frame
foreach v of varlist *_1960 {
    gen `v'_share = `v'/wt
}

* Calculate Immigrant group shares
frame create Shares
frame Shares {
    use "$Data/StateAnalysisFile.dta", clear          // Note that this loads a file which contains tabluation across migrant groups
    replace foreign = 1 if ImmigrantGroup != "United States" & foreign == 0
    keep ImmigrantGroup foreign statefip year BodiesSupplied StateName
    egen L = total(BodiesSupplied), by(statefip year) // Create the labor force size (wt above)
    gen share = BodiesSupplied / L
    keep if foreign
    encode ImmigrantGroup, gen(group)
    levelsof group, loc(groups)
    qui summ group
    loc maxgroup = `r(max)' // useful for later (is just 10)
    /*
    The encode key is as follows
    1 = Africa
    2 = Canada-Australia-New-Zealand
    3 = China
    4 = India
    5 = Latin America
    6 = Mexico
    7 = Other
    8 = Rest of Asia
    9 = Russia and Eastern Europe
    10 = Western Europe
    */
    ren share share_
    keep share_ group statefip year StateName
    duplicates drop statefip year StateName group, force
    reshape wide share_, i(statefip year StateName) j(group) // Obervation is now (state, year)

    foreach grp of numlist `groups' { // Create 1994 shares
        bysort StateName (year): gen share1994_`grp' = share_`grp'[1]
        replace share1994_`grp' = 0 if mi(share1994_`grp')
    }
}

* Retrieve the shares
frlink m:1 statefip year, frame(Shares)
frget share_* share1994_*, from(Shares)
drop Shares
frame drop Shares

* Calculate Aggregate Shifts
frame create Shifts
frame Shifts {
    use "$Data/StateAnalysisFile.dta", clear // An observation in this dataset is (state, year, migrant group)
    replace foreign = 1 if ImmigrantGroup != "United States" & foreign == 0
    keep ImmigrantGroup foreign statefip year BodiesSupplied StateName
    keep if foreign
    encode ImmigrantGroup, gen(group)
    qui levelsof group, loc(groups)
    encode StateName, gen(StateGroup)
    qui levelsof StateGroup, loc(states)
    foreach state in `states' { // For the leave one out totals
        gen BodiesSupplied_`state' = BodiesSupplied if StateGroup == `state'
    }
    collapse (sum) National = BodiesSupplied (firstnm) BodiesSupplied_*, by(group year)
    xtset group year
    gen flow = f.National - National
    foreach state in `states' { // Leave one out
        gen NationalLess`state'_ = National - BodiesSupplied_`state'
        gen FlowLess`state'_ = f.NationalLess`state' - NationalLess`state'
    }

    ren flow flow_
    ren National National_
    drop BodiesSupplied_*
    reshape wide National_ flow_ FlowLess* NationalLess*, i(year) j(group)
    foreach i of numlist `groups' { // Create national growth rates
        gen growth_`i' = flow_`i'/National_`i'
        foreach state of numlist `states' { // Leave one out
            gen GrowthLess`state'_`i' = FlowLess`state'_`i' / NationalLess`state'_`i'
        }
        egen growth_loo_`i' = rowmean(GrowthLess*)
    }
}

* Retrieve the shifts
frlink m:1 year, frame(Shifts)
frget growth_* growth_loo_*, from(Shifts)
drop Shifts
frame drop Shifts

* Generate instruments
gen settlement = Africa_1960_share * growth_1 + CaAuNz_1960_share * growth_2 + China_1960_share * growth_3 ///
+ India_1960_share * growth_4 + LA_1960_share * growth_5 + Mexico_1960_share * growth_6 + Other_1960_share * growth_7 ///
+ AsiaOther_1960_share * growth_8 + EastEu_1960_share * growth_9 + WestEu_1960_share * growth_10

gen PredictedFlow = share_1 * growth_1 + share_2 * growth_2 + share_3 * growth_3 + share_4 * growth_4 + share_5 * growth_5 + ///
share_6 * growth_6 + share_7 * growth_7 + share_8 * growth_8 + share_9 * growth_9 + share_10 * growth_10

* Average growth rates and first lag share
gen PredictedFlowl1 = l.share_1 * growth_1 + l.share_2 * growth_2 + l.share_3 * growth_3 + l.share_4 * growth_4 + l.share_5 * growth_5 + ///
l.share_6 * growth_6 + l.share_7 * growth_7 + l.share_8 * growth_8 + l.share_9 * growth_9 + l.share_10 * growth_10

* Leave one out growth rates and first lag share
gen PredictedLooFlow1 = l.share_1 * growth_loo_1 + l.share_2 * growth_loo_2 + l.share_3 * growth_loo_3 + l.share_4 * growth_loo_4 + l.share_5 * growth_loo_5 + ///
l.share_6 * growth_loo_6 + l.share_7 * growth_loo_7 + l.share_8 * growth_loo_8 + l.share_9 * growth_loo_9 + l.share_10 * growth_loo_10

* Leave one out and fixed 1994 shares
gen BartikLoo1994Shares = share1994_1 * growth_loo_1 + share1994_2 * growth_loo_2 + share1994_3 * growth_loo_3 + share1994_4 * growth_loo_4 + ///
share1994_5 * growth_loo_5 + share1994_6 * growth_loo_6 + share1994_7 * growth_loo_7 + share1994_8 * growth_loo_8 + share1994_9 * growth_loo_9 + ///
share1994_10 * growth_loo_10


/****************************************************************************************
OLS
****************************************************************************************/
frame create Estimates
frame Estimates {

    * Horizon of each regression
    gen h = .

    * Ols estimates
    gen rho = .
    gen se = .

    * Iv Estimates
    gen rhoIv = .
    gen seIv = .
    gen FStatIv = .

    * Iv with leave one out
    gen rhoIvLoo = .
    gen seIvLoo = .
    gen FStatIvLoo = .

}

forvalues h = -2/10 {

    di "***********************************************************************************************************"
    di "COMPUTING REGRESSION AT HORIZON `h'"
    di "***********************************************************************************************************"
    if `h' < 0 loc horizon = "Back" + string(abs(`h'))
    if `h' >= 0 loc horizon = "Forward" + string(abs(`h'))

    qui reghdfe z`horizon' f [pw = wt], absorb(State year) vce(cluster year)
    frame Estimates {
        insobs 1
        replace h = `h' if _n == _N
        replace rho = _b[f] if _n == _N
        replace se = _se[f] if _n == _N
    }

    if `h' != 0 { 
        
        * National growth rate IV regs
        qui {
            ivregress 2sls z`horizon' (f = PredictedFlowl1) i.year i.State [pw=wt], vce(cluster year)
            estat firststage
            mat Stats = r(singleresults)
        }
        
        frame Estimates {
            replace FStatIv = Stats[1,4] if _n == _N
            replace rhoIv = _b[f] if _n == _N
            replace seIv = _se[f] if _n == _N
        }
        
        * Leave one out growth rates
        qui {
            ivregress 2sls z`horizon' (f = PredictedLooFlow1) i.year i.State [pw=wt], vce(cluster year)
            estat firststage
            mat Stats = r(singleresults)
        }
        
        frame Estimates {
            replace FStatIvLoo = Stats[1,4] if _n == _N
            replace rhoIvLoo = _b[f] if _n == _N
            replace seIvLoo = _se[f] if _n == _N
        }
        
    }

    else { // Horizon zero won't compute using ivregress
        frame Estimates {
            replace FStatIv = 0 if _n == _N
            replace rhoIv = 0 if _n == _N
            replace seIv = 0 if _n == _N
            replace FStatIvLoo = 0 if _n == _N
            replace rhoIvLoo = 0 if _n == _N
            replace seIvLoo = 0 if _n == _N
        }
    }
    //cap ivregress 2sls z`horizon' (f = PredictedFlowl1) i.year i.State [pw=wt], vce(cluster year)
    //estat firststage
    //loc FStat = r(singleresults)[1,4]
    /*if _rc != 0 {
        frame Estimates {
            replace rhoIv = 0 if _n == _N
            replace seIv = 0 if _n == _N
        }
    }
    else {
        frame Estimates {
            replace FStatIv = `FStat'
            replace rhoIv = _b[f] if _n == _N
            replace seIv = _se[f] if _n == _N
        }
    }

    cap ivregress 2sls z`horizon' (f = PredictedLooFlow1) i.year i.State [pw=wt], vce(cluster year)
    estat firststage
    loc FStat = r(singleresults)[1,4]
    if _rc != 0 {
        frame Estimates {
            replace rhoIvLoo = 0 if _n == _N
            replace seIvLoo = 0 if _n == _N
        }
    }
    else {
        frame Estimates {
            replace FStatIvLoo = `FStat'
            replace rhoIvLoo = _b[f] if _n == _N
            replace seIvLoo = _se[f] if _n == _N
        }
    }
    */ 
}

/***************************
Impulse Responses
****************************/
frame Estimates {

    gen top = 1.645 * se + rho
    gen bottom = -1.645 * se + rho
    gen topIv  = 1.645 * seIv + rhoIv
    gen bottomIv = -1.645 * seIv + rhoIv
    gen topIvLoo = 1.645 * seIvLoo + rhoIvLoo
    gen bottomIvLoo = - 1.645 * seIvLoo + rhoIvLoo

    * Overlay the OLS and Iv estimates for Tfp
    tw line rho h, lc(ebblue) || line rhoIv h, lc(orange) || rarea top bottom h, fcolor(ebblue%30) lwidth(none) ///
    || rarea topIv bottomIv h, fcolor(orange%30) lwidth(none) ///
    legend(label(1 "OLS") label(2 "IV") order(1 2) pos(5) ring(0)) name(IvOlsTfp) xlab(-2(1)10,nogrid labsize(small)) ylab(,nogrid labsize(small)) ///
    yline(0,lc(black%70) lp(solid)) ytitle("{&beta}{subscript:h}") ///
    note("Standard errors clustered by year, 90% confidence. Shift-share IV constructed using j = 1.")

    * Overlay the OLS and Iv estimates for Tfp using Loo flows
    tw line rho h, lc(ebblue) || line rhoIvLoo h, lc(orange) || rarea top bottom h, fcolor(ebblue%30) lwidth(none) ///
    || rarea topIvLoo bottomIvLoo h, fcolor(orange%30) lwidth(none) ///
    legend(label(1 "OLS") label(2 "IV") order(1 2) pos(5) ring(0)) name(IvLooOlsTfp) xlab(-2(1)10,nogrid labsize(small)) ylab(,nogrid labsize(small)) ///
    yline(0,lc(black%70) lp(solid)) ytitle("{&beta}{subscript:h}")
}

/***************************
Robustness Checks and Testing
****************************/
frame Estimates {
    graph bar FStatIv FStatIvLoo if h > 0, over(h) bar(1, color(ebblue%70) fcolor(ebblue%70)) bar(2, color(orange%70) fcolor(orange%70)) ///
    ylab(0(20)200,nogrid labsize(small)) yline(10,lc(black%70) lp(dash)) ///
    legend(label(1 "Baseline Iv") label(2 "Leave One Out Iv") order(1 2) pos(1) ring(0)) b1title("Horizon") ytitle("First Stage F-Stat") ///
    name(FirstStageF)
}
/**************************
Export Graphs
**************************/
foreach g in `graphs' {
    graph export "$Graphs/`g'.pdf", replace name(`g')
}