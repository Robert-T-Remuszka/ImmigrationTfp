clear all
set mem 10g
do Globals
loc graphs "OlsOutput OlsTfp OlsCapital IvOlsTfp"
loc samp1 !inlist(statefip, "11")

use "$Data/StateAnalysisFileTfp.dta", clear

/*********************************************************************************************
Cleaning
**********************************************************************************************/
encode StateName, gen(State)
gen wt = BodiesSupplied1 + BodiesSupplied0 // Employment
order State year
xtset State year
gen f = (f1.BodiesSupplied1 - BodiesSupplied1)/wt
order State year
forvalues h = -2/10 { // LHS of LPs
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
    use "$Data/StateAnalysisFile.dta", clear
    replace foreign = 1 if ImmigrantGroup != "United States" & foreign == 0
    keep ImmigrantGroup foreign statefip year BodiesSupplied StateName
    egen L = total(BodiesSupplied), by(statefip year)
    gen share = BodiesSupplied / L
    keep if foreign
    encode ImmigrantGroup, gen(group)

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
    reshape wide share_, i(statefip year StateName) j(group)
}

* Retrieve the shares
frlink m:1 statefip year, frame(Shares)
frget share_*, from(Shares)
drop Shares
frame drop Shares

* Calculate Aggregate Shifts
frame create Shifts
frame Shifts {
    use "$Data/StateAnalysisFile.dta", clear
    replace foreign = 1 if ImmigrantGroup != "United States" & foreign == 0
    keep ImmigrantGroup foreign statefip year BodiesSupplied StateName
    keep if foreign
    encode ImmigrantGroup, gen(group)
    collapse (sum) National = BodiesSupplied, by(group year)
    xtset group year
    gen flow = f.National - National
    ren flow flow_
    ren National National_
    reshape wide National_ flow_, i(year) j(group)
}

* Retrieve the shifts
frlink m:1 year, frame(Shifts)
frget flow_* National_*, from(Shifts)
drop Shifts
frame drop Shifts

forvalues i = 1/10 {
    gen growth_`i' = flow_`i'/National_`i'
}

* Generate instruments
gen settlement = Africa_1960_share * growth_1 + CaAuNz_1960_share * growth_2 + China_1960_share * growth_3 ///
+ India_1960_share * growth_4 + LA_1960_share * growth_5 + Mexico_1960_share * growth_6 + Other_1960_share * growth_7 ///
+ AsiaOther_1960_share * growth_8 + EastEu_1960_share * growth_9 + WestEu_1960_share * growth_10

gen PredictedFlow = share_1 * growth_1 + share_2 * growth_2 + share_3 * growth_3 + share_4 * growth_4 + share_5 * growth_5 + ///
share_6 * growth_6 + share_7 * growth_7 + share_8 * growth_8 + share_9 * growth_9 + share_10 * growth_10

gen PredictedFlowl1 = l.share_1 * growth_1 + l.share_2 * growth_2 + l.share_3 * growth_3 + l.share_4 * growth_4 + l.share_5 * growth_5 + ///
l.share_6 * growth_6 + l.share_7 * growth_7 + l.share_8 * growth_8 + l.share_9 * growth_9 + l.share_10 * growth_10

gen PredictedFlowl2 = l2.share_1 * growth_1 + l2.share_2 * growth_2 + l2.share_3 * growth_3 + l2.share_4 * growth_4 + l2.share_5 * growth_5 + ///
l2.share_6 * growth_6 + l2.share_7 * growth_7 + l2.share_8 * growth_8 + l2.share_9 * growth_9 + l2.share_10 * growth_10
/****************************************************************************************
OLS
****************************************************************************************/
frame create Estimates
frame Estimates {
    gen rho = .
    gen h = .
    gen se = .
    gen rhoY = .
    gen seY = .
    gen rhoK = .
    gen seK = .
    gen rhoIv = .
    gen seIv = .
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
    
    qui reghdfe y`horizon' f [pw = wt], absorb(State year) vce(cluster year)
    frame Estimates {
        replace rhoY = _b[f] if _n == _N
        replace seY = _se[f] if _n == _N
    }
    qui reghdfe k`horizon' f [pw = wt], absorb(State year) vce(cluster year)
    frame Estimates {
        replace rhoK = _b[f] if _n == _N
        replace seK = _se[f] if _n == _N
    }
    
    cap ivregress 2sls z`horizon' (f = PredictedFlowl1) i.year i.State [pw=wt], vce(cluster year)
    if _rc != 0 {
        frame Estimates {
            replace rhoIv = 0 if _n == _N
            replace seIv = 0 if _n == _N
        }
    }
    else {
        estat firststage
        frame Estimates {
            replace rhoIv = _b[f] if _n == _N
            replace seIv = _se[f] if _n == _N
        }
    }
    
}

/***************************
OLS Impulse Plots
****************************/
frame Estimates {

    gen top = 1.645 * se + rho
    gen bottom = -1.645 * se + rho
    gen topY = 1.645 * seY + rhoY
    gen bottomY = -1.645 * se + rhoY
    gen topK = 1.645 * seK + rhoK
    gen bottomK = -1.645 * seK + rhoK
    gen topIv  = 1.645 * seIv + rhoIv
    gen bottomIv = -1.645 * seIv + rhoIv

    tw line rho h, lc(ebblue) || rarea top bottom h, fcolor(ebblue%30) lcolor(ebblue%50) lwidth(thin) ///
    legend(off) xlab(-2(1)10,nogrid) ylab(,nogrid) yline(0,lc(black%70) lp(solid)) ///
    name(OlsTfp) ytitle("{&beta}{subscript:h}") note("Error bands correspond to 90% level of confidence.")
    
    tw line rhoY h, lc(ebblue) || rarea topY bottomY h, fcolor(ebblue%30) lcolor(ebblue%50) lwidth(thin) ///
    legend(off) xlab(-2(1)10,nogrid) ylab(,nogrid) yline(0,lc(black%70) lp(solid)) ///
    name(OlsOutput) ytitle("{&beta}{subscript:h}") note("Error bands correspond to 90% level of confidence.")

    tw line rhoK h, lc(ebblue) || rarea topK bottomK h, fcolor(ebblue%30) lcolor(ebblue%50) lwidth(thin) ///
    legend(off) xlab(-2(1)10,nogrid) ylab(,nogrid) yline(0,lc(black%70) lp(solid)) ///
    name(OlsCapital) ytitle("{&beta}{subscript:h}") note("Error bands correspond to 90% level of confidence.")

    tw line rhoIv h, lc(ebblue) || rarea topIv bottomIv h, fcolor(ebblue%30) lcolor(ebblue%50) lwidth(thin) ///
    legend(off) xlab(-2(1)10,nogrid) ylab(,nogrid) yline(0,lc(black%70) lp(solid)) ///
    name(IvTfp) ytitle("{&beta}{subscript:h}") note("Error bands correspond to 90% level of confidence.")

    * Overlay the OLS and Iv estimates for Tfp
    tw line rho h, lc(ebblue) || line rhoIv h, lc(orange) || rarea top bottom h, fcolor(ebblue%30) lwidth(none) ///
    || rarea topIv bottomIv h, fcolor(orange%30) lwidth(none) ///
    legend(label(1 "OLS") label(2 "IV") order(1 2) pos(5) ring(0)) name(IvOlsTfp) xlab(-2(1)10,nogrid labsize(small)) ylab(,nogrid labsize(small)) ///
    yline(0,lc(black%70) lp(solid)) ytitle("{&beta}{subscript:h}") ///
    note("Standard errors clustered by year, 90% confidence. Shift-share IV constructed using j = 1.")
    
}

foreach g in `graphs' {
    graph export "$Graphs/`g'.pdf", replace name(`g')
}