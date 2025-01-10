set mem 10g
clear all
do Globals.do
do Functions.do

use "$Data/StateAnalysisFile.dta", clear
order statefip StateAbb year ImmigrantGroup foreign
sort statefip StateAbb year ImmigrantGroup foreign

/***********
A few touch ups
************/
* These were likely missing CITIZEN, or born in US territories other than PR
replace foreign = 1 if ImmigrantGroup != "United States" & foreign == 0 

* Calculate real GDP and log
gen Y = NGdp * 100 / P
la var Y "Real GDP"
gen logY = log(Y)

* I converted the units to millions of hours in Python
replace HoursSupplied = HoursSupplied * 1e+6

/***************************
TFP Estimation
****************************/
frame copy default TfpEstimation
frame TfpEstimation {
    
    collapse (sum) HoursSupplied BodiesSupplied ///
    (firstnm) NGdp P *_1920 *_1930 *_1940 *_1950 *_1960 Y logY, by(statefip StateAbb year foreign)

    * Create separate variables for each value of foreign (0/1)
    foreach var in BodiesSupplied HoursSupplied {
    
        separate `var', by(foreign)

        * Fill in the missing values (this should really be an option @ StataCorp)
        qui levelsof foreign
        foreach lvl in `r(levels)' {
            egen `var'`lvl'Fill = mean(`var'`lvl'), by(statefip year)
            replace `var'`lvl' = `var'`lvl'Fill / 1e+3 // GDP is in thousands of dollars
            drop `var'`lvl'Fill
        }
    }

    * Every second observation is redundant now and so are the non-separated variables
    bysort statefip year (foreign): drop if [_n==2]
    drop foreign HoursSupplied BodiesSupplied

    * For state fixed effects
    encode statefip, gen(state)
    
    * How many paramters are there excluding the FEs?
    qui tab state
    glo NRegions = `r(r)'
    glo M = 5                            // Constant term + 4 CES parameters
    glo FirstYear = 1994
    glo LastYear = 2023
    loc NYearFe = $LastYear - $FirstYear // Number of time FEs in the model
    loc NStateFe = $NRegions - 1         // Number of State FEs

    matrix InitVals = [1, 0.5, 1, 1, 0, J(1,`NYearFe',0), J(1,`NStateFe', 0)]
    loc k = colsof(InitVals)
    nl CesFe @ logY BodiesSupplied1 BodiesSupplied0, initial(InitVals) eps(1e-5) iter(1000) nparam(`k') ///
    variables(BodiesSupplied0 BodiesSupplied1) hasconstant(b5)

    predict double resid, res
    gen Z = exp(res + _b[/b5]) if state == 1
    forvalues fip = 2/$NRegions {
        loc feindex = $M + `NYearFe' + `fip' - 1
        replace Z = exp(res + _b[/b`feindex'] + _b[/b5]) if state == `fip'
    }

    drop resid state
    save "$Data/StateAnalysisFileTfp.dta", replace
}
