/******************
REQUIRED PACKAGES: You will need to have xframappend installed to run this code.
> ssc install xframeappend
******************/

set mem 10g
clear all
do Globals.do
do Functions.do

use "$Data/CountyData.dta", clear

loc SampleDef01 "ForeignVar == 2 & HoursVar == 2" // Foreigners exclude naturalized citizens, Full time if >= 35 hrs last week
loc SampleDef02 "ForeignVar == 1 & HoursVar == 2" // Foreigners include naturalized citizens, Full time if >= 35 hrs last week
loc SampleDef03 "ForeignVar == 1 & HoursVar == 1" // Foreigners include naturalized citizens, Full time if >= 40 hrs last week
loc SampleDef04 "ForeignVar == 2 & HoursVar == 1" // Foreigners exclude naturalized citizens, Full time if >= 40 hrs last week

order fipcode year HoursVar ForeignVar foreign
sort fipcode year HoursVar ForeignVar foreign

/***********
A few touch ups
************/
gen Y = YNom * 100 / P
la var Y "Real GDP"
gen logY = log(Y)
gen StateAbb = substr(CountyName,-2,.)
replace StateAbb = substr(CountyName,-3,2) if substr(StateAbb,2,1) == "*"
bysort fipcode year HoursVar ForeignVar (foreign): gen obs = _N
drop if obs == 1
drop obs
encode StateAbb, gen(State)

* Creates separate variables for foreign and domestic quantities
* I think this can also be done with reshape...
foreach var in BodiesSupplied HoursSupplied {
    
    * Create separate variables for each value of foreign (0/1)
    separate `var', by(foreign)

    * Fill in the missing values (this should really be an option @ StataCorp)
    qui levelsof foreign
    foreach lvl in `r(levels)' {
        egen `var'`lvl'Fill = mean(`var'`lvl'), by(fipcode year HoursVar ForeignVar)
        replace `var'`lvl' = `var'`lvl'Fill / 1000 // GDP is in thousands of dollars
        drop `var'`lvl'Fill
    }
}

bysort fipcode year HoursVar ForeignVar (foreign): drop if [_n==2] // These are all redundant now

/**********
Estimation
***********/

* Save the number of distinct counties
encode fipcode, gen(fipcodenum)
qui tab fipcodenum
loc NCounties = `r(r)' - 1 // Leaving one out to avoid collinearity

* Save the number of years
qui tab year
qui tab year
loc NYears = `r(r)' - 1 // Leaving one out to avoid colllinearity

* Initialize and run
matrix InitVals = [1, 0.5, 1, 1, 0, J(1,`NYears',1), J(1,`NCounties', 1)]
loc k = colsof(InitVals)
nl CesFe @ logY BodiesSupplied1 BodiesSupplied0 if `SampleDef01', initial(InitVals) eps(1e-5) iter(1000) nparam(`k') ///
variables(BodiesSupplied0 BodiesSupplied1) hasconstant(b5)

* Calculate TFP and Save Data
frame copy default TfpEstimates, replace
frame TfpEstimates { // Save a new dataset with these estimates

    keep if `SampleDef01'

    predict double resid if `SampleDef01', res
    gen Z = exp(res + _b[/b5]) if fipcodenum == 1
    forvalues fip = 2/508 { // There are 508 fipcode nums but only 507 fe coeffecients
        loc feindex = 5 + 18 + `fip' - 1 // 5 CES params +  18 year fes and - 1 because fipcodenum==1 is left out of the estimation
        replace Z = exp(res + _b[/b`feindex'] + _b[/b5]) if fipcodenum == `fip'
    }

    keep fipcode year CountyName YNom HoursSupplied0 HoursSupplied1 BodiesSupplied0 ///
    BodiesSupplied1 met2013 P Y StateAbb Z
    save "$Data/CountyDataTfpEstimates.dta", replace

}
