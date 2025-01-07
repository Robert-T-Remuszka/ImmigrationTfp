/******************
REQUIRED PACKAGES: You will need to have xframappend installed to run this code.
> ssc install xframeappend
******************/

set mem 10g
clear all
do Globals.do

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
        replace `var'`lvl' = `var'`lvl'Fill
        drop `var'`lvl'Fill
    }
}

bysort fipcode year HoursVar ForeignVar (foreign): drop if [_n==2] // These are all redundant now

/**********
Estimation
***********/
* Create fixed effects local for estimation
qui tab State, gen(I)
loc FEs ""
foreach s of numlist 1/`r(r)' {
    loc FEs "`FEs' I`s'"
}

* Initialize data frame for estimates
frame create TfpEstimates
frame TfpEstimates {
    gen State = .
    gen year = .
    gen Z = .
}

foreach year of numlist 2005(1)2023 {
    
    * Estimate TFP by Nonlinear Least Squares
    nl (logY = {fe: `FEs'} + {beta} * log( ///
    {lambda}^(1 - 1/{beta})*({alphaF}*BodiesSupplied1)^(1/{beta}) + (1-{lambda})^(1-1/{beta})*({alphaD}*BodiesSupplied0)^(1/{beta}) ///
    )) if year == `year' & `SampleDef01', initial(beta 1 alphaF 1 alphaD 1 lambda 0.5) nocons eps(1e-5) iter(1000)

    foreach state of numlist 1/50 {
        
        frame create AppendThis
        frame AppendThis { // Store the estimates
            insobs 1
            gen State = `state'
            gen year = `year'
            gen Z = exp(_b[/fe_I`state'])
            replace Z = . if Z == 1
        }

        frame TfpEstimates: xframeappend AppendThis, drop
    }
    
}

frame TfpEstimates {
    use "$Data/TfpByState.dta", clear
    la def State ///
        1 "AK" ///
        2 "AL" ///
        3 "AR" ///
        4 "AZ" ///
        5 "CA" ///
        6 "CO" ///
        7 "CT" ///
        8 "DC" ///
        9 "DE" /// 
        10 "FL" ///
        11 "GA" ///
        12 "HI" ///
        13 "IA" ///
        14 "ID" ///
        15 "IL" ///
        16 "IN" ///
        17 "KS" ///
        18 "KY" ///
        19 "LA" ///
        20 "MA" ///
        21 "MD" ///
        22 "ME" ///
        23 "MI" ///
        24 "MN" ///
        25 "MO" ///
        26 "MS" ///
        27 "MT" ///
        28 "NC" ///
        29 "ND" ///
        30 "NE" ///
        31 "NH" ///
        32 "NJ" ///
        33 "NM" ///
        34 "NV" ///
        35 "NY" /// 
        36 "OH" ///
        37 "OK" ///
        38 "OR" ///
        39 "PA" ///
        40 "RI" ///
        41 "SC" ///
        42 "TN" ///
        43 "TX" ///
        44 "UT" ///
        45 "VA" ///
        46 "VT" ///
        47 "WA" ///
        48 "WI" ///
        49 "WV" ///
        50 "WY"

    la values State State
    decode State, gen(StateAbb)
    drop State
    order StateAbb year
    sort StateAbb year

    save "$Data/TfpByState.dta", replace
    
}