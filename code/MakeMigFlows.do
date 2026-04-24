clear all
do Globals

/*********
Creating a dataset that reports migration flows between US and ROW
**********/
frame create appended
frame create outflows_appended

****** Appending the ACS files first
loc files: dir "${Data}/acs" files "*.dta"
foreach f in `files' {

    if substr("`f'", 4, 4) != "2000" { // Use the CPS for this year

        use "${Data}/acs/`f'", clear

        **** Apply the same sample restrictions as in the previous data
        * Sample restrictions
        drop if AGE < 16 | AGE == 999                               // AGE == 999 if missing
        drop if UHRSWORK < 35                                       // UHRSWORK < 35 also gets rid of missings which are coded as 0

        * Prepare the birthplace variable
        tostring BPLD, replace
        replace BPLD = "00" + BPLD if strlen(BPLD) == 3
        replace BPLD = "0" + BPLD  if strlen(BPLD) == 4
        drop if inlist(substr(BPLD,1,1), "8", "9")                  // Drop those who don't have a country of origin

        * Some additional cleaning
        drop if CITIZEN == 9
        replace INCWAGE = . if inlist(INCWAGE, 999999, 999998)      // Drop missing wage data

        * Assigning origins based on Peri (2012)
        gen ImmigrantGroup = "United States" if substr(BPLD,1,1) == "0"
        replace ImmigrantGroup = "Mexico" if BPLD == "20000"
        replace ImmigrantGroup = "Latin America" if BPLD == "11000" | (inlist(substr(BPLD,1,2), "21", "25", "26", "30") & BPLD != "26030")
        replace ImmigrantGroup = "Western Europe" if inlist(BPLD,"45300","45000") | inlist(substr(BPLD,1,2), "40", "41", "42", "43")
        replace ImmigrantGroup = "Russia and Eastern Europe" if inlist(substr(BPLD,1,2), "45", "46") & !inlist(BPLD,"45300", "45000")
        replace ImmigrantGroup = "Canada-Australia-New Zealand" if inlist(BPLD, "15000", "70020", "70010") | substr(BPLD, 1, 2) == "15"
        replace ImmigrantGroup = "China" if BPLD == "50000"
        replace ImmigrantGroup = "India" if BPLD == "52100"
        replace ImmigrantGroup = "Rest of Asia" if inlist(substr(BPLD,1,2), "50", "51", "52","53", "54", "55") & !inlist(BPLD,"52100","50000")
        replace ImmigrantGroup = "Africa" if substr(BPLD,1,2) == "60"
        drop if mi(ImmigrantGroup)                                 // Drop unassigned origins - primarily US terriories

        * Create a domestic born indicator
        gen Domestic = ImmigrantGroup == "United States"
        drop ImmigrantGroup

        * Keeping what we need
        keep PERWT STATEFIP YEAR Domestic MIGPLAC1
        drop if inlist(MIGPLAC1, 999, 997, 988)

        * Set up regions - in the US last year then extract the state fip, else 'ROW'
        tostring MIGPLAC1, replace
        tostring STATEFIP, replace
        replace MIGPLAC1 = "00" + MIGPLAC1 if strlen(MIGPLAC1) == 1
        replace MIGPLAC1 = "0"  + MIGPLAC1 if strlen(MIGPLAC1) == 2
        replace STATEFIP = "0"  + STATEFIP if strlen(STATEFIP) == 1

        gen Origin = "ROW" if substr(MIGPLAC1, 1, 1) != "0"
        replace Origin = substr(MIGPLAC, 2, 2) if mi(Origin)
        replace Origin = STATEFIP if Origin == "00"
        drop MIGPLAC1
        
        * Count flows by origin and destination using PERWT
        collapse (sum) PERWT, by(YEAR STATEFIP Origin Domestic)
        ren PERWT Flow
        ren YEAR year
        ren STATEFIP statefip
        reshape wide Flow, i(year statefip Origin) j(Domestic)
        ren Flow1 Domestic_
        ren Flow0 Foreign_
        replace Foreign = 0 if mi(Foreign)
        replace Domestic = 0 if mi(Domestic)
        fillin statefip Origin
        drop _fillin
        qui summ year
        replace year = `r(mean)' if mi(year)
        replace Foreign = 0 if mi(Foreign)
        replace Domestic = 0 if mi(Domestic)
        
        * Append
        frame appended: xframeappend default
    }
}

****** Now append the CPS files
loc files: dir "${Data}/cps/" files "*.dta"
foreach f in `files' {

    use "${Data}/cps/`f'", clear

    * Setup statefip codes as strings
    tostring STATEFIP, replace
    replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1
    
    * Sample restrictions
    drop if AGE < 16 | AGE == 999                               // AGE == 999 if missing
    drop if UHRSWORK < 35                                       // UHRSWORK < 35 also gets rid of missings which are coded as 0

    * Prepare the birthplace variable
    tostring BPL, replace
    replace BPL = "00" + BPL if strlen(BPL) == 3
    replace BPL = "0" + BPL  if strlen(BPL) == 4
    drop if inlist(substr(BPL,1,1), "8", "9")                  // Drop those who don't have a country of origin

    * Some additional cleaning
    drop if CITIZEN == 9
    replace INCWAGE = . if inlist(INCWAGE, 999999, 999998)          // Drop missing wage data
    egen incwage_censor = pctile(INCWAGE), p(99) by(STATEFIP)       // Make the censoring as consistent as possible with ACS
    replace INCWAGE = incwage_censor if INCWAGE >= incwage_censor


    * Assigning origins based on Peri (2012)
    gen ImmigrantGroup = "United States" if substr(BPL,1,1) == "0"
    replace ImmigrantGroup = "Mexico" if BPL == "20000"
    replace ImmigrantGroup = "Latin America" if BPL == "11000" | (inlist(substr(BPL,1,2), "21", "25", "26", "30") & BPL != "26030")
    replace ImmigrantGroup = "Western Europe" if inlist(BPL,"45300","45000") | inlist(substr(BPL,1,2), "40", "41", "42", "43")
    replace ImmigrantGroup = "Russia and Eastern Europe" if inlist(substr(BPL,1,2), "45", "46") & !inlist(BPL,"45300", "45000")
    replace ImmigrantGroup = "Canada-Australia-New Zealand" if inlist(BPL, "15000", "70020", "70010") | substr(BPL, 1, 2) == "15"
    replace ImmigrantGroup = "China" if BPL == "50000"
    replace ImmigrantGroup = "India" if BPL == "52100"
    replace ImmigrantGroup = "Rest of Asia" if inlist(substr(BPL,1,2), "50", "51", "52","53", "54", "55") & !inlist(BPL,"52100","50000")
    replace ImmigrantGroup = "Africa" if substr(BPL,1,2) == "60"
    drop if mi(ImmigrantGroup)                                 // Drop unassigned origins - primarily US terriories

    * Create a domestic born indicator
    gen Domestic = ImmigrantGroup == "United States"
    drop ImmigrantGroup

    * Keeping what we need
    keep ASECWT STATEFIP YEAR Domestic MIGSTA1

    * Set up regions
    tostring MIGSTA1, replace
    tostring STATEFIP, replace
    replace MIGSTA1 = "0" + MIGSTA1 if strlen(MIGSTA1) == 1
    replace STATEFIP = "0"  + STATEFIP if strlen(STATEFIP) == 1

    * generate origins
    gen Origin = "ROW" if MIGSTA1 == "91"
    replace Origin = STATEFIP if MIGSTA1 == "99"
    replace Origin = MIGSTA1  if mi(Origin)
    drop if MIGSTA1 == "00"
    drop MIGSTA1

    * Calculate the flows
    collapse (sum) ASECWT, by(YEAR STATEFIP Origin Domestic)
    ren ASECWT Flow
    reshape wide Flow, i(YEAR STATEFIP Origin) j(Domestic)
    fillin STATEFIP Origin
    drop _fillin
    ren YEAR year
    ren STATEFIP statefip
    ren Flow0 Foreign_
    ren Flow1 Domestic_
    qui summ year
    replace year = `r(mean)' if mi(year)
    replace Foreign = 1e-4 if mi(Foreign)
    replace Domestic = 1e-4 if mi(Domestic)

    frame appended: xframeappend default
    
}

*** Calculate flows from each state to ROW. The following identity must hold for each state
* dStock = Total Inflows - Outflows <--> 
*        = Total Inflows - (Outflows to rest of US + Outflow to ROW) <-->
* Outflow to ROW = Total Inflows - OutFlows to rest of US - dStock

* Records inflows from each state and ROW into any US region
frame copy appended default, replace

* Total inflows to each state (excluding stayers)
preserve
    keep if Origin != statefip
    collapse (sum) Foreign_ Domestic_, by(year statefip)
    ren Foreign_ TotalIn_For
    ren Domestic_ TotalIn_Dom
    tempfile inflows
    save `inflows'
restore

* Total outflows from each state to other US states
preserve
    keep if Origin != "ROW" & Origin != statefip
    collapse (sum) Foreign_ Domestic_, by(year Origin)
    ren Origin statefip
    ren Foreign_ TotalOut_For
    ren Domestic_ TotalOut_Dom
    tempfile outflows
    save `outflows'
restore

* Merge stock changes with inflows and outflows; compute ROW outflows
frame create row_calc
frame row_calc {
    use "${Data}/StateAnalysis.dta", clear
    keep statefip year Supply_Foreign Supply_Domestic
    bys statefip (year): gen dFor = Supply_Foreign - Supply_Foreign[_n-1]
    bys statefip (year): gen dDom = Supply_Domestic - Supply_Domestic[_n-1]
    keep statefip year dFor dDom
    drop if mi(dFor) | mi(dDom)
    merge 1:1 year statefip using `inflows', nogen keep(3)
    merge 1:1 year statefip using `outflows', nogen keep(3)
    gen Foreign_  = max(1e-4, TotalIn_For - TotalOut_For - dFor)
    gen Domestic_ = max(1e-4, TotalIn_Dom - TotalOut_Dom - dDom)
    keep year statefip Foreign_ Domestic_
    ren statefip Origin
    gen statefip = "ROW"
}

* use the stock-flow relation to calculate flows from ROW to ROW - need some UN pop counts for this
frame create row_to_row
frame row_to_row {
    import delimited "${UnPopCounts}/WPP2024_Demographic_Indicators_Medium.csv", clear

    keep if loctypename == "World" | location == "United States of America"
    keep tpopulation1jan time location
    keep if time > 1993
    replace location = "US" if location == "United States of America"
    ren tpopulation pop_
    ren time year
    reshape wide pop_, i(year) j(location) s
    gen pop_row = pop_World - pop_US
    replace pop_row = pop_row * 1000        // UN data is in thousands
    drop pop_World pop_US
    sort year
    gen pop_row_L = pop_row[_n-1]
}

* ROW→US Foreign flows by year
preserve
    keep if Origin == "ROW"
    collapse (sum) Foreign_, by(year)
    ren Foreign_ ROW_to_US_For
    tempfile row_to_us
    save `row_to_us'
restore

frame row_to_row {
    merge 1:1 year using `row_to_us', nogen keep(3)
    gen Foreign_  = max(1e-4, pop_row_L - ROW_to_US_For)   // Those of the ROW stock that didn't come to US must have stayed
    gen Domestic_ = 1e-4                                   // This is an approximation, this is going to be a small number
    gen statefip  = "ROW"
    gen Origin    = "ROW"
    keep year statefip Origin Foreign_ Domestic_
}

frame appended {
    drop if year == 1994 // don't have ROW inflows in this year becuase of the stock approach taken above
    xframeappend row_calc, drop
    xframeappend row_to_row, drop
}

*** Winsorizing and saving
frame copy appended default, replace
sort year statefip

qui summ Foreign_ if !(statefip == "ROW" & Origin == "ROW"), d
replace Foreign_  = min(`r(p95)', Foreign_) if !(statefip == "ROW" & Origin == "ROW")
qui summ Domestic_ if !(statefip == "ROW" & Origin == "ROW"), d
replace Domestic_ = min(`r(p95)', Domestic_) if !(statefip == "ROW" & Origin == "ROW")

reshape wide For Dom, i(year statefip) j(Origin) s
drop if year == 1994 // don't need this one observation

save "${Data}/MigFlows.dta", replace