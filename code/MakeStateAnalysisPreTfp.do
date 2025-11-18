clear all
do Globals

loc skipIpums     0
loc skipBea       0
loc skipEY        0
loc skipPrePeriod 0
loc skipMerge     0

/*********************************
    CLEAN AND APPEND IPUMS DATA
**********************************/
if !`skipIpums' {

    /*********************************
            CLEAN ACS
    **********************************/
    loc files: dir "${Data}/acs/" files "*.dta"
    foreach f in `files' {
        
        frame create appendthis
        frame appendthis {

            use "${Data}/acs/`f'", clear
            
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


            * Keep just what we need for aggregated analysis
            keep PERWT STATEFIP YEAR ImmigrantGroup INCWAGE
            qui ds ImmigrantGroup, not
            foreach v in `r(varlist)' {
                loc newname = strlower("`v'")
                ren `v' `newname'
            }

            * Aggregate up to the state x year x immigrant group level
            collapse (mean) Wage = incwage (count) Supply = incwage [fw = perwt], by(statefip year ImmigrantGroup)

            * Reshape wide - some groups have short enough names that we can use those as their stubs
            gen abbr = "_CaAuNz" if ImmigrantGroup == "Canada-Australia-New Zealand"
            replace abbr = "_US" if ImmigrantGroup == "United States"
            replace abbr = "_WestEu" if ImmigrantGroup == "Western Europe"
            replace abbr = "_LA" if ImmigrantGroup == "Latin America"
            replace abbr = "_EastEu" if ImmigrantGroup == "Russia and Eastern Europe"
            replace abbr = "_AsiaOther" if ImmigrantGroup == "Rest of Asia"
            replace abbr = "_" + ImmigrantGroup if mi(abbr)

            tostring year, replace
            tostring statefip, replace
            replace statefip = "0" + statefip if strlen(statefip) == 1
            drop ImmigrantGroup
            reshape wide Wage Supply, i(statefip year) j(abbr) s

            * Labor wage and quantity variables
            qui ds Wage_*
            foreach v in `r(varlist)' {
                
                * Put longer group name in variable labels
                loc origin = substr("`v'", 6, .)
                if substr("`v'", 6, .) == "CaAuNz" loc origin "Canada-Australia-New Zealand"
                if substr("`v'", 6, .) == "US" loc origin "United States"
                if substr("`v'", 6, .) == "WestEu" loc origin "Western Europe"
                if substr("`v'", 6, .) == "LA" loc origin "Latin America"
                if substr("`v'", 6, .) == "EastEu" loc origin "Russia and Eastern Europe"
                if substr("`v'", 6, .) == "AsiaOther" loc origin "Asia Other"
                
                la var `v' "Migrant wage (2017 dollars), `origin'"
            }

            qui ds Supply_*
            foreach v in `r(varlist)' {

                loc origin = substr("`v'", 8, .)
                if substr("`v'", 8, .) == "CaAuNz" loc origin "Canada-Australia-New Zealand"
                if substr("`v'", 8, .) == "US" loc origin "United States"
                if substr("`v'", 8, .) == "WestEu" loc origin "Western Europe"
                if substr("`v'", 8, .) == "LA" loc origin "Latin America"
                if substr("`v'", 8, .) == "EastEu" loc origin "Russia and Eastern Europe"
                if substr("`v'", 8, .) == "AsiaOther" loc origin "Asia Other"

                la var `v' "Foreign-born quantity, `origin'"
            }
        }

        xframeappend appendthis, drop
    }
    /*********************************
            CLEAN CPS
    **********************************/
    loc files: dir "${Data}/cps/" files "*.dta"
    foreach f in `files' {
        
        frame create appendthis
        frame appendthis {

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


            * Keep just what we need for aggregated analysis
            keep ASECWT STATEFIP YEAR ImmigrantGroup INCWAGE
            qui ds ImmigrantGroup, not
            foreach v in `r(varlist)' {
                loc newname = strlower("`v'")
                ren `v' `newname'
            }

            * Aggregate up to the state x year x immigrant group level
            collapse (mean) Wage = incwage (count) Supply = incwage [pw = asecwt], by(statefip year ImmigrantGroup)

            * Reshape wide - some groups have short enough names that we can use those as their stubs
            gen abbr = "_CaAuNz" if ImmigrantGroup == "Canada-Australia-New Zealand"
            replace abbr = "_US" if ImmigrantGroup == "United States"
            replace abbr = "_WestEu" if ImmigrantGroup == "Western Europe"
            replace abbr = "_LA" if ImmigrantGroup == "Latin America"
            replace abbr = "_EastEu" if ImmigrantGroup == "Russia and Eastern Europe"
            replace abbr = "_AsiaOther" if ImmigrantGroup == "Rest of Asia"
            replace abbr = "_" + ImmigrantGroup if mi(abbr)

            tostring year, replace
            drop ImmigrantGroup
            reshape wide Wage Supply, i(statefip year) j(abbr) s

            * Labor wage and quantity variables
            qui ds Wage_*
            foreach v in `r(varlist)' {
                
                * Put longer group name in variable labels
                loc origin = substr("`v'", 6, .)
                if substr("`v'", 6, .) == "CaAuNz" loc origin "Canada-Australia-New Zealand"
                if substr("`v'", 6, .) == "US" loc origin "United States"
                if substr("`v'", 6, .) == "WestEu" loc origin "Western Europe"
                if substr("`v'", 6, .) == "LA" loc origin "Latin America"
                if substr("`v'", 6, .) == "EastEu" loc origin "Russia and Eastern Europe"
                if substr("`v'", 6, .) == "AsiaOther" loc origin "Asia Other"
                
                la var `v' "Foreign-born wage (2017 Dollars), `origin'"
            }

            qui ds Supply_*
            foreach v in `r(varlist)' {

                loc origin = substr("`v'", 8, .)
                if substr("`v'", 8, .) == "CaAuNz" loc origin "Canada-Australia-New Zealand"
                if substr("`v'", 8, .) == "US" loc origin "United States"
                if substr("`v'", 8, .) == "WestEu" loc origin "Western Europe"
                if substr("`v'", 8, .) == "LA" loc origin "Latin America"
                if substr("`v'", 8, .) == "EastEu" loc origin "Russia and Eastern Europe"
                if substr("`v'", 8, .) == "AsiaOther" loc origin "Asia Other"

                la var `v' "Foreign-born quantity, `origin'"
            }
        }

        xframeappend appendthis, drop
    }

    destring year, replace
}
/*********************************
    CLEAN AND MERGE STATE GDP
**********************************/
if !`skipBea' {

    /*********************************
            CLEAN 1963 to 1997
    **********************************/
    frame create StateGdp63to97 
    frame StateGdp63to97 {

        import delimited "${GdpData}/SAGDP_SIC/SAGDP2S__ALL_AREAS_1963_1997.csv", clear

        * Only interested in industry totals for each state
        keep if description == "All industry total"
        drop if geoname == "United States"

        * Extract the numerical part of the fips code + drop trailing zeros
        replace geofips = substr(geofips, 3, 2)

        * Keep just the identifiers and the values - note units are millions of USD
        keep geofips geoname v*
        qui ds v*
        foreach v in `r(varlist)' { // Variable labels have the year (beacuse years are invalid variable names)

            loc yyyy: var lab `v'
            ren `v' GDP_`yyyy'
            destring GDP_`yyyy', replace
            la var GDP_`yyyy' ""

        }

        * Reshape long
        reshape long GDP_, i(geofips geoname) j(year)
        rename GDP_ GDP

        * Prepare for merge with the main frame
        ren geofips statefip
        ren geoname state

        * Drop 1997 since it is in the Naics-based file too
        drop if year == 1997

    }

    /*********************************
            CLEAN 1998 to 2023
    **********************************/
    frame create StateGdp97to23
    frame StateGdp97to23 {

        import delimited "${GdpData}/SAGDP/SAGDP2N__ALL_AREAS_1997_2023.csv", clear
        
        * Only interested in industry totals for each state
        keep if linecode == 3
        drop if geoname == "United States *"
    
        * Extract the numerical part of the fips code + drop trailing zeros
        replace geofips = substr(geofips, 3, 2)

        * Keep just the identifiers and the values - note units are millions of USD
        keep geofips geoname v*
        qui ds v*
        foreach v in `r(varlist)' { // Variable labels have the year (beacuse years are invalid variable names)

            loc yyyy: var lab `v'
            ren `v' GDP_`yyyy'
            destring GDP_`yyyy', replace
            la var GDP_`yyyy' ""

        }

        * Reshape long
        reshape long GDP_, i(geofips geoname) j(year)
        ren GDP_ GDP

        * Prepare for merge with the main frame
        ren geofips statefip
        ren geoname state

    }

    /*********************************
        Append Frames
    **********************************/    
    frame StateGdp63to97 {
        xframeappend StateGdp97to23, drop
        recast str2 statefip, force
        replace GDP = GDP * 1e+6        // Units were millions
        la var GDP "Real GDP, 2017 Dollars" // Still nominal at this point but will be deflated later
        tempfile StateGdp
        save `StateGdp', replace
    }
}

/*********************************
    CLEAN AND MERGE EL-SHAGI AND YAMARIK
**********************************/
if !`skipEY' {

    * El-Shagi and Yamarik cap stocks are real 2009 dollars. Let's reinflate
    frame create InvDeflator 
    frame InvDeflator {

        import delimited "${Data}/Price Deflators/InvestmentDeflator.csv", clear
        gen daten = date(observation_date,"YMD")
        format daten %td

        gen year = yofd(daten)
        ren a pk
        drop daten observation
        drop if mi(pk)
        keep if inrange(year, 1947, 2021)

        * Change base year to 2009 so that we can reinflate
        qui summ pk if year == 2009
        replace pk = pk / `r(mean)'
        keep year pk
        label var pk ""
    }

    frame create CapStocks
    frame CapStocks {
        
        use "${CapStock}/state_capital_yesdata21.dta", clear
    
        * Make statefip consistent with other data
        tostring fips, gen(statefip)
        drop fips
        replace statefip = "0" + statefip if strlen(statefip) == 1
        keep statefip year cap
        replace cap = cap * 1e+6 // Eyeballing the graph on the website, national stock should be about 21 trillion -> units here are in millions

        * Reinflate
        frlink m:1 year, frame(InvDeflator)
        frget pk, from(InvDeflator)
        drop InvDeflator
        frame drop InvDeflator
        gen CapStock = cap * pk
        drop pk cap

        tempfile capstock
        save `capstock', replace
    }
    
}

/*********************************
    CLEAN AND MERGE STATE GDP
**********************************/
if !`skipPrePeriod' {

    frame create PrePeriod
    frame PrePeriod {
        
        use "${Data}/PrePeriod.dta", clear
        tostring statefip, replace
        replace statefip = "0" + statefip if strlen(statefip) == 1
        
        qui ds statefip, not
        foreach v in `r(varlist)' {

            loc yyyy = substr("`v'", -4, .)
            loc origin = subinstr("`v'", "`yyyy'", "", 1)
            if "`origin'" == "CaAuNz" loc origin "Canada-Australia-New Zealand"
            if "`origin'" == "US" loc origin "United States"
            if "`origin'" == "WestEu" loc origin "Western Europe"
            if "`origin'" == "LA" loc origin "Latin America"
            if "`origin'" == "EastEu" loc origin "Russia and Eastern Europe"
            if "`origin'" == "AsiaOther" loc origin "Asia Other"

            la var `v' "Preperiod quantity (`yyyy'), `origin'"
        }

        tempfile preperiod
        save `preperiod', replace
    }
}

/*********************************
        MERGE IT ALL TOGETHER
**********************************/
if !`skipMerge' {
    
    * IMPUS + STATEGDP
    * The only year in main only was 2024 (51 observations)
    merge 1:1 statefip year using `StateGdp', nogen keep(3)
    order state statefip year
    sort statefip year

    * ADD IN EL-SHAGI YAMARIK
    merge 1:1 statefip year using `capstock', nogen keep(3)

    * ADD IN THE PRE PERIOD SHARES
    merge m:1 statefip using `preperiod', nogen

}

* Label some stuff
la var state ""
la var CapStock "Real Capital Stock, 2017 Dollars"
la var Supply_US "Domestic-born quantity"
la var Wage_US "Domestic-born wage (2017 Dollars)"

* Convert to common base year (2017)
frame create InvDeflator 
frame InvDeflator {

    import delimited "${Data}/Price Deflators/InvestmentDeflator.csv", clear
    gen daten = date(observation_date,"YMD")
    format daten %td

    gen year = yofd(daten)
    gen pk = a / 100
    la var pk ""
    drop daten observation a
    keep if inrange(year, 1994, 2021)

    tempfile pk
    save `pk', replace
}

frame create GdpDeflator
frame GdpDeflator {

    import delimited "${Data}/Price Deflators/GdpDeflator.csv", clear
    gen daten = date(observation_date,"YMD")
    format daten %td

    gen year = yofd(daten)
    gen p = a / 100
    la var p ""
    drop daten observation a
    keep if inrange(year, 1994, 2021)

    tempfile p
    save `p', replace

}

merge m:1 year using `pk', nogen
merge m:1 year using `p', nogen

qui ds Wage_*
foreach v in `r(varlist)' {
    
    replace `v' = `v' / p
}

replace GDP = GDP / p
replace CapStock = CapStock / pk
drop p pk

/************************************
    CALCULATE FOREIGN AND DOMESTIC SUPPLIES/WAGES
************************************/
* Quantities
egen Supply_Total    = rowtotal(Supply_*)
gen  Supply_Foreign  = Supply_Total - Supply_US
gen  Supply_Domestic = Supply_US

* Wages (quantitiy weighted)
qui ds Wage_*
foreach v in `r(varlist)' {

    loc origin = substr("`v'", 6, .)
    gen Wage_weighted_`origin' = Supply_`origin' * Wage_`origin'

}

egen Wage_Total    = rowtotal(Wage_weighted_*)
gen  Wage_Foreign  = (Wage_Total - Wage_weighted_US) / Supply_Foreign
gen  Wage_Domestic = Wage_US
gen  Wage = Wage_Total / Supply_Total

drop Wage_weighted_* Wage_Total

* Add some labels
la var Supply_Total    "Quantity of labor, total"
la var Supply_Foreign  "Quantity of labor, foreign-born"
la var Supply_Domestic "Quantity of labor, domestic-born"
la var Wage_Foreign    "Wage, foreign-born"
la var Wage_Domestic   "Wage, domestic"
la var Wage            "Wage, average"

/***********
    SAVING
***********/
label save using VarLabels.do, replace
sort statefip year
save "${Data}/StateAnalysisPreTfp.dta", replace
