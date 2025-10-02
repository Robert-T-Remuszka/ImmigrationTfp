clear all
do Globals

* Looping trhough each decade and implementing a set of uniform cleaning code
foreach yyyy of numlist 1920(10)1990 {
    
    frame create PrePeriod

    frame PrePeriod {

        import delimited "${PrePeriod}/PrePeriod`yyyy'.csv", clear
        tostring bpld, replace

        * Make bpld a string of uniform length
        replace bpld = "0" + bpld if strlen(bpld) == 4
        replace bpld = "00" + bpld if strlen(bpld) == 3
        replace bpld = "000" + bpld if strlen(bpld) == 2
        replace bpld = "0000" + bpld if strlen(bpld) == 1

        * Using the first number to identify those we don't have an origin for --- drop these
        drop if inlist(substr(bpld,1,1), "8", "9")

        * Drop those who don't have an age or are younger than 16
        drop if age == 999
        drop if age < 16

        * Assigning origins based on Peri (2012)
        gen ImmigrantGroup = "United States" if substr(bpld,1,1) == "0"
        replace ImmigrantGroup = "Mexico" if bpld == "20000"
        replace ImmigrantGroup = "Latin America" if bpld == "11000" | (inlist(substr(bpld,1,2), "21", "25", "26", "30") & bpld != "26030")
        replace ImmigrantGroup = "Western Europe" if inlist(bpld,"45300","45000") | inlist(substr(bpld,1,2), "40", "41", "42", "43")
        replace ImmigrantGroup = "Russia and Eastern Europe" if inlist(substr(bpld,1,2), "45", "46") & !inlist(bpld,"45300", "45000")
        replace ImmigrantGroup = "Canada-Australia-New Zealand" if inlist(bpld, "15000", "70020", "70010") | substr(bpld, 1, 2) == "15"
        replace ImmigrantGroup = "China" if bpld == "50000"
        replace ImmigrantGroup = "India" if bpld == "52100"
        replace ImmigrantGroup = "Rest of Asia" if inlist(substr(bpld,1,2), "50", "51", "52","53", "54", "55") & !inlist(bpld,"52100","50000")
        replace ImmigrantGroup = "Africa" if substr(bpld,1,2) == "60"

        * Drop unassigned origins
        drop if mi(ImmigrantGroup)

    }
    
    xframeappend PrePeriod, drop

}

* Collapse
gen counter = 1
collapse (sum) counter [iw = perwt], by(ImmigrantGroup year statefip)

* Prepare for reshape
gen abbr = "_CaAuNz" if ImmigrantGroup == "Canada-Australia-New Zealand"
replace abbr = "_US" if ImmigrantGroup == "United States"
replace abbr = "_WestEu" if ImmigrantGroup == "Western Europe"
replace abbr = "_LA" if ImmigrantGroup == "Latin America"
replace abbr = "_EastEu" if ImmigrantGroup == "Russia and Eastern Europe"
replace abbr = "_AsiaOther" if ImmigrantGroup == "Rest of Asia"
replace abbr = "_" + ImmigrantGroup if mi(abbr)

* Reshape
tostring year, replace
drop ImmigrantGroup
reshape wide counter, i(statefip year) j(abbr) s
reshape wide counter_*, i(statefip) j(year) s

* Drop counter prefix in varnames
qui ds counter_*
foreach var in `r(varlist)' {
    loc newname = subinstr("`var'", "counter_", "", 1)
    ren `var' `newname'

    replace `newname' = 0 if mi(`newname')
}

save "${Data}/PrePeriod.dta", replace