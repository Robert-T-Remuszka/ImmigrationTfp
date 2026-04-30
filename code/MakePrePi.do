clear all
do Globals

loc MakeUsFlows 1
loc MakeUs2RowFlows 1

* Use 1995-1996 to initialize Pi's because UN migrant stock snapshot is for 1995
loc tailyear = 1996
loc headyear = `tailyear' - 1


/******************* US DESTINATIONS *************************/    
****** Go through each ACS/CPS download and construct flows into US destinations
frame create US_Destination

* CPS
loc files: dir "${Data}/cps/" files "*.dta"
foreach f in `files' {

    if inlist(substr("`f'", 4, 4), "`tailyear'") {
        
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
        ren STATE Destination
        collapse (sum) ASECWT, by(YEAR Destination Origin Domestic)
        ren ASECWT Flow
        reshape wide Flow, i(YEAR Destination Origin) j(Domestic)
        fillin Destination Origin
        drop _fillin
        ren YEAR Year
        ren Flow0 Foreign
        ren Flow1 Domestic
        qui summ Year
        replace Year = `r(mean)' if mi(Year)
        replace Foreign = 0 if mi(Foreign)
        replace Domestic = 0 if mi(Domestic)

        frame US_Destination {
            xframeappend default
            tempfile UsFlows
            save `UsFlows', replace
        }

    }
}

/******************* Merge in Stocks *************************/
use statefip year Supply_Foreign Supply_Domestic using "${Data}/StateAnalysis.dta", clear    
keep if year == `tailyear' - 1
drop year
ren statefip Origin

merge 1:m Origin using "`UsFlows'", keep(2 3) nogen
order Origin Destination Foreign Domestic Supply_Foreign Supply_Domestic
sort Origin Destination
drop Year

/********************* Fill In ROW Stocks *************************/

* Calculate the total US labor force
gen Total_Labor_Orig = Supply_For + Supply_Dom
egen Total_Labor_US = total(Total_Labor_Orig), by(Destination)
loc Total_Labor_US = Total_Labor_US[1]
drop Tot*

* Get the UN Population estimate for the headyear
frame create UN_Pop
frame UN_Pop {

    import delimited "${UN}/WPP2024_Demographic_Indicators_Medium.csv", clear varn(1)
    keep if loctypename == "World" & time == `headyear'
    keep tpopulation1july
    gen GlobalPop = tpop * 1000 // convert to levels
    qui summ GlobalPop
    loc GlobalPop`headyear' = `r(mean)'
}
frame drop UN_Pop

* Calculate the US-born migrant stock outside of US
frame create UN_Stocks
frame UN_Stocks {

    import excel "${UN}/undesa_pd_2024_ims_stock_by_sex_destination_and_origin.xlsx", clear sheet("Table 1") cellrange(A11:AE28041) first
    tostring Locationcodeofdestination, gen(Destination_code)
    tostring Locationcodeoforigin, gen(Origin_code)
    ren Regiondevelopmentgroupcount Destination
    ren F Origin
    ren I Stock
    keep Origin Destination Origin_code Destination_code Stock
    replace Destination_code = "000" + Destination_code if strlen(Destination_code) == 1
    replace Destination_code = "00" + Destination_code if strlen(Destination_code) == 2
    replace Destination_code = "0"  + Destination_code if strlen(Destination_code) == 3
    replace Origin_code      = "000" + Origin_code if strlen(Origin_code)  == 1
    replace Origin_code      = "00" + Origin_code if strlen(Origin_code)  == 2
    replace Origin_code      = "0" + Origin_code if strlen(Origin_code)  == 3
    order Destination Destination_code Origin Origin_code

    keep if Origin == "United States of America*"
    keep if substr(Destination_code, 1, 1) == "0" & substr(Destination_code, 2, 1) != "9"

    collapse (sum) Tot_Domestic_ROW = Stock
    qui summ Tot
    loc Domestic_ROW = `r(mean)'
    
}
frame drop UN_Stocks

replace Supply_Domestic = `Domestic_ROW' if mi(Supply_Domestic)
replace Supply_Foreign  = `GlobalPop`headyear'' - `Total_Labor_US' - `Domestic_ROW' if mi(Supply_Foreign)

* Add ROW destination
preserve
    keep Origin
    duplicates drop
    gen Destination = "ROW"
    tempfile newobs
    save `newobs'
restore
append using `newobs'
sort Origin Destination

********************** Calculate the Pi's *************************
gen pi_Foreign  = Foreign  / Supply_Foreign 
gen pi_Domestic = Domestic / Supply_Domestic

* Floor for natural zeros (no observed flow to a US destination)
replace pi_Foreign  = 1e-4 if pi_Foreign  == 0 
replace pi_Domestic = 1e-4 if pi_Domestic == 0

* If observed outflow rates sum to > 0.99 for an origin, rescale proportionally
* so the sum equals 0.99, leaving at least 0.01 probability for ROW and probabilities sum to one
egen pi_tot_Foreign  = total(pi_Foreign),  by(Origin)
egen pi_tot_Domestic = total(pi_Domestic), by(Origin)
replace pi_Foreign  = pi_Foreign  * min(0.99 / pi_tot_Foreign,  1) if !mi(pi_Foreign)
replace pi_Domestic = pi_Domestic * min(0.99 / pi_tot_Domestic, 1) if !mi(pi_Domestic)
drop pi_tot_Foreign pi_tot_Domestic
egen pi_tot_Foreign  = total(pi_Foreign),  by(Origin)
egen pi_tot_Domestic = total(pi_Domestic), by(Origin)
replace pi_Foreign  = 1 - pi_tot_Foreign  if mi(pi_Foreign)
replace pi_Domestic = 1 - pi_tot_Domestic if mi(pi_Domestic)
drop pi_tot*

* Clean up the names and labels
ren Foreign DestinationIn_Foreign
ren Domestic DestinationIn_Domestic

la var DestinationIn_F "Foreign flow into destination from origin"
la var DestinationIn_D "Domestic flow into destination from origin"
la var pi_F "Foreigh choice prb of moving to destination from origin"
la var pi_D "Domestic choice prb of moving to destination from origin"

* Fill in the supplies for origign-destination (XX,ROW) pairs - recall that I added ROW on lines 166 to 175
bys Origin (Destination): replace Supply_Foreign = Supply_Foreign[_n-1] if Destination == "ROW"
bys Origin (Destination): replace Supply_Domestic = Supply_Domestic[_n-1] if Destination == "ROW"

* Name the supplies indicating the years
ren Supply_Foreign Supply_Foreign_`headyear'
ren Supply_Domestic Supply_Domestic_`headyear'
la var Supply_Foreign_ "Foreign Stock in Origin, `headyear'"
la var Supply_Domestic_ "Domestic Stock in Origin, `headyear'"

* Fill in the model-consistent flows to ROW
replace DestinationIn_F = pi_F * Supply_F if mi(DestinationIn_F)
replace DestinationIn_D = pi_D * Supply_D if mi(DestinationIn_D)

* Calculate total flows out of an origin
egen TotalOut_Domestic = total(DestinationIn_D), by(Origin)
egen TotalOut_Foreign  = total(DestinationIn_F), by(Origin)

* Calculate total flows into an origin
preserve
    
    collapse (sum) DestinationIn_D DestinationIn_F, by(Destination)
    
    ren DestinationIn_D TotalIn_Domestic
    ren DestinationIn_F TotalIn_Foreigin
    ren Destination Origin
    
    tempfile OriginInflows
    save `OriginInflows', replace

restore

merge m:1 Origin using "`OriginInflows'", nogen
gen Supply_Foreign_`tailyear'  = Supply_Foreign_`headyear'  + TotalIn_F - TotalOut_F
gen Supply_Domestic_`tailyear' = Supply_Domestic_`headyear' + TotalIn_D - TotalOut_D
la var Supply_Foreign_`tailyear' "Foreign Stock in Origin, `tailyear'"
la var Supply_Domestic_`tailyear' "Domestic Stock in Origin, `tailyear'"
drop Total* DestinationIn*

ren pi_F pi_F_
ren pi_D pi_D_
ren Supply_Foreign_`headyear' Foreign_`headyear'
ren Supply_Domestic_`headyear' Domestic_`headyear'
ren Supply_Foreign_`tailyear' Foreign_`tailyear'
ren Supply_Domestic_`tailyear' Domestic_`tailyear'
sort Origin Destination
reshape wide pi_F_ pi_D_, i(Origin Foreign_* Domestic_*) j(Destination) s

* Bring in State names for readablity
frame create FipNameCross
frame FipNameCross {
    use state statefip using "${Data}/StateAnalysis.dta", clear
    duplicates drop
    insobs 1
    replace state = "ROW" if mi(state)
    replace statefip = "ROW" if mi(statefip)
    ren state State
    ren statefip Origin
}

frlink 1:1 Origin, frame(FipNameCross)
frget State, from(FipNameCross)
drop FipNameCross
frame drop FipNameCross
order Origin State

save "${Data}/PiMat.dta", replace
