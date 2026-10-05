clear all
do Globals

/*================================================================
Loops once over data/acs/Acs2001.dta-Acs2024.dta to build three files:
  1. AcsPiPanel.dta          -- state stock/migration-flow panel; feeds
     EstimateScaleBetaAcs.do's nu^D/nu^F and EstimateBilateralCosts.do's f^n_ll'
  2. IndividualCpAnalysis.dta -- individual-level file for the choice-
     probability probit (its state-year relative wage is the ratio of file
     3's composition-adjusted wages); feeds EstimateCp.do
  3. StateWageSupplyPanel.dta -- state x year x nativity labor supply,
     hourly wages (unresidualized and composition-adjusted/residualized),
     gross in-migration from MIGPLAC1 (Inflow_Domestic is interstate-only;
     Inflow_Foreign also includes new arrivals from abroad, being the sum of
     the per-nationality Inflow<code> columns below), and, by nationality
     (BPL): current-period arrivals (Inflow<code>) and stock (Stock<code>)
     -- the numerator and denominator for a shift-share instrument's
     growth-rate term -- plus the subset of arrivals coming from abroad only
     (NewArr<code>, and NewArr_Foreign as their state total), which excludes
     foreign-born movers between US states. Feeds MakeLpivPanel.do.

All branches share UHRSWORK>=35. The individual-level probit file (2) and
the state wage/supply panel (3) additionally require a full-time-full-year
(FTFY) restriction (weeks worked>=40); the migration-flow branches (1, and
the nationality-specific arrivals) do not, since a weeks-worked floor
mechanically penalizes people who just moved (a move-year employment gap is
common) and distorts the mover/stayer counts nu is identified from. Branch
(3) builds the only wage measure, an hours/weeks-adjusted hourly wage, CPI-U
deflated to real 2009 dollars; branch (2) takes its relative wage from it,
and branch (1) carries no wages.

The earliest flow-panel year is 2002, since a Year=2001 flow needs a
Year=2000 origin-stock denominator that this extract doesn't include.
================================================================*/

frame create FlowPanel
frame create StockWagePanel
frame create IndividualPanel
frame create StatePersonPanel
frame create NationalityFlowPanel
frame create NationalityStockPanel
frame create NewArrivalsPanel

forval yr = 2001/2024 {

    use "${Data}/acs/Acs`yr'.dta", clear

    /*---------------- SHARED BASE RESTRICTIONS (full-time, full-year) -------------*/
    drop if AGE < 16
    drop if UHRSWORK < 35
    drop if BPL >= 900
    gen Domestic = BPL < 100
    gen Foreign  = 1 - Domestic

    * Unify weeks worked across ACS's WKSWORK1/WKSWORK2 coding change (WKSWORK1,
    * continuous, fielded 2001-2007 and 2019+; WKSWORK2, categorical intervals,
    * 2008-2018). WKSWORK2 categories mapped to interval midpoints: 1=1-13wks,
    * 2=14-26, 3=27-39, 4=40-47, 5=48-49, 6=50-52.
    capture confirm variable WKSWORK1
    if !_rc gen Weeks = WKSWORK1
    else    gen Weeks = .

    capture confirm variable WKSWORK2
    if !_rc {
        replace Weeks = 7    if mi(Weeks) & WKSWORK2 == 1
        replace Weeks = 20   if mi(Weeks) & WKSWORK2 == 2
        replace Weeks = 33   if mi(Weeks) & WKSWORK2 == 3
        replace Weeks = 43.5 if mi(Weeks) & WKSWORK2 == 4
        replace Weeks = 48.5 if mi(Weeks) & WKSWORK2 == 5
        replace Weeks = 51   if mi(Weeks) & WKSWORK2 == 6
    }

    * NOTE: the weeks-worked>=40 floor (full-time-full-year, on top of the
    * UHRSWORK>=35 already imposed above) is applied LOCALLY within branches
    * (2) and (3) below, not here -- branch (1) and the nationality-arrivals
    * branch keep the lighter UHRSWORK>=35-only restriction (see header note).

    * STATEFIP is intentionally left numeric here -- IndividualCpAnalysis.dta
    * (branch 2 below) keeps it numeric, while AcsPiPanel.dta (branch 1) and the
    * state panel (branch 3) need it as a zero-padded string; each of those
    * branches converts it locally, inside its own preserve block, so branch 2's
    * output isn't affected.
    replace INCWAGE = . if INCWAGE == 999999
    egen incwage_censor = pctile(INCWAGE), p(99) by(STATEFIP)
    replace INCWAGE = incwage_censor if INCWAGE >= incwage_censor

    /*================================================================
      BRANCH 1: stock cross-section and migration-flow panel
    ================================================================*/
    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        collapse (rawsum) Stock = PERWT, by(STATEFIP YEAR Domestic)
        reshape wide Stock, i(STATEFIP YEAR) j(Domestic)
        ren Stock0 Stock_Foreign
        ren Stock1 Stock_Domestic
        ren STATEFIP State
        ren YEAR Year

        frame StockWagePanel: xframeappend default
    restore

    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        capture confirm variable MIGPLAC1
        if !_rc {
            keep PERWT STATEFIP YEAR Domestic MIGPLAC1

            drop if MIGPLAC1 >= 100
            gen Origin = STATEFIP if MIGPLAC1 == 0
            replace Origin = string(MIGPLAC1, "%02.0f") if mi(Origin)
            drop MIGPLAC1

            ren STATEFIP Destination
            collapse (sum) PERWT, by(YEAR Destination Origin Domestic)
            ren PERWT Flow
            reshape wide Flow, i(YEAR Destination Origin) j(Domestic)
            fillin Destination Origin
            drop _fillin
            ren YEAR Year
            ren Flow0 Flow_Foreign
            ren Flow1 Flow_Domestic
            qui summ Year
            replace Year = `r(mean)' if mi(Year)
            replace Flow_Foreign  = 0 if mi(Flow_Foreign)
            replace Flow_Domestic = 0 if mi(Flow_Domestic)

            frame FlowPanel: xframeappend default
        }
    restore

    /*---------------- FOREIGN-BORN ARRIVALS BY NATIONALITY (BPL) -----------------*/
    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        capture confirm variable MIGPLAC1
        if !_rc {
            keep if Foreign == 1

            * "Arrived" = anything other than same-state (no move, or moved
            * within the same state); origin can be another US state or
            * abroad (MIGPLAC1>=100) -- both count as a new arrival to this
            * state for this nationality group.
            gen Origin = STATEFIP if MIGPLAC1 == 0
            replace Origin = string(MIGPLAC1, "%02.0f") if mi(Origin)
            drop if Origin == STATEFIP

            keep PERWT STATEFIP YEAR BPL
            ren STATEFIP Destination
            ren YEAR Year

            collapse (sum) PERWT, by(Destination Year BPL)
            ren PERWT Inflow

            frame NationalityFlowPanel: xframeappend default
        }
    restore

    /*---------------- NEW FOREIGN-BORN ARRIVALS FROM ABROAD, BY NATIONALITY ------*/
    * Same as the arrivals branch above but restricted to MIGPLAC1>=100 (residence
    * one year ago was outside the US), so interstate moves of the foreign born
    * are excluded.
    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        capture confirm variable MIGPLAC1
        if !_rc {
            keep if Foreign == 1 & MIGPLAC1 >= 100

            keep PERWT STATEFIP YEAR BPL
            ren STATEFIP Destination
            ren YEAR Year

            collapse (sum) PERWT, by(Destination Year BPL)
            ren PERWT NewArr

            frame NewArrivalsPanel: xframeappend default
        }
    restore

    /*---------------- FOREIGN-BORN STOCK BY NATIONALITY (BPL) --------------------*/
    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        keep if Foreign == 1
        collapse (sum) PERWT, by(STATEFIP YEAR BPL)
        ren PERWT Stock
        ren YEAR Year

        frame NationalityStockPanel: xframeappend default
    restore

    /*================================================================
      BRANCH 2: individual rows for the choice-probability probit
    ================================================================*/
    preserve
        drop if mi(Weeks) | Weeks < 40

        keep STATEFIP YEAR OCC1990 Foreign PERWT
        ren YEAR Year

        frame IndividualPanel: xframeappend default
    restore

    /*================================================================
      BRANCH 3: state-level wage/supply panel -- hours/weeks-adjusted
      hourly wage, full-time-full-year (FTFY) sample
    ================================================================*/
    preserve
        drop if mi(Weeks) | Weeks < 40

        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        gen HourlyWage = INCWAGE / (UHRSWORK * Weeks)
        drop if mi(HourlyWage) | HourlyWage <= 0

        keep STATEFIP YEAR Domestic PERWT HourlyWage AGE EDUC SEX RACE BPL
        ren YEAR Year

        frame StatePersonPanel: xframeappend default
    restore

}

/*================================================================
  BRANCH 1 OUTPUT: AcsPiPanel.dta
================================================================*/
frame change StockWagePanel
capture frame drop default
tempfile StockWageFile
save `StockWageFile'

frame change FlowPanel
frame drop StockWagePanel

gen t = Year - 1

preserve
    use `StockWageFile', clear
    keep State Year Stock_Domestic Stock_Foreign
    ren State Origin
    ren Year t
    ren Stock_Domestic Supply_Domestic_Origin
    ren Stock_Foreign  Supply_Foreign_Origin
    tempfile OriginStock
    save `OriginStock'
restore

merge m:1 Origin t using `OriginStock', keep(1 3) nogen

preserve
    use `StockWageFile', clear
    keep State Year Stock_Domestic Stock_Foreign
    ren State Origin
    ren Stock_Domestic Supply_Domestic_Origin_lead
    ren Stock_Foreign  Supply_Foreign_Origin_lead
    tempfile OriginStockLead
    save `OriginStockLead'
restore

merge m:1 Origin Year using `OriginStockLead', keep(1 3) nogen

preserve
    use `StockWageFile', clear
    keep State Year Stock_Domestic Stock_Foreign
    ren State Destination
    ren Stock_Domestic Supply_Domestic_Dest
    ren Stock_Foreign  Supply_Foreign_Dest
    tempfile DestStock
    save `DestStock'
restore

merge m:1 Destination Year using `DestStock', keep(1 3) nogen

la var Origin                    "Origin state (l), FIPS"
la var Destination                "Destination state (l'), FIPS"
la var Year                       "Survey year Y = t+1 (destination/arrival year)"
la var t                          "Origin period, t = Y-1"
la var Flow_Domestic               "Weighted domestic-born flow, Origin -> Destination, realized in Year"
la var Flow_Foreign                "Weighted foreign-born flow, Origin -> Destination, realized in Year"
la var Supply_Domestic_Origin      "Origin domestic-born labor stock at t (ACS)"
la var Supply_Foreign_Origin       "Origin foreign-born labor stock at t (ACS)"
la var Supply_Domestic_Origin_lead "Origin domestic-born labor stock at t+1 (ACS)"
la var Supply_Foreign_Origin_lead  "Origin foreign-born labor stock at t+1 (ACS)"
la var Supply_Domestic_Dest        "Destination domestic-born labor stock at t+1 (ACS)"
la var Supply_Foreign_Dest         "Destination foreign-born labor stock at t+1 (ACS)"

order Origin Destination Year t Flow_Domestic Flow_Foreign ///
      Supply_Domestic_Origin Supply_Foreign_Origin ///
      Supply_Domestic_Origin_lead Supply_Foreign_Origin_lead ///
      Supply_Domestic_Dest Supply_Foreign_Dest
sort Origin Destination Year

save "${Data}/AcsPiPanel.dta", replace

* State-level gross in-migration (destination total across all origins,
* domestic and foreign born separately), for the state panel below -- a
* direct MIGPLAC1-based flow measure, rather than a first-difference of
* stocks. Excludes Origin==Destination (non-movers/in-state movers), which
* the bilateral panel above includes as a same-state "flow" for stay-rate
* calculations elsewhere but which isn't a real migration flow here.
preserve
    drop if Origin == Destination
    collapse (sum) Flow_Domestic, by(Destination Year)
    ren Destination STATEFIP
    ren Flow_Domestic Inflow_Domestic
    tempfile InflowFile
    save `InflowFile'
restore

* Foreign-born arrivals by nationality (BPL), wide by BPL code -- one column
* per country/region-of-birth code, e.g. Inflow200 for BPL==200 (Mexico).
* Inflow_Foreign (the state-level total) is built here, as the sum across all
* BPL codes, rather than from the interstate-only FlowPanel above -- it needs
* to include arrivals from abroad, which is exactly what this branch's own
* per-BPL columns already capture and FlowPanel's Flow_Foreign does not.
frame NationalityFlowPanel {
    capture frame drop default
    reshape wide Inflow, i(Destination Year) j(BPL)
    ds Inflow*
    foreach v of varlist `r(varlist)' {
        replace `v' = 0 if mi(`v')
    }
    egen double Inflow_Foreign = rowtotal(Inflow*)
    ren Destination STATEFIP
}
frame NationalityFlowPanel: tempfile NationalityInflowFile
frame NationalityFlowPanel: save `NationalityInflowFile'

* Arrivals from abroad only, wide by BPL code (NewArr200 for Mexico, etc.);
* NewArr_Foreign is the state total across BPL codes.
frame NewArrivalsPanel {
    reshape wide NewArr, i(Destination Year) j(BPL)
    ds NewArr*
    foreach v of varlist `r(varlist)' {
        replace `v' = 0 if mi(`v')
    }
    egen double NewArr_Foreign = rowtotal(NewArr*)
    ren Destination STATEFIP
}
frame NewArrivalsPanel: tempfile NewArrivalsFile
frame NewArrivalsPanel: save `NewArrivalsFile'

* Foreign-born stock by nationality (BPL), wide by BPL code -- one column
* per country/region-of-birth code, e.g. Stock200 for BPL==200 (Mexico).
* This is L^F_{m,l,t}: the growth-rate denominator for the eventual
* shift-share instrument (paired with the Inflow<BPL> columns above as the
* numerator), on the same sample restriction (UHRSWORK>=35, no FTFY floor)
* as the arrivals themselves, so numerator and denominator are consistent.
frame NationalityStockPanel {
    capture frame drop default
    reshape wide Stock, i(STATEFIP Year) j(BPL)
    ds Stock*
    foreach v of varlist `r(varlist)' {
        replace `v' = 0 if mi(`v')
    }
}
frame NationalityStockPanel: tempfile NationalityStockFile
frame NationalityStockPanel: save `NationalityStockFile'

/*================================================================
  BRANCH 3 OUTPUT: StateWageSupplyPanel.dta -- unresidualized and residualized
  hourly wages by state x year x nativity

  Residualized wage: one pooled regression,
    ln(HourlyWage) ~ AGE + AGE^2 + i.EDUC + i.SEX + i.RACE + year FE + birthplace FE,
  where the birthplace categories are domestic-born plus the world regions
  the shift-share instrument is built from (same BPL grouping as
  MakeLpivPanel.do). The composition-adjusted log wage is year FE +
  birthplace FE + individual residual, plus the sample-mean contribution of
  the demographic controls (one constant, common to everyone, so wages stay
  in dollars at average characteristics). Returns to demographics are common
  across nativities and states, so the domestic/foreign ratio of the
  collapsed series compares workers at the same characteristics. Each
  person's adjusted wage is exponentiated before averaging: the state x year
  x nativity series is a mean of wage levels, like the unresidualized one.
================================================================*/
frame change StatePersonPanel
capture frame drop default

gen lnHourlyWage = ln(HourlyWage)
gen AGE2 = AGE^2

* Birthplace category: 0 = domestic-born, 1-12 = world regions, 13 = any
* foreign BPL code outside those regions.
loc region_1  "200"
loc region_2  "210 250 260 299"
loc region_3  "300"
loc region_4  "150 700 710 100 105 110 115 160 199"
loc region_5  "400 401 402 404 405 410 411 412 413 414 420 421 425 426 429"
loc region_6  "430 433 434 436 438"
loc region_7  "450 451 452 453 454 455 456 457 460 461 462 465 499"
loc region_8  "500 501 502 509"
loc region_9  "511 512 513 514 515 516 517 518 519"
loc region_10 "520 521 522 524"
loc region_11 "531 532 534 535 536 537 540 541 542 543 544 548 549 599"
loc region_12 "600"

gen byte BplRegion = cond(Domestic == 1, 0, 13)
forval r = 1/12 {
    foreach code of local region_`r' {
        qui replace BplRegion = `r' if BPL == `code'
    }
}
la def BplRegion 0 "Domestic" 1 "Mexico" 2 "CentAmCarib" 3 "SouthAmerica" 4 "CanadaOceania" 5 "NWEurope" ///
                 6 "SEurope" 7 "EEurope" 8 "EAsia" 9 "SEAsia" 10 "SAsia" 11 "MidEast" 12 "Africa" 13 "Other"
la val BplRegion BplRegion

reghdfe lnHourlyWage AGE AGE2 i.EDUC i.SEX i.RACE [pw = PERWT], absorb(Year BplRegion) resid
predict double xb  if e(sample), xb
predict double d   if e(sample), d
predict double res if e(sample), residuals

qui summ xb [aw = PERWT], meanonly
gen double AdjLnWage = r(mean) + d + res
drop xb d res

gen double AdjWageLevel = exp(AdjLnWage)
tabstat HourlyWage AdjWageLevel [aw = PERWT], by(Domestic)
tabstat HourlyWage AdjWageLevel [aw = PERWT], by(BplRegion)
drop AdjWageLevel

/*---------------- UNRESIDUALIZED SERIES ------------------------------*/
preserve
    collapse (mean) HourlyWage [pw = PERWT], by(STATEFIP Year Domestic)
    reshape wide HourlyWage, i(STATEFIP Year) j(Domestic)
    ren HourlyWage0 Wage_Foreign_Unresid
    ren HourlyWage1 Wage_Domestic_Unresid
    tempfile UnresidFile
    save `UnresidFile'
restore

/*---------------- RESIDUALIZED SERIES ---------------------------------*/
preserve
    drop if mi(AdjLnWage)
    gen double AdjWage = exp(AdjLnWage)
    collapse (mean) AdjWage [pw = PERWT], by(STATEFIP Year Domestic)
    reshape wide AdjWage, i(STATEFIP Year) j(Domestic)
    ren AdjWage0 Wage_Foreign_Resid
    ren AdjWage1 Wage_Domestic_Resid
    tempfile ResidFile
    save `ResidFile'
restore

/*---------------- LABOR SUPPLY (FTFY headcount, for reference) -------*/
collapse (rawsum) Supply = PERWT, by(STATEFIP Year Domestic)
reshape wide Supply, i(STATEFIP Year) j(Domestic)
ren Supply0 Supply_Foreign
ren Supply1 Supply_Domestic

merge 1:1 STATEFIP Year using `UnresidFile', nogen
merge 1:1 STATEFIP Year using `ResidFile', nogen
merge 1:1 STATEFIP Year using `InflowFile', nogen
merge 1:1 STATEFIP Year using `NationalityInflowFile', nogen
merge 1:1 STATEFIP Year using `NationalityStockFile', nogen
merge 1:1 STATEFIP Year using `NewArrivalsFile', nogen

ds NewArr*
foreach v of varlist `r(varlist)' {
    replace `v' = 0 if mi(`v')
}

loc inflowdomfor "Inflow_Domestic Inflow_Foreign"
ds Inflow*
loc allinflowvars `r(varlist)'
loc bplinflowvars: list allinflowvars - inflowdomfor
foreach v of local bplinflowvars {
    replace `v' = 0 if mi(`v')
}

ds Stock*
foreach v of varlist `r(varlist)' {
    replace `v' = 0 if mi(`v')
}

/*---------------- CPI-U DEFLATION (2009 base, matches rest of pipeline) */
preserve
    import delimited "${Data}/Price Deflators/CpiDeflator.csv", clear varnames(1)
    ren year Year
    qui summ cpiaucns if Year == 2009
    gen p = cpiaucns / r(mean)
    keep Year p
    tempfile CpiFile
    save `CpiFile'
restore

merge m:1 Year using `CpiFile', keep(1 3) nogen

foreach v of varlist Wage_*_Unresid Wage_*_Resid {
    replace `v' = `v' / p
}
drop p

la var STATEFIP              "State FIPS code"
la var Year                  "Survey year"
la var Supply_Domestic       "Domestic-born FTFY labor supply (person-weighted headcount)"
la var Supply_Foreign        "Foreign-born FTFY labor supply (person-weighted headcount)"
la var Wage_Domestic_Unresid "Domestic-born hourly wage, PERWT-weighted mean, real 2009 USD (CPI-U deflated), FTFY sample"
la var Wage_Foreign_Unresid  "Foreign-born hourly wage, PERWT-weighted mean, real 2009 USD (CPI-U deflated), FTFY sample"
la var Wage_Domestic_Resid   "Domestic-born hourly wage, composition-adjusted (pooled regression on age/age2/educ/sex/race), real 2009 USD, FTFY sample"
la var Wage_Foreign_Resid    "Foreign-born hourly wage, composition-adjusted (pooled regression on age/age2/educ/sex/race), real 2009 USD, FTFY sample"
la var Inflow_Domestic       "Domestic-born gross in-migration from another US state, realized in Year (ACS MIGPLAC1, UHRSWORK>=35 sample, no FTFY floor, interstate only); missing for Year=2001 (needs a Year=2000 origin stock not in this extract)"
la var Inflow_Foreign        "Foreign-born gross in-migration, realized in Year (ACS MIGPLAC1, UHRSWORK>=35 sample, no FTFY floor; sum of the Inflow<BPL> columns below, so includes both interstate moves and new arrivals from abroad)"

ds Inflow*
loc allinflowvars `r(varlist)'
loc bplinflowvars: list allinflowvars - inflowdomfor
foreach v of local bplinflowvars {
    loc bplcode = substr("`v'", 7, .)
    la var `v' "Foreign-born arrivals with BPL==`bplcode', realized in Year (ACS MIGPLAC1, UHRSWORK>=35 sample, no FTFY floor; includes both interstate moves and new arrivals from abroad; Inflow_Foreign above is the sum of these)"
}

la var NewArr_Foreign "Foreign-born arrivals from abroad only (MIGPLAC1>=100), realized in Year (ACS, UHRSWORK>=35 sample, no FTFY floor); sum of the NewArr<BPL> columns, excludes foreign-born interstate movers"
ds NewArr*
loc allnewvars `r(varlist)'
loc newarrtot "NewArr_Foreign"
loc bplnewvars: list allnewvars - newarrtot
foreach v of local bplnewvars {
    loc bplcode = substr("`v'", 7, .)
    la var `v' "Foreign-born arrivals from abroad only (MIGPLAC1>=100) with BPL==`bplcode', realized in Year (ACS, UHRSWORK>=35 sample, no FTFY floor)"
}

ds Stock*
foreach v of varlist `r(varlist)' {
    loc bplcode = substr("`v'", 6, .)
    la var `v' "Foreign-born stock (L^F_m,l,t) with BPL==`bplcode' in Year (ACS, UHRSWORK>=35 sample, no FTFY floor -- same sample as the matching Inflow<BPL> column, for a consistent growth-rate numerator/denominator)"
}

order STATEFIP Year Supply_Domestic Supply_Foreign ///
      Wage_Domestic_Unresid Wage_Foreign_Unresid Wage_Domestic_Resid Wage_Foreign_Resid ///
      Inflow_Domestic Inflow_Foreign
sort STATEFIP Year

save "${Data}/StateWageSupplyPanel.dta", replace
di as text "Wrote ${Data}/StateWageSupplyPanel.dta"

/*================================================================
  BRANCH 2 OUTPUT: IndividualCpAnalysis.dta
================================================================*/
* State-year log relative wage for the probit, from the composition-adjusted
* wages of the state panel just saved (still in memory).
keep STATEFIP Year Wage_Domestic_Resid Wage_Foreign_Resid
gen w_tilde = ln(Wage_Domestic_Resid / Wage_Foreign_Resid)
la var w_tilde "ln(Wage_Domestic_Resid/Wage_Foreign_Resid), state-year composition-adjusted hourly wages (StateWageSupplyPanel.dta)"
destring STATEFIP, replace
keep STATEFIP Year w_tilde
tempfile StateWageFile
save `StateWageFile'

use "${PeriSparber}/Abilities_occ1990.dta", clear
gen manual   = (1/19)*(im_a22+im_a23+im_a24+im_a25+im_a26+im_a27+im_a28+im_a29 ///
                       +im_a30+im_a31+im_a32+im_a33+im_a34+im_a35+im_a36+im_a37 ///
                       +im_a38+im_a39+im_a40)
gen language = (1/4)*(im_a1+im_a2+im_a3+im_a4)
gen cm_ratio = language / manual
keep occ1990 manual language cm_ratio
ren occ1990 OCC1990
la var cm_ratio "Peri-Sparber basic communication/manual ratio (language/manual), from their own ICPSR replication code"
tempfile TaskAbility
save `TaskAbility'

frame change IndividualPanel
capture frame drop default

merge m:1 STATEFIP Year using `StateWageFile', keep(1 3) nogen
merge m:1 OCC1990 using `TaskAbility', keep(1 3)
di as text "--- OCC1990 merge with Peri-Sparber's occupation-ability index ---"
tab _merge
count if _merge == 1
drop if _merge == 1
drop _merge

preserve
    collapse (rawsum) TotalEmp = PERWT (mean) cm_ratio, by(OCC1990)
    gsort cm_ratio
    gen CumEmp = sum(TotalEmp)
    qui summ CumEmp, meanonly
    loc grand = r(max)
    gen tau = (CumEmp - 0.5*TotalEmp) / `grand'
    gen phi_inv_tau = invnormal(tau)
    la var tau         "Employment-weighted percentile rank of Peri-Sparber's c/m ratio across occ1990 categories, pooled over all sample years"
    la var phi_inv_tau "Phi^-1(tau) -- the model's task-index regressor"
    keep OCC1990 tau phi_inv_tau cm_ratio TotalEmp
    tempfile TauLookup
    save `TauLookup'
restore

merge m:1 OCC1990 using `TauLookup', keep(1 3) nogen

drop manual language

la var STATEFIP    "State FIPS code"
la var Year        "Survey year"
la var OCC1990     "IPUMS occupation, 1990 basis"
la var Foreign     "1 = foreign-born (BPL>=100), 0 = domestic-born -- the probit outcome, Pr(Foreign)=phi(tau)"
la var PERWT       "ACS person weight"

order STATEFIP Year OCC1990 Foreign PERWT w_tilde tau phi_inv_tau cm_ratio TotalEmp
sort STATEFIP Year OCC1990

count
tab Foreign
summ w_tilde tau phi_inv_tau cm_ratio

save "${Data}/IndividualCpAnalysis.dta", replace
di as text "Wrote ${Data}/IndividualCpAnalysis.dta"
