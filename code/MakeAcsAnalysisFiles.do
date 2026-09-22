clear all
do Globals

/*================================================================
Loops once over data/acs/Acs2001.dta-Acs2024.dta to build three files:
  1. AcsPiPanel.dta          -- state stock/wage/migration-flow panel; feeds
     EstimateScaleBetaAcs.do's nu^D/nu^F and EstimateBilateralCosts.do's f^n_ll'
  2. IndividualCpAnalysis.dta -- individual-level file for the choice-
     probability probit; feeds EstimateCp.do
  3. StateWageSupplyPanel.dta -- state x year x nativity labor supply,
     hourly wages (unresidualized and composition-adjusted/residualized),
     gross in-migration (from MIGPLAC1, aggregated across all US origins),
     and, by nationality (BPL): current-period arrivals (Inflow<code>) and
     stock (Stock<code>) -- the numerator and denominator for a shift-share
     instrument's growth-rate term. Not currently read by any other step in
     this pipeline.

All branches share UHRSWORK>=35. The individual-level probit file (2) and
the state wage/supply panel (3) additionally require a full-time-full-year
(FTFY) restriction (weeks worked>=40); the migration-flow branches (1, and
the nationality-specific arrivals) do not, since a weeks-worked floor
mechanically penalizes people who just moved (a move-year employment gap is
common) and distorts the mover/stayer counts nu is identified from. Branches
(1) and (2) use raw INCWAGE as their wage measure; branch (3) uses an
hours/weeks-adjusted hourly wage, CPI-U deflated to real 2009 dollars.

The earliest flow-panel year is 2002, since a Year=2001 flow needs a
Year=2000 origin-stock denominator that this extract doesn't include.
================================================================*/

frame create FlowPanel
frame create StockWagePanel
frame create IndividualStateWagePanel
frame create IndividualPanel
frame create StatePersonPanel
frame create NationalityFlowPanel
frame create NationalityStockPanel

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
      BRANCH 1: stock/wage cross-section and migration-flow panel
    ================================================================*/
    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        collapse (mean) Wage = INCWAGE (rawsum) Stock = PERWT [pw = PERWT], ///
            by(STATEFIP YEAR Domestic)
        reshape wide Wage Stock, i(STATEFIP YEAR) j(Domestic)
        ren Wage0  Wage_Foreign
        ren Wage1  Wage_Domestic
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
      BRANCH 2: state-year wage cross-section and individual rows for the
      choice-probability probit
    ================================================================*/
    preserve
        drop if mi(Weeks) | Weeks < 40

        collapse (mean) Wage = INCWAGE [pw = PERWT], by(STATEFIP YEAR Domestic)
        reshape wide Wage, i(STATEFIP YEAR) j(Domestic)
        ren Wage0 Wage_Foreign
        ren Wage1 Wage_Domestic
        ren YEAR Year

        frame IndividualStateWagePanel: xframeappend default
    restore

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

        keep STATEFIP YEAR Domestic PERWT HourlyWage AGE EDUC SEX RACE
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
    keep State Year Wage_Domestic Wage_Foreign
    ren State Origin
    ren Wage_Domestic Wage_Domestic_Origin
    ren Wage_Foreign  Wage_Foreign_Origin
    tempfile OriginWage
    save `OriginWage'
restore

merge m:1 Origin Year using `OriginWage', keep(1 3) nogen

preserve
    use `StockWageFile', clear
    keep State Year Wage_Domestic Wage_Foreign Stock_Domestic Stock_Foreign
    ren State Destination
    ren Wage_Domestic Wage_Domestic_Dest
    ren Wage_Foreign  Wage_Foreign_Dest
    ren Stock_Domestic Supply_Domestic_Dest
    ren Stock_Foreign  Supply_Foreign_Dest
    tempfile DestWage
    save `DestWage'
restore

merge m:1 Destination Year using `DestWage', keep(1 3) nogen

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
la var Wage_Domestic_Origin        "Origin domestic-born mean wage at t+1 (ACS, nominal)"
la var Wage_Foreign_Origin         "Origin foreign-born mean wage at t+1 (ACS, nominal)"
la var Wage_Domestic_Dest          "Destination domestic-born mean wage at t+1 (ACS, nominal)"
la var Wage_Foreign_Dest           "Destination foreign-born mean wage at t+1 (ACS, nominal)"

order Origin Destination Year t Flow_Domestic Flow_Foreign ///
      Supply_Domestic_Origin Supply_Foreign_Origin ///
      Supply_Domestic_Origin_lead Supply_Foreign_Origin_lead ///
      Supply_Domestic_Dest Supply_Foreign_Dest ///
      Wage_Domestic_Origin Wage_Foreign_Origin Wage_Domestic_Dest Wage_Foreign_Dest
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
    collapse (sum) Flow_Domestic Flow_Foreign, by(Destination Year)
    ren Destination STATEFIP
    ren Flow_Domestic Inflow_Domestic
    ren Flow_Foreign  Inflow_Foreign
    tempfile InflowFile
    save `InflowFile'
restore

* Foreign-born arrivals by nationality (BPL), wide by BPL code -- one column
* per country/region-of-birth code, e.g. Inflow200 for BPL==200 (Mexico).
frame NationalityFlowPanel {
    capture frame drop default
    reshape wide Inflow, i(Destination Year) j(BPL)
    ds Inflow*
    foreach v of varlist `r(varlist)' {
        replace `v' = 0 if mi(`v')
    }
    ren Destination STATEFIP
}
frame NationalityFlowPanel: tempfile NationalityInflowFile
frame NationalityFlowPanel: save `NationalityInflowFile'

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
  BRANCH 2 OUTPUT: IndividualCpAnalysis.dta
================================================================*/
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

frame IndividualStateWagePanel {
    capture frame drop default
    gen w_tilde = ln(Wage_Domestic / Wage_Foreign)
    la var w_tilde "ln(Wage_Domestic/Wage_Foreign), state-year, ACS, PERWT-weighted mean, top-1%-by-state winsorized (NOT composition-adjusted yet)"
    keep STATEFIP Year w_tilde
}
frame IndividualStateWagePanel: tempfile StateWageFile
frame IndividualStateWagePanel: save `StateWageFile'

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

/*================================================================
  BRANCH 3 OUTPUT: StateWageSupplyPanel.dta -- unresidualized and residualized
  hourly wages by state x year x nativity

  Residualized wage: ln(HourlyWage) ~ AGE + AGE^2 + i.EDUC + i.SEX + i.RACE,
  year fixed effects only (not state x year -- returns to demographics are
  assumed common across states, so beta is identified off pooled/between-
  state variation and the composition-adjusted state-year level is
  recovered by collapsing the residual, not by saturating with state x year
  dummies), estimated SEPARATELY by nativity (Domestic, Foreign) since
  returns to age/education plausibly differ by nativity and pooling them
  would let compositional differences contaminate the residual.
================================================================*/
frame change StatePersonPanel
capture frame drop default

gen lnHourlyWage = ln(HourlyWage)
gen AGE2 = AGE^2

gen AdjLnWage = .

* By linearity, the PERWT-weighted mean of X_i'beta-hat across individuals
* equals Xbar'beta-hat (the fitted value AT reference/mean demographics) --
* no need for a synthetic reference observation or margins call.
qui levelsof Domestic, local(nativities)
foreach n of local nativities {
    qui reghdfe lnHourlyWage AGE AGE2 i.EDUC i.SEX i.RACE if Domestic == `n' [pw = PERWT], absorb(Year) resid
    predict double xb_`n'  if e(sample), xb
    predict double d_`n'   if e(sample), d
    predict double res_`n' if e(sample), residuals

    qui summ xb_`n' [aw = PERWT] if Domestic == `n', meanonly
    loc xbbar = r(mean)

    * Adjusted log wage = Xbar'beta-hat (reference demographics, nativity-
    * specific, fixed across the whole panel) + year FE + individual
    * residual -- the year FE (not state x year) is the only nationally
    * common component, so all remaining state-to-state variation in the
    * collapsed series comes through the residual, exactly as intended.
    replace AdjLnWage = `xbbar' + d_`n' + res_`n' if Domestic == `n' & e(sample)
    drop xb_`n' d_`n' res_`n'
}

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
    collapse (mean) AdjLnWage [pw = PERWT], by(STATEFIP Year Domestic)
    gen AdjWage = exp(AdjLnWage)
    drop AdjLnWage
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
la var Wage_Domestic_Resid   "Domestic-born hourly wage, composition-adjusted (age/age2/educ/sex/race, year FE), real 2009 USD, FTFY sample"
la var Wage_Foreign_Resid    "Foreign-born hourly wage, composition-adjusted (age/age2/educ/sex/race, year FE), real 2009 USD, FTFY sample"
la var Inflow_Domestic       "Domestic-born gross in-migration from another US state, realized in Year (ACS MIGPLAC1, UHRSWORK>=35 sample, no FTFY floor, interstate only); missing for Year=2001 (needs a Year=2000 origin stock not in this extract)"
la var Inflow_Foreign        "Foreign-born gross in-migration from another US state, realized in Year (ACS MIGPLAC1, UHRSWORK>=35 sample, no FTFY floor, interstate only -- excludes new arrivals from abroad, unlike the Inflow<BPL> columns below); missing for Year=2001"

ds Inflow*
loc allinflowvars `r(varlist)'
loc bplinflowvars: list allinflowvars - inflowdomfor
foreach v of local bplinflowvars {
    loc bplcode = substr("`v'", 7, .)
    la var `v' "Foreign-born arrivals with BPL==`bplcode', realized in Year (ACS MIGPLAC1, UHRSWORK>=35 sample, no FTFY floor; includes both interstate moves and new arrivals from abroad, unlike Inflow_Foreign above)"
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
