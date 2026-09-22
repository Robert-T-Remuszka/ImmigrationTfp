clear all
do Globals

/*================================================================
Consolidates the three ACS-derived analysis files that previously each ran
their own independent loop over data/acs/Acs*.dta:
  1. MakeAcsPi.do              -> AcsPiPanel.dta        (state stock/wage/
     migration-flow panel; feeds EstimateScaleBetaAcs.do's nu^D/nu^F and
     EstimateBilateralCosts.do's f^n_ll')
  2. MakeIndividualAnalysis.do -> IndividualCpAnalysis.dta (individual-level
     file for the choice-probability probit; feeds EstimateCp.do)
  3. NEW state-level wage/supply panel -> StateWageSupplyPanel.dta (feeds
     AggSupply_Estimate.jl's rho estimation; replaces
     MakeStateAnalysisPreTfp.do now that state-level capital/output are no
     longer needed)

Branches (1) and (2) are reproduced BIT-FOR-BIT on their existing sample
restriction (UHRSWORK>=35, no weeks-worked floor) and wage definition (raw
INCWAGE, no hours/weeks adjustment) -- nu, bilateral costs, and the probit
estimates should not move as a side effect of this consolidation. Branch (3)
uses a full-time-full-year (FTFY) restriction (UHRSWORK>=35 AND weeks
worked>=40) and an hours/weeks-adjusted hourly wage, CPI-U deflated. This
split is deliberate and temporary: branches (1)/(2) are to be raised to the
same FTFY/hourly-wage standard in a later pass, at which point nu, bilateral
costs, and the probit will need re-estimating.

Sample range is 2001-2024 (ACS's "residence 1 year ago" question, MIGPLAC1,
isn't fielded before 2001; CPS is no longer needed to fill in 2023-2024 now
that ACS itself covers them). One consequence of dropping the pre-2001
extract entirely: the earliest flow-panel year shifts from 2001 to 2002,
since a Year=2001 flow needs a Year=2000 origin-stock denominator that's no
longer pulled.
================================================================*/

frame create FlowPanel
frame create StockWagePanel
frame create IndividualStateWagePanel
frame create IndividualPanel
frame create StatePersonPanel

forval yr = 2001/2024 {

    use "${Data}/acs/Acs`yr'.dta", clear

    /*---------------- SHARED BASE RESTRICTIONS (identical to the old files) -------*/
    drop if AGE < 16
    drop if UHRSWORK < 35
    drop if BPL >= 900
    gen Domestic = BPL < 100
    gen Foreign  = 1 - Domestic

    * NOTE: STATEFIP is intentionally left numeric here -- MakeIndividualAnalysis.do
    * (branch 2 below) always kept it numeric, while MakeAcsPi.do (branch 1) and the
    * new state panel (branch 3) need it as a zero-padded string; each of those
    * branches converts it locally, inside its own preserve block, so branch 2's
    * output isn't affected.
    replace INCWAGE = . if INCWAGE == 999999
    egen incwage_censor = pctile(INCWAGE), p(99) by(STATEFIP)
    replace INCWAGE = incwage_censor if INCWAGE >= incwage_censor

    /*================================================================
      BRANCH 1: MakeAcsPi.do's stock/wage cross-section and flow panel,
      unchanged
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

    /*================================================================
      BRANCH 2: MakeIndividualAnalysis.do's state-year wage cross-section
      and individual rows, unchanged
    ================================================================*/
    preserve
        collapse (mean) Wage = INCWAGE [pw = PERWT], by(STATEFIP YEAR Domestic)
        reshape wide Wage, i(STATEFIP YEAR) j(Domestic)
        ren Wage0 Wage_Foreign
        ren Wage1 Wage_Domestic
        ren YEAR Year

        frame IndividualStateWagePanel: xframeappend default
    restore

    preserve
        keep STATEFIP YEAR OCC1990 Foreign PERWT
        ren YEAR Year

        frame IndividualPanel: xframeappend default
    restore

    /*================================================================
      BRANCH 3: NEW state-level wage/supply panel -- FTFY restriction and
      hours/weeks-adjusted hourly wage
    ================================================================*/
    preserve
        tostring STATEFIP, replace
        replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

        * Unify weeks worked across ACS's WKSWORK1/WKSWORK2 coding change
        * (WKSWORK1, continuous, fielded 2001-2007 and 2019+; WKSWORK2,
        * categorical intervals, 2008-2018). WKSWORK2 categories mapped to
        * interval midpoints: 1=1-13wks, 2=14-26, 3=27-39, 4=40-47, 5=48-49,
        * 6=50-52.
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

        * Full-time, full-year (FTFY): UHRSWORK>=35 already imposed above;
        * add the weeks-worked floor.
        drop if mi(Weeks) | Weeks < 40

        gen HourlyWage = INCWAGE / (UHRSWORK * Weeks)
        drop if mi(HourlyWage) | HourlyWage <= 0

        keep STATEFIP YEAR Domestic PERWT HourlyWage AGE EDUC SEX RACE
        ren YEAR Year

        frame StatePersonPanel: xframeappend default
    restore

}

/*================================================================
  BRANCH 1 OUTPUT: AcsPiPanel.dta (identical construction to MakeAcsPi.do)
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

/*================================================================
  BRANCH 2 OUTPUT: IndividualCpAnalysis.dta (identical construction to
  MakeIndividualAnalysis.do)
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

order STATEFIP Year Supply_Domestic Supply_Foreign ///
      Wage_Domestic_Unresid Wage_Foreign_Unresid Wage_Domestic_Resid Wage_Foreign_Resid
sort STATEFIP Year

save "${Data}/StateWageSupplyPanel.dta", replace
di as text "Wrote ${Data}/StateWageSupplyPanel.dta"
