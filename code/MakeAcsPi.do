clear all
do Globals

/*================================================================
ACS analog of MakePi.do. Builds a LONG (Origin, Destination, Year)
state-to-state migration FLOW panel from ACS 2000-2022, for the same
eq.-(27) analog (Nu_Derivation.md) that MakePi.do builds from CPS-ASEC.

WHY THIS FILE EXISTS (2026-08-20): the CPS-ASEC version of this panel
gave correctly-specified but statistically unusable estimates of
(nu^D, nu^F) -- confirmed by direct inspection (sign conventions, merge
lead/lag directions, pairid directionality all correct; Hansen J passes)
-- because CPS-ASEC's ~75k households/year gives only 1-5 respondents
per (Origin,Destination,Year) migration corridor, so ~80% of Domestic
and ~93% of Foreign cells are exact zeros and the wage-ratio regressor's
within-pair signal (pair-demeaned SD ~0.054) is swamped by clustered
noise. CDP themselves flag exactly this problem (p.763 fn.36: March CPS
is too thin to build interregional migration flows) and sidestep it by
using ACM's pre-aggregated flow matrices instead. ACS gives a structural
fix rather than a workaround: ~2-3M person records/year (30-40x CPS-ASEC),
so corridor cells are far less frequently zero and far less noisy.

Variable mapping from MakePi.do (CPS) to here (ACS), confirmed against
Acs2019.dta's actual coding before writing this file:
  ASECWT      -> PERWT       (person weight; no separate migration weight in ACS)
  MIGSTA1     -> MIGPLAC1    (prior-year state; SAME coding convention: 1-56 are
                               FIPS state codes, but ACS's "did not move" sentinel
                               is 0, not CPS's "99", and "abroad" is >=100, not
                               CPS's single code "91")
  UHRSWORKT   -> UHRSWORK    (ACS ranges 0-99 with 0 itself as the not-in-universe
                               sentinel -- no separate 999/998-style code the way
                               CPS's UHRSWORKT has, so `< 35` alone already excludes
                               NIU without CPS's extra `>= 997` clause)
  INCWAGE     -> INCWAGE     (same variable name; topcode/missing sentinel is
                               999999 only -- ACS has no second code analogous to
                               CPS's 999998)
  CITIZEN==9  -> (dropped)   (ACS's CITIZEN has no NIU/unknown code the way CPS's
                               does -- confirmed by tabulation, no value 9 exists)
  BPL         -> BPL         (same variable name, same first-digit-of-100 logic:
                               <100 = US state, [100,900) = territory or foreign
                               country, >=900 = not reported (only 133/3.24M in a
                               spot check -- negligible, dropped like CPS's "8"/"9"
                               BPL prefix))
  HFLAG fix   -> (not needed) (ACS has no CPS-ASEC-style split-sample year)

Acs2000.dta has no MIGPLAC1 at all -- IPUMS confirms the "residence 1 year
ago" question wasn't fielded in that transitional year (the standard rolling
ACS design, and this question, start in 2001). Still used for the stock/wage
panel (its own labor stock/wage cross-section doesn't need MIGPLAC1, and
Year=2000 stock is needed as the Origin-at-t denominator for Year=2001 flows)
-- just skipped for the flow panel itself.

No stock-source mismatch risk here (unlike the original StateAnalysis.dta bug):
stocks and wages are built from this same ACS extract, exactly as MakePi.do
built them entirely from CPS rather than merging in StateAnalysis.dta.
================================================================*/

/******************* BUILD BOTH PANELS, ALL ACS YEARS *************************/
frame create FlowPanel
frame create StockWagePanel

loc files: dir "${Data}/acs/" files "Acs*.dta"
foreach f in `files' {

    use "${Data}/acs/`f'", clear

    * Sample restrictions -- shared by both the flow and stock/wage constructions below
    drop if AGE < 16
    * full-time workers only. ACS's UHRSWORK NIU sentinel is 0 (not CPS's 0-then-999
    * convention shift), and 0 < 35 already excludes it -- no second clause needed.
    drop if UHRSWORK < 35

    * Domestic/foreign-born flag from birthplace only, same cutoff logic as MakePi.do:
    * BPL < 100 is a US state (matches STATEFIP numbering), [100,900) is a territory or
    * foreign country, >= 900 is not reported (drop, negligible share of the sample)
    drop if BPL >= 900
    gen Domestic = BPL < 100

    tostring STATEFIP, replace
    replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

    /*---------------- STOCK + WAGE CROSS SECTION (current-state snapshot) ----------*/
    preserve
        replace INCWAGE = . if INCWAGE == 999999
        egen incwage_censor = pctile(INCWAGE), p(99) by(STATEFIP)
        replace INCWAGE = incwage_censor if INCWAGE >= incwage_censor

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

    /*---------------- FLOW (Origin -> Destination) ----------------------------------*/
    capture confirm variable MIGPLAC1
    if !_rc {
        keep PERWT STATEFIP YEAR Domestic MIGPLAC1

        * Interstate only: drop anyone whose prior-year residence was abroad or a US
        * territory (MIGPLAC1 >= 100, l not in S) -- mirrors MakePi.do's MIGSTA1==91 drop
        drop if MIGPLAC1 >= 100
        gen Origin = STATEFIP if MIGPLAC1 == 0                        // did not move: origin = current state
        replace Origin = string(MIGPLAC1, "%02.0f") if mi(Origin)     // moved: origin = reported prior state
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

}

frame change StockWagePanel
frame drop default
tempfile StockWageFile
save `StockWageFile'

frame change FlowPanel
frame drop StockWagePanel

/******************* MERGE IN ORIGIN STOCK (denominator, dated t = Year-1) *******/
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

/******************* MERGE IN ORIGIN'S STOCK AT t+1 (cross-origin ratio's stock term) */
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

/******************* MERGE IN ORIGIN AND DESTINATION WAGES AND DESTINATION STOCK (dated t+1 = Year) ****/
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

/******************* CLEAN UP AND SAVE *************************/
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
