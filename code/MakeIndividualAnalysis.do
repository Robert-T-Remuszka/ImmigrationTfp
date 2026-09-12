clear all
do Globals

/*================================================================
Builds the individual-level (person-year) analysis file for the "Choice Probabilities" probit 

SAMPLE: ACS 2000-2022

SAMPLE RESTRICTIONS -- deliberately IDENTICAL to MakeAcsPi.do's, so this
individual-level file's sample matches the state-level panel, (AGE>=16, full-time workers only,
BPL<900, same Domestic/Foreign definition).

TASK INDEX (tau): rebuilt directly from Peri and Sparber (2009)'s own
replication code (ICPSR project 113564, AEJApp-2008-0057.R3_data_codes/
Extraction/TempONetVals.do) -- their code
defines, on Abilities_occ1990.dta (52 O*NET Ability importance scores,
already percentile-rescaled by PS against the 2000 Census workforce and
keyed on occ1990):
    manual   = mean(im_a22..im_a40)   [19 "movement and strength" abilities]
    language = mean(im_a1..im_a4)     [4 oral/written expression/comprehension]
    c/m      = language / manual      [PS's basic communication/manual index]
tau is then THIS project's own construction, not PS's: the employment-weighted
percentile rank of c/m across occupations (weight = total ACS employment in
that occupation, Domestic+Foreign pooled, summed across ALL sample years) --
required because the model's z(tau) is time- and state-invariant (Section 3:
"the comparative advantage schedule does not vary across states"), so tau
is a single fixed ranking, not recomputed by state or year. This is a
RANK-based percentile

w_tilde: ln(Wage_Domestic/Wage_Foreign) at the state-year level, using the
IDENTICAL wage-cleaning procedure as MakeAcsPi.do's stock/wage cross-section
(topcode at 999999, winsorize top 1% by state, weight by PERWT) -- recomputed
here rather than loaded from a saved file since MakeAcsPi.do doesn't persist
that intermediate cross-section on its own.
================================================================*/

/******************* BUILD THE TASK INDEX FROM PERI-SPARBER'S DATA **********/
use "${PeriSparber}/Abilities_occ1990.dta", clear
gen manual   = (1/19)*(im_a22+im_a23+im_a24+im_a25+im_a26+im_a27+im_a28+im_a29 ///
                       +im_a30+im_a31+im_a32+im_a33+im_a34+im_a35+im_a36+im_a37 ///
                       +im_a38+im_a39+im_a40)
gen language = (1/4)*(im_a1+im_a2+im_a3+im_a4)
gen cm_ratio = language / manual
keep occ1990 manual language cm_ratio
ren occ1990 OCC1990
la var cm_ratio "Peri-Sparber basic communication/manual ratio (language/manual), from their own ICPSR replication code"
count
tempfile TaskAbility
save `TaskAbility'

/******************* BUILD BOTH PANELS (INDIVIDUAL ROWS + STATE-YEAR WAGES), ONE PASS ***/
frame create IndividualPanel
frame create StateWagePanel

loc files: dir "${Data}/acs/" files "Acs*.dta"
foreach f in `files' {

    use "${Data}/acs/`f'", clear

    * Sample restrictions -- IDENTICAL to MakeAcsPi.do's (see its header for the
    * ACS-specific coding notes these restrictions rely on).
    drop if AGE < 16
    drop if UHRSWORK < 35
    drop if BPL >= 900
    gen Domestic = BPL < 100
    gen Foreign  = 1 - Domestic

    /*---------------- STATE-YEAR WAGE CROSS SECTION (same cleaning as MakeAcsPi.do) --*/
    preserve
        replace INCWAGE = . if INCWAGE == 999999
        egen incwage_censor = pctile(INCWAGE), p(99) by(STATEFIP)
        replace INCWAGE = incwage_censor if INCWAGE >= incwage_censor

        collapse (mean) Wage = INCWAGE [pw = PERWT], by(STATEFIP YEAR Domestic)
        reshape wide Wage, i(STATEFIP YEAR) j(Domestic)
        ren Wage0 Wage_Foreign
        ren Wage1 Wage_Domestic
        ren YEAR Year

        frame StateWagePanel: xframeappend default
    restore

    /*---------------- INDIVIDUAL ROWS (for the task-occupation side) ----------------*/
    keep STATEFIP YEAR OCC1990 Foreign PERWT
    ren YEAR Year

    frame IndividualPanel: xframeappend default
}

/******************* STATE-YEAR WAGE RATIO (w_tilde) ****************************/
frame change StateWagePanel
frame drop default
gen w_tilde = ln(Wage_Domestic / Wage_Foreign)
la var w_tilde "ln(Wage_Domestic/Wage_Foreign), state-year, ACS, PERWT-weighted mean, top-1%-by-state winsorized (NOT composition-adjusted yet)"
keep STATEFIP Year w_tilde
tempfile StateWageFile
save `StateWageFile'

/******************* ASSEMBLE THE INDIVIDUAL ANALYSIS FILE ***********************/
frame change IndividualPanel

merge m:1 STATEFIP Year using `StateWageFile', keep(1 3) nogen
* Occupation-level task-ability data is only available for civilian occupied
* workers (339 occ1990 categories) -- military/unclassified codes won't match.
merge m:1 OCC1990 using `TaskAbility', keep(1 3)
di as text "--- OCC1990 merge with Peri-Sparber's occupation-ability index ---"
tab _merge
count if _merge == 1
drop if _merge == 1
drop _merge

/******************* TASK PERCENTILE (tau): EMPLOYMENT-WEIGHTED, POOLED ACROSS YEARS **/
preserve
    collapse (rawsum) TotalEmp = PERWT (mean) cm_ratio, by(OCC1990)
    * z(tau) is assumed time- and state-invariant (Section 3), so tau is a single,
    * fixed national ranking -- not recomputed by state or year.
    gsort cm_ratio
    gen CumEmp = sum(TotalEmp)
    qui summ CumEmp, meanonly
    loc grand = r(max)
    * Midpoint percentile: centers each occupation's mass rather than placing it
    * at the top of its own cumulative bucket.
    gen tau = (CumEmp - 0.5*TotalEmp) / `grand'
    gen phi_inv_tau = invnormal(tau)
    la var tau         "Employment-weighted percentile rank of Peri-Sparber's c/m ratio across occ1990 categories, pooled over all sample years"
    la var phi_inv_tau "Phi^-1(tau) -- the model's task-index regressor"
    keep OCC1990 tau phi_inv_tau cm_ratio TotalEmp
    tempfile TauLookup
    save `TauLookup'
restore

merge m:1 OCC1990 using `TauLookup', keep(1 3) nogen

/******************* CLEAN UP AND SAVE *************************/
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
