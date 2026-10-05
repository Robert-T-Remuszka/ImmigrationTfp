clear all
do Globals

/*================================================================
Builds the state-year panel that MakeIrf.do estimates local projections
on, from three sources:
  - StateAggSupplyAcs.csv    (AggSupply_Estimate.jl's rho-estimation panel):
                              Z, L, lambda, plus its own Supply/Wage columns
  - StateWageSupplyPanel.dta (MakeAcsAnalysisFiles.do): Supply/Wage levels,
                              Inflow_Domestic/Foreign, and per-BPL-code
                              Inflow<BPL>/Stock<BPL> arrivals and stocks
  - EnclaveSharesBpl.dta     (MakeEnclaves.do): 1960-2000 decennial-census
                              migrant enclave counts by state x BPL code

Enclave counts are reshaped from one-row-per-state-x-BPL to one-row-per-state
(wide across BPL codes, matching the Inflow<BPL>/Stock<BPL> naming convention)
and merged in as m:1 on STATEFIP -- each state's historical settlement
pattern is fixed, so the same values broadcast across all 24 years. Both
population-based (count_bpl_pop_*) and full-time-worker-based
(count_bpl_ft_*, available 1980 onward only) counts are kept for every
available decade.

Also builds the Bartik/shift-share instrument, Bartik_1990: the 79 foreign-
nationality (BPL) groups are first aggregated into 11 world regions (BPL's own
"ns"/"nec" sub-blocks already group most small-sending countries, so this
isn't an arbitrary aggregation on top -- see the region_* locals below for the
exact BPL-to-region mapping, built from usa_00408.xml's BPL value labels).
Regions, not the 79 raw BPL codes, are used because the exact Rotemberg-weight
decomposition needs one instrument column per group x year, and 79 groups x
~20 years would exceed the ~1,173 usable state-year observations (11 regions
x ~20 years does not).

For each region, its full-time-employment share of state l's FOREIGN-BORN
workforce (RegionEmp<year>_<region>/TotalForeignEmp<year>, base year 1990 or
2000) times its NATIONAL growth rate (RegionInflow_{r,t}/L.RegionStock_{r,t}
-- national, not state-level, so the zero-national-stock case the shares'
small-cell counterpart would risk is not expected to bind), summed across
regions. Shares are of the foreign-born stock specifically (not total
employment across all birthplaces), so the 11 regions' shares sum to 1 by
construction -- this matches fg's own foreign-stock normalization (below) and
gives year fixed effects the complete-shares property needed to isolate
shock variation (Borusyak-Hull-Jaravel 2022, sec. 4.2-4.3). The growth rate is
lagged because ACS's MIGPLAC1 (source of Inflow<BPL>) is backward-looking: a
flow dated Year=t is the move between t-1 and t, so it must be scaled by the
stock as of t-1, not the contemporaneous (post-flow) stock at t. The
individual per-region share and growth-rate columns are kept (not collapsed
away) so a Rotemberg-weight decomposition can be run on them.

Alongside these, a new-arrivals version is built from NewArr<BPL> (arrivals
from abroad only, excluding foreign-born interstate movers): BartikNew_1990
and BartikNew_2000 (plus BartikNewLoo_1990/BartikNewLoo_2000, whose national
growth rates leave out the state's own arrivals and stock) use the same shares with national growth rates counting only arrivals to the
US as a whole (NatGrowthNew_*), and fgNew scales the state's arrivals from
abroad by its lagged foreign-born stock. Bartik_1990/fg are unchanged.

No population/employment weighting is used anywhere in this file, per
Chodorow-Reich (2020, "Regional Data in Macroeconomics: Some Advice for
Practitioners", NBER WP 26501, sec. 5): with only 51 cross-sectional units,
his own state-level-calibrated Monte Carlo shows population weighting can
increase both bias and variance in exactly this kind of small-N IV setting.

Saves data/StateAnalysisPanel.dta.
================================================================*/

/******************* Z, L, LAMBDA (FROM THE RHO ESTIMATION PANEL) ************/
* Only Z/L/lambda are kept -- StateAggSupplyAcs.csv's own Supply_*/Wage_* columns
* are just a passthrough of StateWageSupplyPanel.dta's, already in AcsPanel below.
import delimited "${Data}/StateAggSupplyAcs.csv", clear stringcols(1) case(preserve)
ren statefip STATEFIP
ren year Year
keep STATEFIP Year Z L lambda
tempfile ZLambda
save `ZLambda'

/******************* ACS STATE PANEL (SUPPLY, WAGES, INFLOW<BPL>, STOCK<BPL>) */
use "${Data}/StateWageSupplyPanel.dta", clear
merge 1:1 STATEFIP Year using `ZLambda', nogen

tempfile AcsPanel
save `AcsPanel'

/******************* ENCLAVE SHARES: LONG-BY-BPL -> WIDE-BY-BPL **************/
use "${Data}/EnclaveSharesBpl.dta", clear

loc enclavestubs "count_bpl_pop_1960 count_bpl_pop_1970 count_bpl_pop_1980 count_bpl_pop_1990 count_bpl_pop_2000 count_bpl_ft_1980 count_bpl_ft_1990 count_bpl_ft_2000"

reshape wide `enclavestubs', i(STATEFIP) j(BPL)

foreach stub of local enclavestubs {
    ds `stub'*
    foreach v of varlist `r(varlist)' {
        replace `v' = 0 if mi(`v')
    }
}

tempfile EnclaveWide
save `EnclaveWide'

/******************* MERGE *****************************************************/
use `AcsPanel', clear
merge m:1 STATEFIP using `EnclaveWide', nogen

/******************* BPL -> WORLD REGION CROSSWALK (11 REGIONS) ***************/
* Built directly from usa_00408.xml's BPL value labels (real country names).
loc region_Mexico          "200"
loc region_CentAmCarib     "210 250 260 299"
loc region_SouthAmerica    "300"
loc region_CanadaOceania   "150 700 710 100 105 110 115 160 199"
loc region_NWEurope        "400 401 402 404 405 410 411 412 413 414 420 421 425 426 429"
loc region_SEurope         "430 433 434 436 438"
loc region_EEurope         "450 451 452 453 454 455 456 457 460 461 462 465 499"
loc region_EAsia           "500 501 502 509"
loc region_SEAsia          "511 512 513 514 515 516 517 518 519"
loc region_SAsia           "520 521 522 524"
loc region_MidEast         "531 532 534 535 536 537 540 541 542 543 544 548 549 599"
loc region_Africa          "600"
loc regions "Mexico CentAmCarib SouthAmerica CanadaOceania NWEurope SEurope EEurope EAsia SEAsia SAsia MidEast Africa"

/******************* NATIONAL REGION GROWTH RATES (THE BARTIK "SHIFT") ********/
loc inflowdomfor "Inflow_Domestic Inflow_Foreign"
ds Inflow*
loc allinflowvars `r(varlist)'
loc bplinflowvars: list allinflowvars - inflowdomfor

* NewArr<BPL>: arrivals from abroad only (no foreign-born interstate movers)
ds NewArr*
loc allnewvars `r(varlist)'
loc newarrtot "NewArr_Foreign"
loc bplnewvars: list allnewvars - newarrtot

ds Stock*
loc bplstockvars `r(varlist)'

* A handful of BPL codes have a Stock<code> column but no Inflow<code> column
* -- that nationality group had zero recorded ACS arrivals across the entire
* 2001-2024 sample (too rare to ever appear in the arrivals collapse), so the
* wide reshape never created the column. Fill in an explicit all-zero series
* rather than silently dropping the group from the growth-rate calculation.
loc stockcodes ""
foreach v of local bplstockvars {
    loc stockcodes "`stockcodes' `=substr("`v'", 6, .)'"
}
loc inflowcodes ""
foreach v of local bplinflowvars {
    loc inflowcodes "`inflowcodes' `=substr("`v'", 7, .)'"
}
loc missinginflow: list stockcodes - inflowcodes
foreach code of local missinginflow {
    gen double Inflow`code' = 0
    loc bplinflowvars "`bplinflowvars' Inflow`code'"
}

loc newcodes ""
foreach v of local bplnewvars {
    loc newcodes "`newcodes' `=substr("`v'", 7, .)'"
}
loc missingnew: list stockcodes - newcodes
foreach code of local missingnew {
    gen double NewArr`code' = 0
    loc bplnewvars "`bplnewvars' NewArr`code'"
}

* Not every BPL code in the region crosswalk actually appears among the 79
* realized Stock<BPL> codes (e.g. 519 "Southeast Asia, ns" is a valid IPUMS
* code but never realized in this ACS sample) -- intersect each region's list
* with the codes that actually exist before referencing Stock<code>/Inflow<code>.
foreach r of local regions {
    loc region_`r': list region_`r' & stockcodes
}

preserve
    keep Year `bplinflowvars' `bplstockvars' `bplnewvars'
    collapse (sum) `bplinflowvars' `bplstockvars' `bplnewvars', by(Year)
    tsset Year

    foreach r of local regions {
        loc inflowvars ""
        loc stockvars ""
        loc newvars ""
        foreach code of local region_`r' {
            loc inflowvars "`inflowvars' Inflow`code'"
            loc stockvars  "`stockvars' Stock`code'"
            loc newvars    "`newvars' NewArr`code'"
        }
        egen double RegionInflow_`r' = rowtotal(`inflowvars')
        egen double RegionStock_`r'  = rowtotal(`stockvars')
        egen double RegionNewArr_`r' = rowtotal(`newvars')
        gen double NatGrowthRegion_`r' = RegionInflow_`r' / L.RegionStock_`r'
        gen double NatGrowthNew_`r'    = RegionNewArr_`r' / L.RegionStock_`r'
    }

    keep Year NatGrowthRegion_* NatGrowthNew_*
    tempfile NatGrowth
    save `NatGrowth'
restore

merge m:1 Year using `NatGrowth', nogen

/******************* BARTIK INSTRUMENT (1990 SHARES OF THE FOREIGN-BORN STOCK) */
* Share denominator is the state's total FOREIGN-born FT employment (sum of
* the 11 regions' own counts), not total employment across all birthplaces --
* this makes the shares complete (sum to 1 by construction, verified below),
* matching fg's own foreign-stock normalization.
loc bartikterms1990 ""
foreach r of local regions {
    loc ftvars ""
    foreach code of local region_`r' {
        loc ftvars "`ftvars' count_bpl_ft_1990`code'"
    }
    egen double RegionEmp1990_`r' = rowtotal(`ftvars')
    loc bartikterms1990 "`bartikterms1990' RegionEmp1990_`r'"
}
egen double TotalForeignEmp1990 = rowtotal(`bartikterms1990')

loc bartikterms1990 ""
foreach r of local regions {
    gen double RegionShare1990_`r' = RegionEmp1990_`r' / TotalForeignEmp1990
    gen double Bartik_1990_`r' = RegionShare1990_`r' * NatGrowthRegion_`r'
    loc bartikterms1990 "`bartikterms1990' Bartik_1990_`r'"

    la var RegionEmp1990_`r'   "1990 FT-employment count, `r' birthplaces, state l"
    la var RegionShare1990_`r' "1990 FT-employment share of `r' birthplaces among ALL foreign-born in state l (RegionEmp1990_`r'/TotalForeignEmp1990); these 11 shares sum to 1 by construction"
    la var NatGrowthRegion_`r' "National ACS growth rate, `r' birthplaces (RegionInflow_t/L.RegionStock_t)"
    la var NatGrowthNew_`r'    "National growth rate from new arrivals from abroad only, `r' birthplaces (RegionNewArr_t/L.RegionStock_t)"
    la var Bartik_1990_`r'     "`r' region's contribution to Bartik_1990 (RegionShare1990_`r' x NatGrowthRegion_`r')"
}
egen double Bartik_1990 = rowtotal(`bartikterms1990'), missing

la var TotalForeignEmp1990 "Total 1990 full-time-worker population, foreign-born only, state l (Bartik_1990 share denominator)"
la var Bartik_1990 "Bartik/shift-share instrument (complete shares): sum over 11 world regions of (1990 FT-employment share of state l's FOREIGN-born workforce) x (national growth rate of the region's ACS stock, aggregate BPL flow/aggregate BPL stock); missing for Year=2001 (no lag)"

/******************* BARTIK, NEW-ARRIVALS SHIFT (1990 SHARES) ****************/
* Same 1990 shares as Bartik_1990, but the national growth rate counts only
* arrivals from abroad to the US as a whole (NatGrowthNew_*), not the sum of
* state inflows -- which also counts foreign-born movers between states, who
* leave the national stock unchanged. Kept alongside Bartik_1990 for comparison.
loc bartiknew1990 ""
foreach r of local regions {
    gen double BartikNew_1990_`r' = RegionShare1990_`r' * NatGrowthNew_`r'
    loc bartiknew1990 "`bartiknew1990' BartikNew_1990_`r'"
    la var BartikNew_1990_`r' "`r' region's contribution to BartikNew_1990 (RegionShare1990_`r' x NatGrowthNew_`r')"
}
egen double BartikNew_1990 = rowtotal(`bartiknew1990'), missing
la var BartikNew_1990 "Bartik with 1990 foreign-born shares (complete) x national growth of arrivals from abroad only, summed over 11 regions; missing for Year=2001 (no lag)"

/******************* BARTIK INSTRUMENT (2000 SHARES OF THE FOREIGN-BORN STOCK) */
* Same construction as Bartik_1990 above, with the base-year shares moved from
* 1990 to 2000 -- the national growth rates (NatGrowthRegion_*) are unchanged,
* since they come from the ACS panel (2001-2024), not the enclave shares.
loc bartikterms2000 ""
foreach r of local regions {
    loc ftvars ""
    foreach code of local region_`r' {
        loc ftvars "`ftvars' count_bpl_ft_2000`code'"
    }
    egen double RegionEmp2000_`r' = rowtotal(`ftvars')
    loc bartikterms2000 "`bartikterms2000' RegionEmp2000_`r'"
}
egen double TotalForeignEmp2000 = rowtotal(`bartikterms2000')

loc bartikterms2000 ""
foreach r of local regions {
    gen double RegionShare2000_`r' = RegionEmp2000_`r' / TotalForeignEmp2000
    gen double Bartik_2000_`r' = RegionShare2000_`r' * NatGrowthRegion_`r'
    loc bartikterms2000 "`bartikterms2000' Bartik_2000_`r'"

    la var RegionEmp2000_`r'   "2000 FT-employment count, `r' birthplaces, state l"
    la var RegionShare2000_`r' "2000 FT-employment share of `r' birthplaces among ALL foreign-born in state l (RegionEmp2000_`r'/TotalForeignEmp2000); these 11 shares sum to 1 by construction"
    la var Bartik_2000_`r'     "`r' region's contribution to Bartik_2000 (RegionShare2000_`r' x NatGrowthRegion_`r')"
}
egen double Bartik_2000 = rowtotal(`bartikterms2000'), missing

la var TotalForeignEmp2000 "Total 2000 full-time-worker population, foreign-born only, state l (Bartik_2000 share denominator)"
la var Bartik_2000 "Bartik/shift-share instrument (complete shares): sum over 11 world regions of (2000 FT-employment share of state l's FOREIGN-born workforce) x (national growth rate of the region's ACS stock, aggregate BPL flow/aggregate BPL stock); missing for Year=2001 (no lag)"

/******************* BARTIK, NEW-ARRIVALS SHIFT (2000 SHARES) ****************/
* Same as BartikNew_1990 above with the base-year shares moved to 2000.
loc bartiknew2000 ""
foreach r of local regions {
    gen double BartikNew_2000_`r' = RegionShare2000_`r' * NatGrowthNew_`r'
    loc bartiknew2000 "`bartiknew2000' BartikNew_2000_`r'"
    la var BartikNew_2000_`r' "`r' region's contribution to BartikNew_2000 (RegionShare2000_`r' x NatGrowthNew_`r')"
}
egen double BartikNew_2000 = rowtotal(`bartiknew2000'), missing
la var BartikNew_2000 "Bartik with 2000 foreign-born shares (complete) x national growth of arrivals from abroad only, summed over 11 regions; missing for Year=2001 (no lag)"

/******************* LEAVE-ONE-OUT NEW-ARRIVALS BARTIKS ***********************/
* For state l, the national growth rate of each region drops l's own arrivals
* from the numerator and l's own lagged stock from the denominator, so a large
* state does not contribute to the shift it is exposed to. Same 1990 / 2000
* shares as BartikNew_1990 / BartikNew_2000.
sort STATEFIP Year
loc bartiknewloo1990 ""
loc bartiknewloo2000 ""
foreach r of local regions {
    loc newvars ""
    loc stockvars ""
    foreach code of local region_`r' {
        loc newvars   "`newvars' NewArr`code'"
        loc stockvars "`stockvars' Stock`code'"
    }
    egen double RegionNewArrSt_`r' = rowtotal(`newvars')
    egen double RegionStockSt_`r'  = rowtotal(`stockvars')
    bys Year: egen double NatNewArr_`r' = total(RegionNewArrSt_`r')
    bys Year: egen double NatStock_`r'  = total(RegionStockSt_`r')
    bys STATEFIP (Year): gen double NatGrowthNewLoo_`r' = (NatNewArr_`r' - RegionNewArrSt_`r') / (NatStock_`r'[_n-1] - RegionStockSt_`r'[_n-1])
    la var NatGrowthNewLoo_`r' "Leave-one-out national growth rate from new arrivals from abroad, `r' birthplaces: (national arrivals - state l's) / (lagged national stock - state l's lagged stock)"

    gen double BartikNewLoo_1990_`r' = RegionShare1990_`r' * NatGrowthNewLoo_`r'
    gen double BartikNewLoo_2000_`r' = RegionShare2000_`r' * NatGrowthNewLoo_`r'
    la var BartikNewLoo_1990_`r' "`r' region's contribution to BartikNewLoo_1990 (RegionShare1990_`r' x NatGrowthNewLoo_`r')"
    la var BartikNewLoo_2000_`r' "`r' region's contribution to BartikNewLoo_2000 (RegionShare2000_`r' x NatGrowthNewLoo_`r')"
    loc bartiknewloo1990 "`bartiknewloo1990' BartikNewLoo_1990_`r'"
    loc bartiknewloo2000 "`bartiknewloo2000' BartikNewLoo_2000_`r'"

    drop RegionNewArrSt_`r' RegionStockSt_`r' NatNewArr_`r' NatStock_`r'
}
egen double BartikNewLoo_1990 = rowtotal(`bartiknewloo1990'), missing
egen double BartikNewLoo_2000 = rowtotal(`bartiknewloo2000'), missing
la var BartikNewLoo_1990 "BartikNew_1990 with leave-one-out national growth rates (state l's own arrivals and stock removed); missing for Year=2001 (no lag)"
la var BartikNewLoo_2000 "BartikNew_2000 with leave-one-out national growth rates (state l's own arrivals and stock removed); missing for Year=2001 (no lag)"

/******************* MIGRATION SHOCK (FOREIGN-STOCK-NORMALIZED) ***************/
* fg = Inflow_Foreign scaled by the state's own lagged FOREIGN-BORN stock (sum
* of the Stock<BPL> columns) -- an actual percent-growth rate in the local
* foreign-born population, matching Bartik_v2's own foreign-only universe and
* closer to the elasticity object the structural model's parameters are
* defined on. (Superseded the original total-employment-normalized version,
* kept below as fg_TotalEmp for reference/robustness.)
sort STATEFIP Year
ds Stock*
egen double StockForeignTot = rowtotal(`r(varlist)')
bys STATEFIP (Year): gen double fg = Inflow_Foreign / StockForeignTot[_n-1]

la var StockForeignTot "State l's total foreign-born stock (sum of Stock<BPL>), lagged for use as fg's denominator"
la var fg "Migration shock: Inflow_Foreign_t / StockForeignTot_{t-1} -- percent growth in the local foreign-born stock; missing for Year=2001 (no lag)"

/******************* MIGRATION SHOCK (ARRIVALS FROM ABROAD ONLY) **************/
bys STATEFIP (Year): gen double fgNew = NewArr_Foreign / StockForeignTot[_n-1]
la var fgNew "Migration shock, arrivals from abroad only: NewArr_Foreign_t / StockForeignTot_{t-1}; excludes foreign-born interstate movers; missing for Year=2001 (no lag)"

/******************* MIGRATION SHOCK (EMPLOYMENT-DENOMINATED, LAGGED) *********/
* fg_TotalEmp = the original definition, scaled by the state's total
* employment as of t-1 instead of its foreign-born stock -- kept for
* comparison/robustness, no longer the default.
gen double TotalEmp = Supply_Domestic + Supply_Foreign
bys STATEFIP (Year): gen double fg_TotalEmp = Inflow_Foreign / TotalEmp[_n-1]

la var TotalEmp "State l's total (domestic + foreign) FTFY labor supply"
la var fg_TotalEmp "Migration shock (employment-denominated): Inflow_Foreign_t / TotalEmp_{t-1}; missing for Year=2001 (no lag)"

sort STATEFIP Year
save "${Data}/StateAnalysisPanel.dta", replace
di as text "Wrote ${Data}/StateAnalysisPanel.dta"
