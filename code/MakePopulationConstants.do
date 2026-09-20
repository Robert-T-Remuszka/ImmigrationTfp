clear all
do Globals

/*================================================================
Builds the external population constants needed to scale the steady-state
labor-supply normalization (see SsSolve_Functions.jl's load_total_population):
US-born people living outside the US, and total world population, both 2015.

Input: data/UnEstimates/undesa_pd_2024_ims_stock_by_sex_destination_and_origin.xlsx
(UN International Migrant Stock 2024 revision, Table 1) and
data/UnEstimates/WPP2024_Demographic_Indicators_Medium.csv (UN World
Population Prospects 2024 revision)
Output: data/PopulationConstants.dta, one row (us_abroad_2015, world_pop_2015).
================================================================*/

/******************* US-BORN PEOPLE LIVING ABROAD, 2015 ***************************/
* import excel with no header row names variables by Excel column letter
* ("2015" appears three times, once per sex, so positional columns are more
* robust than relying on import excel to de-duplicate text headers): F origin
* name, E destination location code, M both-sexes-combined stock in 2015.
import excel "${Data}/UnEstimates/undesa_pd_2024_ims_stock_by_sex_destination_and_origin.xlsx", ///
    sheet("Table 1") cellrange(A11) clear
drop in 1
destring E M, replace force

* Real countries only (location codes 900 = World and 1829-1836 = UN regions
* would otherwise double-count); origin = United States excludes the (absent)
* USA-to-USA row automatically.
keep if F == "United States of America*" & E < 900 & !missing(E)

qui summ M, meanonly
loc us_abroad = r(sum)
di as text "US-born people living abroad, 2015: " as result `us_abroad'

/******************* WORLD POPULATION, 2015 ****************************************/
import delimited "${Data}/UnEstimates/WPP2024_Demographic_Indicators_Medium.csv", clear
keep if location == "World" & time == 2015

* TPopulation1July is in thousands (standard UN WPP convention).
qui summ tpopulation1july, meanonly
loc world_pop = r(mean) * 1000

di as text "World population, 2015: " as result `world_pop'

/******************* SAVE ***********************************************************/
clear
set obs 1
gen double us_abroad_2015 = `us_abroad'
gen double world_pop_2015 = `world_pop'
la var us_abroad_2015 "US-born people living outside the US, 2015 (UN Intl Migrant Stock 2024 rev., Table 1)"
la var world_pop_2015  "World population, 2015 (UN World Population Prospects 2024 rev., mid-year)"

save "${Data}/PopulationConstants.dta", replace
di as text "Wrote ${Data}/PopulationConstants.dta"
