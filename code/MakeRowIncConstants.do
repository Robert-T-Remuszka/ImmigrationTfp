clear all
do Globals

/*================================================================
Builds Rest-of-World wage constants, Wᵈ_ROW (Americans-abroad wage) and
Wᶠ_ROW (rest-of-world wage), both in real 2009 USD to match the rest of the
model's real-dollar units (SsSolve_Functions.jl's A calibration).

Inputs (fetched by code/FetchRowIncConstants.py -- re-run that script to
refresh from the live sources):
- data/IrsForeignEarnedIncome.csv: IRS Statistics of Income Table 2
  ("Individual Income Tax Returns With Form 2555: Foreign-Earned Income ...
  by Country or Region"), "All geographic areas" row, tax years 2011 and
  2016 (Total foreign-earned income / Number of returns).
- data/WorldBankGdpLaborForce2015.csv: World Bank GDP (current US$,
  indicator NY.GDP.MKTP.CD) and total labor force (SL.TLF.TOTL.IN), World
  and United States, 2015.
Also uses data/GdpPriceDeflator.csv (already in the repo) to convert nominal
dollars to real 2009 USD.
Output: data/RowIncomeConstants.dta, one row (wd_row_2009usd, wf_row_2009usd).
================================================================*/

/******************* GDP PRICE DEFLATOR LOOKUP *******************/
import delimited "${Data}/GdpPriceDeflator.csv", clear varnames(1)
gen year = real(substr(v1, 1, 4))
foreach y in 2009 2011 2015 2016 {
    qui summ pricedeflator if year == `y', meanonly
    loc def_`y' = r(mean)
}
di as text "GDP price deflator: 2009=`def_2009' 2011=`def_2011' 2015=`def_2015' 2016=`def_2016'"

/******************* Wᵈ_ROW: AMERICANS-ABROAD WAGE ********************************/
* IRS Table 2's per-filer foreign-earned income at 2011 and 2016, each deflated
* to real 2009 USD, then linearly interpolated to 2015 (4/5 of the way from
* 2011 to 2016) since IRS only publishes this table every ~5 years.
import delimited "${Data}/IrsForeignEarnedIncome.csv", clear
sort tax_year
loc wd_2011 = weighted_average_usd[1] * (`def_2009' / `def_2011')
loc wd_2016 = weighted_average_usd[2] * (`def_2009' / `def_2016')
loc wd_row  = `wd_2011' + (`wd_2016' - `wd_2011') * (2015 - 2011) / (2016 - 2011)
di as text "Wᵈ_ROW candidates: 2011=`wd_2011' 2016=`wd_2016' -> 2015 (interpolated) = `wd_row'"

/******************* Wᶠ_ROW: REST-OF-WORLD WAGE ***********************************/
* GDP per labor-force member, World excluding the US, 2015 (current US$),
* deflated to real 2009 USD. Labor force (not population) so this isn't
* diluted by dependency-ratio differences across countries the way GDP per
* capita would be.
import delimited "${Data}/WorldBankGdpLaborForce2015.csv", clear
qui summ value if indicator == "NY.GDP.MKTP.CD" & country == "World", meanonly
loc gdp_world = r(mean)
qui summ value if indicator == "NY.GDP.MKTP.CD" & country == "United States", meanonly
loc gdp_us = r(mean)
qui summ value if indicator == "SL.TLF.TOTL.IN" & country == "World", meanonly
loc lf_world = r(mean)
qui summ value if indicator == "SL.TLF.TOTL.IN" & country == "United States", meanonly
loc lf_us = r(mean)

loc row_gdp_per_worker = (`gdp_world' - `gdp_us') / (`lf_world' - `lf_us')
loc wf_row = `row_gdp_per_worker' * (`def_2009' / `def_2015')
di as text "World GDP=`gdp_world' US GDP=`gdp_us' World LF=`lf_world' US LF=`lf_us'"
di as text "Wᶠ_ROW (nominal 2015 USD) = `row_gdp_per_worker' -> real 2009 USD = `wf_row'"

/******************* SAVE **********************************************************/
clear
set obs 1
gen double wd_row_2009usd = `wd_row'
gen double wf_row_2009usd = `wf_row'
la var wd_row_2009usd "Americans-abroad wage, real 2009 USD (IRS SOI Table 2 foreign-earned income per filer, 2011/2016 interpolated to 2015)"
la var wf_row_2009usd "Rest-of-world wage, real 2009 USD (World Bank GDP per labor-force member, World ex-US, 2015)"

save "${Data}/RowIncomeConstants.dta", replace
di as text "Wrote ${Data}/RowIncomeConstants.dta"
