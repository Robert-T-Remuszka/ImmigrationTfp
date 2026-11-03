clear all
do Globals

/*================================================================
Builds the ROW-origin ACS flow counts used in the ROW calibration, all from
the ACS 2016 wave reporting MIGPLAC1 >= 100 (living abroad in 2015 -- the
same origin year as every other piece of this calibration):
  - immigration_count_2015: new foreign-born immigrants, the entry-rate
    target's numerator;
  - return_migration_count_2015: returning domestic-born natives;
  - adult_immigration_count_2015: new foreign-born immigrants aged 15+ of any
    work status, the population Gallup's desire-to-migrate count covers. It
    is the numerator of the admission probability ξ.

The first two use the same sample restrictions as AcsPiPanel.dta's
construction (age 16+, full-time workers, BPL<900), but keep exactly the rows
that panel drops: MIGPLAC1 >= 100 is the ROW-origin counterpart to its
interstate-only panel. Residents of DC are dropped throughout, since the
model's US locations are the 50 states.

Input: data/acs/Acs2016.dta.
Output: data/RowFlowCounts.dta, one row.
================================================================*/

use "${Data}/acs/Acs2016.dta", clear

drop if BPL >= 900
drop if STATEFIP == 11
gen Domestic = BPL < 100

* ROW-origin only: MIGPLAC1 >= 100 means living abroad (or in a US
* territory) a year ago, i.e. in 2015.
keep if MIGPLAC1 >= 100

qui summ PERWT if Domestic == 0 & AGE >= 15, meanonly
loc adult_immigration_count = r(sum)

drop if AGE < 16
drop if UHRSWORK < 35

qui summ PERWT if Domestic == 0, meanonly
loc immigration_count = r(sum)
qui summ PERWT if Domestic == 1, meanonly
loc return_migration_count = r(sum)

di as text "Immigration count (foreign-born, abroad in 2015): " as result `immigration_count'
di as text "Return migration count (domestic-born, abroad in 2015): " as result `return_migration_count'
di as text "Adult immigration count (foreign-born, 15+, any work status): " as result `adult_immigration_count'

clear
set obs 1
gen double immigration_count_2015 = `immigration_count'
gen double return_migration_count_2015 = `return_migration_count'
gen double adult_immigration_count_2015 = `adult_immigration_count'
la var immigration_count_2015 "Weighted count of foreign-born ACS 2016 respondents who lived abroad in 2015 (new immigrants)"
la var return_migration_count_2015 "Weighted count of domestic-born ACS 2016 respondents who lived abroad in 2015 (returning natives)"
la var adult_immigration_count_2015 "Weighted count of foreign-born ACS 2016 respondents aged 15+, any work status, who lived abroad in 2015"

save "${Data}/RowFlowCounts.dta", replace
di as text "Wrote ${Data}/RowFlowCounts.dta"
