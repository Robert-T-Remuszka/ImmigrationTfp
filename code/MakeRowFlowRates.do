clear all
do Globals

/*================================================================
Builds the two ROW-origin ACS flow counts needed to identify the asymmetric
ROW migration costs: the number of new foreign-born immigrants and the
number of returning domestic-born natives, both observed in the ACS 2016
wave reporting MIGPLAC1 >= 100 (living abroad in 2015 -- the same origin
year as every other piece of this calibration).

Same sample restrictions as AcsPiPanel.dta's construction (age 16+, full-time
workers, BPL<900), but keeping exactly the rows that panel drops: MIGPLAC1
>= 100 is the ROW-origin counterpart to its interstate-only panel.

Input: data/acs/Acs2016.dta.
Output: data/RowFlowCounts.dta, one row (immigration_count_2015,
return_migration_count_2015).
================================================================*/

use "${Data}/acs/Acs2016.dta", clear

drop if AGE < 16
drop if UHRSWORK < 35
drop if BPL >= 900
gen Domestic = BPL < 100

* ROW-origin only: MIGPLAC1 >= 100 means living abroad (or in a US
* territory) a year ago, i.e. in 2015.
keep if MIGPLAC1 >= 100

collapse (sum) PERWT, by(Domestic)

qui summ PERWT if Domestic == 0, meanonly
loc immigration_count = r(sum)
qui summ PERWT if Domestic == 1, meanonly
loc return_migration_count = r(sum)

di as text "Immigration count (foreign-born, abroad in 2015): " as result `immigration_count'
di as text "Return migration count (domestic-born, abroad in 2015): " as result `return_migration_count'

clear
set obs 1
gen double immigration_count_2015 = `immigration_count'
gen double return_migration_count_2015 = `return_migration_count'
la var immigration_count_2015 "Weighted count of foreign-born ACS 2016 respondents who lived abroad in 2015 (new immigrants)"
la var return_migration_count_2015 "Weighted count of domestic-born ACS 2016 respondents who lived abroad in 2015 (returning natives)"

save "${Data}/RowFlowCounts.dta", replace
di as text "Wrote ${Data}/RowFlowCounts.dta"
