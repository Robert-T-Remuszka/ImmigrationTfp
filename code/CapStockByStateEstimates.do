clear all
set mem 10g
do Globals

/******************************************************************************
Clean Fixed Asset file - Note figures are in MILLIONS
- The BEA only lists the industry names in the fixed asset tables,
so I had to download a concordance table from them. This is also in the
"$Data/FixedAssets" directory as well.
******************************************************************************/
import delimited "$Data/FixedAssets/NetFixedAssetsTable.csv", clear varnames(1)
drop if _n <= 2
drop if _n == 2
ren v2 IndustryName
drop table31
forvalues i = 3/32 { // Rename with corresponding year
    loc newname = string(v`i'[1])
    ren v`i' FixedAssets_`newname'
}
drop if inlist(_n,1,2) // no longer needed
drop if _n >= 96 // blank rows

/*
GOAL: Construct a consistent time series of the capital stock by state based on
industries reported in Peri (2012):
1. Agriculture, Forestry, Fishing and Hunting
2. Mining 
3. Utilities
4. Construction
5. Manufacturing
6. Wholesale Trade
7. Retail Trade
8. Transportation and Warehousing
9. Information
10. Finance and Insurance
11. Real Estate and Rental and Leasing
12. Professional, Scientific, and Technical Services
13. Management of Companies and Enterprises
14. Administrative and Waste Management Services
15. Educational Services
16. Health Care and Social Assistance
17. Arts, Entertainment, and Recreation
18. Accommodation and Food Services;
19. Other Services, except Government
*/
keep if substr(IndustryName,1,1) != " "

* drop the footnote number in "Management of companies"
replace IndustryName = subinstr(IndustryName," 5", "",1)

* Create NAICS codes
gen Naics = ""
replace Naics = "11" if IndustryName == "Agriculture, forestry, fishing, and hunting"
replace Naics = "21" if IndustryName == "Mining"
replace Naics = "22" if IndustryName == "Utilities"
replace Naics = "23" if IndustryName == "Construction"
replace Naics = "31-33" if IndustryName == "Manufacturing"
replace Naics = "42" if IndustryName == "Wholesale trade"
replace Naics = "44-45" if IndustryName == "Retail trade"
replace Naics = "48-49" if IndustryName == "Transportation and warehousing"
replace Naics = "51" if IndustryName == "Information"
replace Naics = "52" if IndustryName == "Finance and insurance"
replace Naics = "53" if IndustryName == "Real estate and rental and leasing"
replace Naics = "54" if IndustryName == "Professional, scientific, and technical services"
replace Naics = "55" if IndustryName == "Management of companies and enterprises"
replace Naics = "56" if IndustryName == "Administrative and waste management services"
replace Naics = "61" if IndustryName == "Educational services"
replace Naics = "62" if IndustryName == "Health and social assistance"
replace Naics = "71" if IndustryName == "Arts, entertainment, and recreation"
replace Naics = "72" if IndustryName == "Accommodation and food services"
replace Naics = "81" if IndustryName == "Other services, except government"
order IndustryName Naics

reshape long FixedAssets_, i(IndustryName Naics) j(year)
ren FixedAssets_ FixedAssets
tempfile FixedAssetData
save `FixedAssetData'

/******************************************************************************
Clean SAGDP - The GDP industry table from 1997 to 2023
******************************************************************************/
frame create Gdp97to23
frame Gdp97to23 {
    import delimited "$Data/gdp_state_industry_BEA/SAGDP/SAGDP2N__ALL_AREAS_1997_2023.csv", clear
    gen statefip = subinstr(geofips,`"""',"",2)
    replace statefip = substr(statefip,2,2) // there is a leading space for some reason
    drop geofips
    drop if mi(geoname)
    drop if inlist(geoname, "Far West", "Rocky Mountain", "Southwest", "Southeast", ///
    "Plains", "Great Lakes", "Mideast", "New England")
    drop if industryclassification == "..."
    drop if inlist(statefip,"00")
    gen third = substr(industryclassification,3,1)
    gen len = strlen(industryclassification)
    keep if len == 2 | third == "-"
    drop if inlist(industryclassification,"92","31-33,51")
    drop len third unit tablename description linecode region
    order statefip
    ren geoname StateName
    tostring v31, force replace // for some reason this get read in just fine...
    loc year = 1997
    forvalues i = 9/35 {
        destring v`i', gen(Va_`year') force
        //drop v`i'
        loc year = `year' + 1
    }
    
    reshape long Va_, i(statefip StateName industryclassification) j(year)
    ren Va_ Va
    ren industryclassification Naics
}
/********************************************************************************
NAICS-SIC Crosswalks from
From ICPSR, Schaller, Zachary, and DeCelles, Paul
********************************************************************************/
frame create SicToNaics
frame SicToNaics {

    /*
    Description of the task this section accomplishes.
    We must use the provided crosswalk to create two-digit SIC categories which correspond
    to those found in the BEA data from 63to96. The crosswalk will then also give us the
    corresponding two-digit NAICS code.
    */
    use "$IndustryCrosswalks/SIC4_to_NAICS6.dta", clear
    
    gen SIC2 = substr(SIC4,1,2)
    gen SIC3 = substr(SIC4,1,3)
    gen industryclassification = "[01-02]" if inlist(SIC2,"01","02")
    replace industryclassification = "[07-09]" if inlist(SIC2,"07","08","09")
    replace industryclassification = "84, 87, 89" if inlist(SIC2,"84", "87", "89")
    replace industryclassification = "371" if SIC3 == "371"
    replace industryclassification = "372-379" if inlist(SIC3,"372","373","374","375","376","377","378","379")
    replace industryclassification = SIC2 if mi(industryclassification)
    drop SIC3 SIC2 Establishments Est_weight Annual_Payroll Pay_weight Emp_weight
    gen NAICS2 = substr(NAICS6,1,2)
    collapse (sum) Employees, by(industryclassification NAICS2)
    egen EmpTot = total(Employees), by(industryclassification)
    gen Weight = Employees / EmpTot
    drop EmpTot Employees
    reshape wide Weight, i(industryclassification) j(NAICS2) string
    tempfile mergethis
    save `mergethis'
}

/******************************************************************************
Clean SAGDP SIC - The GDP industry table from 1963 to 1997
******************************************************************************/
frame create Gdp63to96 
frame Gdp63to96 {
    import delimited "$Data/gdp_state_industry_BEA/SAGDP_SIC/SAGDP2S__ALL_AREAS_1963_1997.csv", clear
    gen statefip = subinstr(geofips,`"""',"",2)
    replace statefip = substr(statefip,2,2) // there is a leading space for some reason
    drop geofips
    drop if mi(geoname)
    drop if inlist(geoname, "Far West", "Rocky Mountain", "Southwest", "Southeast", ///
    "Plains", "Great Lakes", "Mideast", "New England")
    
    drop if industryclassification == "..."
    drop if inlist(statefip,"00")
    drop unit tablename description linecode region
    loc year = 1963
    forvalues i = 9/43 {
        destring v`i', gen(Va_`year') force
        drop v`i'
        loc year = `year' + 1
    }
    drop Va_1997 // already available in the 97 to 23 frame
    drop if inlist(industryclassification,"A","B","C","D","E","F","G","H","I")
    drop if inlist(industryclassification,"73, 84, 87, 89", "60, 61", "36, 38")
    order statefip
    ren geoname StateName
    reshape long Va_, i(statefip StateName industryclassification) j(year)
    ren Va_ Va

    merge m:1 industryclassification using `mergethis', nogen keep(3)
    sort StateName year industryclassification
    foreach v of varlist Weight* {
        loc Naics2 = substr("`v'",-2,.)
        gen Va_`Naics2' = Va * Weight`Naics2'
        drop Weight`Naics2'
    }

    /*
    NOTES on Dropped Observations:
    On page 94 of the included Census Document by From ICPSR, Schaller, Zachary, and DeCelles, Paul,
    there is no mapping for Sic code 40 (railroad transportation) because "separate establishments primarily engaged in long distance trucking, 
    stevedoring, water transportation, railroad transportation, or pipeline transportation for other establishments of the same enterprise 
    are classified in the corresponding transportation industry."
    */

    /*
    Now reshape and then collapse on the Naics codes and we are done.
    */
    drop Va
    reshape long Va_, i(StateName industryclassification year) j(NAICS2)
    ren Va_ Va
    collapse (sum) Va, by(year StateName statefip NAICS2)
    tostring NAICS2, replace
    gen Naics = ""
    replace Naics = "31-33" if inlist(NAICS2,"31","32","33")
    replace Naics = "44-45" if inlist(NAICS2,"44","45")
    replace Naics = "48-49" if inlist(NAICS2,"48","49")
    replace Naics = NAICS2 if mi(Naics)
    drop NAICS2
    collapse (sum) Va, by(StateName year statefip Naics)
    la var Va ""
    la var StateName ""
}

frame create GdpByIndustry
frame GdpByIndustry {
    xframeappend Gdp63to96 Gdp97to23, drop // You will need to ssc install xframeappend
    drop if year < 1994                    // This is the first year in the rest of my data
    merge m:1 Naics year using `FixedAssetData', keep(3) nogen // 19 observations did not match...but they also didn't have an associated state of fipcode

}

frame copy GdpByIndustry default, replace
frame drop GdpByIndustry
egen TotalVa = total(Va), by(Naics year)
gen VaShare = Va/TotalVa
gen IndCapStock = VaShare * FixedAssets
collapse (sum) K = IndCapStock, by(year statefip StateName)
la var K ""

******** Merge in El-Shagi and Yamarik
frame create fred
frame fred {
    
    import fred A006RD3Q086SBEA, daterange(1994 2023) aggr(a)
    gen year = yofd(daten)
    ren A pk

    * Change base year to 2009
    qui summ pk if year == 2009
    replace pk = pk * 100 / `r(mean)'
    keep year pk
}

frame create ElshagiYamarik 
frame ElshagiYamarik {

    use "${CapStock}/state_capital_yesdata21.dta", clear
    tostring fip, gen(statefip)
    replace statefip = "0" + statefip if strlen(statefip) == 1
    keep statefip year cap

    * Re-inflate (we will deflate later when estimating TFP)
    frlink m:1 year, frame(fred)
    frget pk, from(fred)
    drop fred
    frame drop fred
    replace cap = cap * (pk/100)
    drop pk
    la var cap ""

    tempfile capital
    save `capital', replace 
}

merge 1:1 statefip year using "`capital'", nogen keep(1 3)
ren K K_old
ren cap K

save "$Data/CapitalStockByState.dta", replace







