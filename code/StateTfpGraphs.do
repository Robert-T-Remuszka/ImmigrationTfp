/******************
REQUIRED PACKAGES: You will need to have shp2dta, geo2xy installed to run this code.
> ssc install shp2dta
> ssc install geo2xy
******************/
set mem 10g
clear all
do Globals.do
loc graphs "TfpEstimates2019 TfpGrowth94to23"

* Generate dta file from US State Shapefile (2013 Tiger Line file from Census)
shp2dta using "$Data/cb_2023_us_state_500k/cb_2023_us_state_500k.shp", ///
data("$Data/ShpDataBase") coor("$Data/ShpCorr") replace

* Resizing and repositioning Alaska and Hawaii
use "$Data/ShpCorr.dta", clear

* Move Hawaii
drop if _X < -165 & _X != . &  _ID == 19
replace _X = _X  + 55  if  _X != .  &  _ID == 39
replace _Y = _Y  + 4  if _Y != .  &  _ID == 39

* Resize Alaska - and drop some Aleutian Islands
replace _X = _X*.4  -55 if  _X !=  . &  _ID == 19
replace _Y = _Y*.4  + 1 if _Y != .  & _ID == 19
drop if _X > -10 & _X != . & _ID == 19

save "$Data/ShpCorr.dta", replace

* Read in the database file
frame create StateDb
frame StateDb { 
    use "$Data/ShpDataBase.dta", clear
    ren NAME StateName
    ren _ID ID                            // For some reason frget won't retrieve this _ in it
}

* Read in Tfp Data
import delimited "$Data/StateAnalysisFileTfp.csv", clear case(preserve)

* Levels by state
frame copy default Collapse
frame Collapse {

    sort StateName year
    collapse Z if year == 2019, by(StateName)
    
    frlink 1:1 StateName, frame(StateDb)
    frget *, from(StateDb)
    ren ID _ID
    
    egen z = std(Z)
    
    * Make map
    format z %12.2f

    spmap z using "$Data/ShpCorr", id(_ID) ///
    fcolor(Blues)    ///
    legstyle(2) legend(pos(7) size(2.8))   ///
    ocolor(black%30) osize(0.05 ..) cln(9) legend(pos(4) title("Standardized Values", size(small))) ///
    name(TfpEstimates2019) 
}


* Growth rates by state
frame copy default Collapse, replace
frame Collapse {
    
    sort StateName year
    bysort StateName (year): gen gz = log(Z[_N]/Z[1]) * 100
    collapse (firstnm) gz, by(StateName)

    frlink 1:1 StateName, frame(StateDb)
    frget *, from(StateDb)
    ren ID _ID

    * Make map
    format gz %12.1f
    spmap gz using "$Data/ShpCorr", id(_ID) ///
    fcolor(Blues) legstyle(2) legend(pos(7) size(2.8)) ocolor(black%30) osize(0.05 ..) cln(9) ///
    legen(pos(4)) name(TfpGrowth94to23)
}

* Export graphs
foreach g in `graphs' {
    graph export "$Graphs/`g'.pdf", replace name(`g')
}