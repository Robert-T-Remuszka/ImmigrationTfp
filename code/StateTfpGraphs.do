/******************
REQUIRED PACKAGES: You will need to have shp2dta, geo2xy installed to run this code.
> ssc install shp2dta
> ssc install geo2xy
******************/
set mem 10g
clear all
do Globals.do

* Generate dta file from US State Shapefile (2013 Tiger Line file from Census)
shp2dta using "$Data/tl_2013_us_state/tl_2013_us_state.shp", ///
data("$Data/ShpDataBase") coor("$Data/ShpCorr") replace

* Resizing and repositioning Alaska and Hawaii
use "$Data/ShpCorr.dta", clear

* Move Hawaii
drop if _X < -165 & _X != . &  _ID == 32
replace _X = _X  + 55  if  _X != .  &  _ID == 32
replace _Y = _Y  + 4  if _Y != .  &  _ID == 32

* Resize Alaska - and drop some Aleutian Islands
replace _X = _X*.4  -55 if  _X !=  . &  _ID == 41
replace _Y = _Y*.4  + 1 if _Y != .  & _ID == 41
drop if _X > -10 & _X != . & _ID == 41

save "$Data/ShpCorr.dta", replace

* Read in the database file
frame create StateDb
frame StateDb{ 
    use "$Data/ShpDataBase.dta", clear
    ren STUSPS StateAbb
    ren _ID ID                            // For some reason frget won't retrieve this _ in it
}

* Read in Tfp Data
use "$Data/CountyDataTfpEstimates.dta", clear

* Time averages by state
frame copy default TimeAvg
frame TimeAvg {

    sort StateAbb year
    collapse (mean) Z, by(StateAbb)
    
    frlink 1:1 StateAbb, frame(StateDb)
    frget *, from(StateDb)
    ren ID _ID

    * Make map
    spmap Z using "$Data/ShpCorr", id(_ID) ///
    fcolor(Blues)    ///
    legstyle(2) legend(pos(7) size(2.8))   ///
    ocolor(black%30) osize(0.05 ..) clnum(9)
}