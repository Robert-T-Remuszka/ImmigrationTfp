clear all
do Globals

use "${Data}/StateAnalysisPreTfp.dta", clear

frame create mergethis
frame mergethis {

    import delimited "${Data}/StateTfpAndTaskAgg.csv", clear
    tostring statefip, replace
    replace statefip = "0" + statefip if strlen(statefip) == 1

    ren z Z
    ren l L
    la var Z "Labor-augmenting productivity"
    la var L "Task aggregate"

    tempfile mergethis
    save `mergethis', replace
}

merge 1:1 statefip year using `mergethis', nogen

save "${Data}/StateAnalysis.dta", replace