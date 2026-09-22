clear all
do Globals

/*================================================================
Builds migrant-enclave weighted counts by state, from the 1960-2000
decennial census extracts in data/census-enclaves/ (Census1960.dta -
Census2000.dta): PERWT-weighted counts by state x birthplace (BPL) and,
where fielded (1980 onward -- not asked in 1960/1970), by state x
self-reported ancestry (ANCESTR1).

Two counts are built for each state x origin-group x decade cell: total
population, and a full-time-worker subset (UHRSWORK>=35, also only
available 1980 onward), so the eventual shift-share instrument can be
normalized by either population or labor-force size without needing to
re-pull anything.

ANCESTR1 codes 995 (Mixture), 996 (Uncodable), 998 (Other), and 999 (Not
Reported) are dropped as non-substantive -- confirmed directly against
IPUMS's own codes/frequencies table. Code 994 (North American) and
everything below it, including state-identity codes like 989 (Wisconsin)
and 993 (Southerner), are real, substantive responses and are kept.

Four long panels (state x origin x decade, for {BPL,ancestry} x
{population,full-time}) are each reshaped wide on decade and merged into
two output files, one row per state x origin-group:
  data/EnclaveSharesBpl.dta:      count_bpl_pop_1960-count_bpl_pop_2000,
                                   count_bpl_ft_1980-count_bpl_ft_2000
  data/EnclaveSharesAncestry.dta: count_ancestr1_pop_1980-count_ancestr1_pop_2000,
                                   count_ancestr1_ft_1980-count_ancestr1_ft_2000
================================================================*/

frame create BplPopPanel
frame create BplFtPanel
frame create AncestryPopPanel
frame create AncestryFtPanel

loc files: dir "${Data}/census-enclaves/" files "Census*.dta"
foreach f in `files' {

    use "${Data}/census-enclaves/`f'", clear

    tostring STATEFIP, replace
    replace STATEFIP = "0" + STATEFIP if strlen(STATEFIP) == 1

    drop if BPL >= 900

    capture confirm variable UHRSWORK
    loc hasUhrs = !_rc

    preserve
        collapse (sum) PERWT, by(STATEFIP YEAR BPL)
        ren PERWT Count
        ren YEAR Year
        frame BplPopPanel: xframeappend default
    restore

    if `hasUhrs' {
        preserve
            keep if UHRSWORK >= 35
            collapse (sum) PERWT, by(STATEFIP YEAR BPL)
            ren PERWT Count
            ren YEAR Year
            frame BplFtPanel: xframeappend default
        restore
    }

    capture confirm variable ANCESTR1
    if !_rc {
        preserve
            drop if inlist(ANCESTR1, 995, 996, 998, 999)
            collapse (sum) PERWT, by(STATEFIP YEAR ANCESTR1)
            ren PERWT Count
            ren YEAR Year
            frame AncestryPopPanel: xframeappend default
        restore

        if `hasUhrs' {
            preserve
                drop if inlist(ANCESTR1, 995, 996, 998, 999)
                keep if UHRSWORK >= 35
                collapse (sum) PERWT, by(STATEFIP YEAR ANCESTR1)
                ren PERWT Count
                ren YEAR Year
                frame AncestryFtPanel: xframeappend default
            restore
        }
    }

}

/*================================================================
  BPL OUTPUT
================================================================*/
frame change BplFtPanel
reshape wide Count, i(STATEFIP BPL) j(Year)
ds Count*
foreach v of varlist `r(varlist)' {
    loc yr = substr("`v'", 6, .)
    ren `v' count_bpl_ft_`yr'
    replace count_bpl_ft_`yr' = 0 if mi(count_bpl_ft_`yr')
}
tempfile BplFtFile
save `BplFtFile'

frame change BplPopPanel
reshape wide Count, i(STATEFIP BPL) j(Year)
ds Count*
foreach v of varlist `r(varlist)' {
    loc yr = substr("`v'", 6, .)
    ren `v' count_bpl_pop_`yr'
    replace count_bpl_pop_`yr' = 0 if mi(count_bpl_pop_`yr')
}
merge 1:1 STATEFIP BPL using `BplFtFile', nogen
la var STATEFIP "State FIPS code"
la var BPL       "IPUMS birthplace code"
sort STATEFIP BPL
save "${Data}/EnclaveSharesBpl.dta", replace
di as text "Wrote ${Data}/EnclaveSharesBpl.dta"

/*================================================================
  ANCESTRY OUTPUT
================================================================*/
frame change AncestryFtPanel
reshape wide Count, i(STATEFIP ANCESTR1) j(Year)
ds Count*
foreach v of varlist `r(varlist)' {
    loc yr = substr("`v'", 6, .)
    ren `v' count_ancestr1_ft_`yr'
    replace count_ancestr1_ft_`yr' = 0 if mi(count_ancestr1_ft_`yr')
}
tempfile AncestryFtFile
save `AncestryFtFile'

frame change AncestryPopPanel
reshape wide Count, i(STATEFIP ANCESTR1) j(Year)
ds Count*
foreach v of varlist `r(varlist)' {
    loc yr = substr("`v'", 6, .)
    ren `v' count_ancestr1_pop_`yr'
    replace count_ancestr1_pop_`yr' = 0 if mi(count_ancestr1_pop_`yr')
}
merge 1:1 STATEFIP ANCESTR1 using `AncestryFtFile', nogen
la var STATEFIP "State FIPS code"
la var ANCESTR1  "IPUMS first-reported ancestry code"
sort STATEFIP ANCESTR1
save "${Data}/EnclaveSharesAncestry.dta", replace
di as text "Wrote ${Data}/EnclaveSharesAncestry.dta"
