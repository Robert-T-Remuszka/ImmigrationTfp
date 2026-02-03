/**********************************************************************
 Generates:
   1) Out/TfpRegs_PanelBuilt.dta
   2) Out/TfpRegs_Rotemberg_full.dta
   3) Out/TfpRegs_Rotemberg_full.csv
**********************************************************************/

version 18
clear all
set more off
do Globals
cap mkdir "{Out}"

*-----------------------------------------------------------
preserve
    use "{Data}/StateAnalysisFile.dta", clear
    keep statefip year ImmigrantGroup BodiesSupplied
    capture confirm string variable statefip
    if _rc tostring statefip, replace
    replace statefip = "0" + statefip if strlen(statefip)==1
    replace ImmigrantGroup = strtrim(ImmigrantGroup)

    tempfile _map
    clear
    input str40 ImmigrantGroup str12 abbr
    "Canada-Australia-New Zealand" "_CaAuNz"
    "United States"                  "_US"
    "Western Europe"                 "_WestEu"
    "Latin America"                  "_LA"
    "Russia and Eastern Europe"      "_EastEu"
    "Rest of Asia"                   "_AsiaOther"
    "China"                          "_China"
    "Africa"                         "_Africa"
    "Mexico"                         "_Mexico"
    "Other"                          "_Other"
    "India"                          "_India"
    end
    save "`_map'", replace

    use "{Data}/StateAnalysisFile.dta", clear
    keep statefip year ImmigrantGroup BodiesSupplied
    capture confirm string variable statefip
    if _rc tostring statefip, replace
    replace statefip = "0" + statefip if strlen(statefip)==1
    replace ImmigrantGroup = strtrim(ImmigrantGroup)

    merge m:1 ImmigrantGroup using "`_map'", keep(match) nogen
    collapse (sum) BodiesSupplied, by(statefip year abbr)
    reshape wide BodiesSupplied, i(statefip year) j(abbr) string

    quietly ds BodiesSupplied_*
    foreach v of varlist `r(varlist)' {
        replace `v' = 0 if missing(`v')
    }

    tempfile WIDE
    save "`WIDE'", replace
restore


capture noisily use "{Data}/StateAnalysisFileTfp.dta", clear
if _rc {
    capture noisily import delimited "{Data}/StateAnalysisFileTfp.csv", clear case(preserve)
    if _rc {
        di as err "Cannot load Data/StateAnalysisFileTfp.{dta,csv}"
        exit 498
    }
}

capture confirm string variable statefip
if _rc tostring statefip, replace
replace statefip = "0" + statefip if strlen(statefip)==1

bysort statefip year: gen byte __row=_n
keep if __row==1
drop __row

merge m:1 statefip year using "`WIDE'", keep(match) nogen

gen byte in_samp = inrange(year, 1994, 2021)

egen double emp = rowtotal(BodiesSupplied_*), missing

* exclude US
local GLIST "Africa AsiaOther CaAuNz China EastEu India LA Mexico Other WestEu US"
foreach g of local GLIST {
    if "`g'"=="US" continue
    capture drop s_`g'
    gen double s_`g' = BodiesSupplied_`g'/emp
}

capture unab v90 : *1990
if _rc==0 & "`v90'"!="" {
    egen double emp1990 = rowtotal(`v90'), missing
    foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico Other WestEu US {
        capture confirm variable `g'1990
        if !_rc & "`g'"!="US" {
            capture drop s_`g'_1990
            gen double s_`g'_1990 = `g'1990/emp1990
        }
    }
    foreach v of varlist s_*_1990 {
        local base = subinstr("`v'","_1990","",1)
        bys statefip (year): egen double __tmp = mean(cond(year==1990, `base', .))
        replace `v' = __tmp if missing(`v')
        drop __tmp
    }
}
else {
    foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico Other WestEu {
        capture drop s_`g'_1990
        bys statefip (year): egen double s_`g'_1990 = mean(cond(year==1990, s_`g', .))
        replace s_`g'_1990 = 0 if missing(s_`g'_1990)
    }
}

* national growth by group: fg_agg<g>
preserve
    keep year BodiesSupplied_*
    collapse (sum) BodiesSupplied_*, by(year)
    egen double emp_nat = rowtotal(BodiesSupplied_*), missing
    tsset year
    foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico Other WestEu US {
        gen double fg_agg`g' = (BodiesSupplied_`g' - L.BodiesSupplied_`g')/L.emp_nat
    }
    keep year fg_agg*
    tempfile NAT
    save "`NAT'", replace
restore
merge m:1 year using "`NAT'", nogen

* Bartik_1990
foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico WestEu {
    capture drop Bartik_1990_`g'
    gen double Bartik_1990_`g' = s_`g'_1990 * fg_agg`g'
}
egen double q = rowtotal(Bartik_1990_*), missing

compress
save "{Out}/TfpRegs_PanelBuilt.dta", replace

/*-----------------------------------------------------------
  Rotemberg weights
-----------------------------------------------------------*/

* foreign totals (exclude US)
gen double f = BodiesSupplied_Africa + BodiesSupplied_AsiaOther + BodiesSupplied_CaAuNz + ///
               BodiesSupplied_China  + BodiesSupplied_EastEu   + BodiesSupplied_India  + ///
               BodiesSupplied_LA     + BodiesSupplied_Mexico   + BodiesSupplied_Other  + ///
               BodiesSupplied_WestEu

bys statefip (year): gen double fg = (f - f[_n-1]) / emp

reghdfe fg if in_samp, absorb(statefip year) resid
predict double Xtilde, resid
reghdfe q  if in_samp, absorb(statefip year) resid
predict double qtilde, resid

reghdfe fg q if in_samp, absorb(statefip year)
predict double Xhat_B, xb
test q
scalar F_fs = r(F)
scalar p_fs = r(p)

summ Xhat_B if in_samp, meanonly
scalar mB = r(mean)
gen double Xhat_B_ct = Xhat_B - mB
summ Xtilde if in_samp, meanonly
scalar mX = r(mean)
gen double Xtilde_ct = Xtilde - mX

gen double __prodB = Xhat_B_ct*Xtilde_ct if in_samp
quietly summ __prodB
scalar denom = r(mean)
drop __prodB

* post Rotemberg rows
tempfile ROT
tempname H_ROT
postfile `H_ROT' str32 group double weight abs_w ///
    double mean_s sd_s mean_z sd_z corr_z_q r2_q_on_z ///
    double n_states_exposed using "`ROT'", replace

preserve
    keep statefip s_*_1990
    collapse (max) s_*_1990, by(statefip)
    tempfile EVER
    save "`EVER'", replace
restore

foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico WestEu {
    local z Bartik_1990_`g'

    * just-identified first stage
    reghdfe fg `z' if in_samp, absorb(statefip year)
    predict double Xhat_k, xb
    summ Xhat_k if in_samp, meanonly
    scalar mk = r(mean)
    gen double Xhat_k_ct = Xhat_k - mk
    gen double __prodk = Xhat_k_ct * Xtilde_ct if in_samp
    quietly summ __prodk
    scalar num = r(mean)
    drop __prodk Xhat_k_ct Xhat_k

    scalar wk = num/denom

    quietly summ s_`g'_1990 if in_samp, meanonly
    scalar ms = r(mean)
    scalar ss = r(sd)

    quietly summ `z' if in_samp, meanonly
    scalar mz = r(mean)
    scalar sz = r(sd)

    corr `z' q if in_samp
    scalar czq = r(rho)

    regress q `z' if in_samp
    scalar r2 = e(r2)

    preserve
        use "`EVER'", clear
        gen byte exposed = (s_`g'_1990>0)
        collapse (max) exposed, by(statefip)
        count if exposed==1
        scalar nexp = r(N)
    restore

    post `H_ROT' ("`g'") (wk) (abs(wk)) (ms) (ss) (mz) (sz) (czq) (r2) (nexp)
}
postclose `H_ROT'

use "`ROT'", clear
gsort -abs_w
gen int    rank       = _n
gen double weight_pct = 100*weight
gen double cum_w_pct  = 100*sum(weight)

order rank group weight weight_pct cum_w_pct mean_s sd_s mean_z sd_z corr_z_q r2_q_on_z n_states_exposed abs_w
compress
save "{Out}/TfpRegs_Rotemberg_full.dta", replace
export delimited using "{Out}/TfpRegs_Rotemberg_full.csv", replace


