/**********************************************************************
Use:
   Out/TfpRegs_PanelBuilt.dta
   Out/TfpRegs_Rotemberg_full.dta
Produces:
   Out/TfpRegs_Table1.xlsx
**********************************************************************/

version 18
clear all
do Globals


use "{Out}/TfpRegs_Rotemberg_full.dta", clear
keep if !missing(weight)

capture confirm variable weight_n
if _rc {
    egen double __Wsum = total(weight)
    egen double __Psum = total(cond(weight>0, weight, .))
    egen double __Nsum = total(cond(weight<0, weight, .))
    replace __Psum = 0 if missing(__Psum)
    replace __Nsum = 0 if missing(__Nsum)

    local Zsum = __Wsum[1]
    if (`Zsum'==. | `Zsum'==0) {
        local Zsum = __Psum[1] + abs(__Nsum[1])
    }
    if (`Zsum'==. | `Zsum'==0) {
        di as err "All weights sum to zero/missing; cannot normalize."
        exit 459
    }

    gen double weight_n = weight/`Zsum'
    gen double abs_w_n  = abs(weight_n)

    drop __Wsum __Psum __Nsum
}

* order by |α| and build presentation columns
gsort -abs_w_n
capture confirm variable rank
if _rc gen int rank = _n
else replace rank = _n

capture confirm variable weight_pct
if _rc gen double weight_pct = 100*weight_n
else   replace    weight_pct = 100*weight_n

capture confirm variable cum_w_pct
if _rc gen double cum_w_pct  = 100*sum(weight_n)
else   replace    cum_w_pct  = 100*sum(weight_n)

tempfile WT
save "`WT'", replace

/* ================= PANEL A ================= */
use "`WT'", clear
quietly summarize weight_n if weight_n<0, meanonly
local A_neg_N    = r(N)
local A_neg_mean = cond(`A_neg_N'>0, r(mean), .)
local A_neg_sum  = cond(`A_neg_N'>0, r(mean)*r(N), 0)

quietly summarize weight_n if weight_n>0, meanonly
local A_pos_N    = r(N)
local A_pos_mean = cond(`A_pos_N'>0, r(mean), .)
local A_pos_sum  = cond(`A_pos_N'>0, r(mean)*r(N), 0)

local A_total     = `A_pos_sum' + `A_neg_sum'
local A_neg_share = cond(`A_total'!=0, `A_neg_sum'/`A_total', .)
local A_pos_share = cond(`A_total'!=0, `A_pos_sum'/`A_total', .)

/* ================= PANEL B inputs: ================= */
tempfile B_GMEANS
preserve
    use "{Out}/TfpRegs_PanelBuilt.dta", clear
    keep if in_samp

    tempname H
    postfile `H' str32 group double gbar varz using "`B_GMEANS'", replace
    foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico WestEu {
        quietly summarize fg_agg`g' if in_samp, meanonly
        local gg = r(mean)
        quietly summarize Bartik_1990_`g' if in_samp, detail
        local vvz = r(Var)
        post `H' ("`g'") (`gg') (`vvz')
    }
    postclose `H'
restore

* β_k and F_k
tempfile B_BETA
preserve
    use "{Out}/TfpRegs_PanelBuilt.dta", clear
    keep if in_samp

    capture confirm variable fg
    if _rc {
        capture drop f
        gen double f = BodiesSupplied_Africa + BodiesSupplied_AsiaOther + BodiesSupplied_CaAuNz + ///
                       BodiesSupplied_China  + BodiesSupplied_EastEu   + BodiesSupplied_India  + ///
                       BodiesSupplied_LA     + BodiesSupplied_Mexico   + BodiesSupplied_Other  + ///
                       BodiesSupplied_WestEu
        bys statefip (year): gen double fg = (f - f[_n-1]) / emp
    }

    tempname HB
    postfile `HB' str32 group double bhat_k Fk using "`B_BETA'", replace
    foreach g in Africa AsiaOther CaAuNz China EastEu India LA Mexico WestEu {
        quietly reghdfe fg Bartik_1990_`g' if in_samp, absorb(statefip year)
        local b  = _b[Bartik_1990_`g']
        test Bartik_1990_`g'
        local Fk = r(F)
        post `HB' ("`g'") (`b') (`Fk')
    }
    postclose `HB'
restore

use "`WT'", clear
keep rank group weight_n
rename weight_n alpha_k
tempfile B_STATS
merge 1:1 group using "`B_GMEANS'", nogen
merge 1:1 group using "`B_BETA'",  nogen
save "`B_STATS'", replace

keep if !missing(alpha_k, gbar, bhat_k, Fk, varz)
count
if r(N) < 2 {
    matrix C = J(5,5,.)
    matrix rownames C = alpha_k gbar bhat_k Fk varz
    matrix colnames C = alpha_k gbar bhat_k Fk varz
}
else {
    capture noisily corr alpha_k gbar bhat_k Fk varz
    if _rc {
        matrix C = J(5,5,.)
        matrix rownames C = alpha_k gbar bhat_k Fk varz
        matrix colnames C = alpha_k gbar bhat_k Fk varz
    }
    else matrix C = r(C)
}

/* ================= PANEL C ================= */
use "`B_STATS'", clear
gen byte pos_g = (gbar>0)
egen double __sum_neg = total(alpha_k*(pos_g==0))
egen double __sum_pos = total(alpha_k*(pos_g==1))
local C_sum_neg = __sum_neg[1]
local C_sum_pos = __sum_pos[1]
drop __sum_neg __sum_pos

/* ================= PANEL D ================= */
use "`WT'", clear
gsort -abs_w_n
capture confirm variable rank
if _rc gen int rank = _n
else replace rank = _n
count
local K = r(N)
if `K' > 10 local K = 10
preserve
    keep in 1/`K'
    tempfile TOPK
    save "`TOPK'", replace
restore

/* ================= PANEL E ================= */
use "`B_STATS'", clear
gen byte sgn = (alpha_k >= 0)
gen double akbk = alpha_k*bhat_k
egen double __Apos   = total(alpha_k*(sgn==1))
egen double __Aneg   = total(alpha_k*(sgn==0))
egen double __sumpos = total(akbk*(sgn==1))
egen double __sumneg = total(akbk*(sgn==0))
local Apos   = __Apos[1]
local Aneg   = __Aneg[1]
local sumpos = __sumpos[1]
local sumneg = __sumneg[1]
local E_pos  = cond(`Apos'!=0 & `Apos'!=., `sumpos'/`Apos', .)
local E_neg  = cond(`Aneg'!=0 & `Aneg'!=., `sumneg'/`Aneg', .)

/* ================= FIRST-STAGE summary ================= */
preserve
    use "{Out}/TfpRegs_PanelBuilt.dta", clear
    keep if in_samp
    capture confirm variable fg
    if _rc {
        capture drop f
        gen double f = BodiesSupplied_Africa + BodiesSupplied_AsiaOther + BodiesSupplied_CaAuNz + ///
                       BodiesSupplied_China  + BodiesSupplied_EastEu   + BodiesSupplied_India  + ///
                       BodiesSupplied_LA     + BodiesSupplied_Mexico   + BodiesSupplied_Other  + ///
                       BodiesSupplied_WestEu
        bys statefip (year): gen double fg = (f - f[_n-1]) / emp
    }
    reghdfe fg q if in_samp, absorb(statefip year)
    test q
    local F_fs = r(F)
    local p_fs = r(p)
restore

/* ================= WRITE EXCEL ================= */
local xlsx = "{Out}/TfpRegs_Table1.xlsx"
cap erase "`xlsx'"
putexcel set "`xlsx'", replace sheet("Table1")

* Panel A
putexcel A1=("Table 1—Summary of Rotemberg Weights: TfpRegs")
putexcel A3=("Panel A. Negative and positive weights")
putexcel A4=("") B4=("Sum") C4=("Mean") D4=("Share")
putexcel A5=("Negative") B5=(`A_neg_sum') C5=(`A_neg_mean') D5=(`A_neg_share')
putexcel A6=("Positive") B6=(`A_pos_sum') C6=(`A_pos_mean') D6=(`A_pos_share')

* Panel B
putexcel A8=("Panel B. Correlations (α_k, g_k, β̂_k, F_k, var(z_k))")
putexcel A9=("") B9=("g_k") C9=("β̂_k") D9=("F_k") E9=("var(z_k)")
putexcel A10=("α_k")
putexcel B10=matrix(C[1,2]) C10=matrix(C[1,3]) D10=matrix(C[1,4]) E10=matrix(C[1,5])
putexcel B11=matrix(C[2,2]) C11=matrix(C[2,3]) D11=matrix(C[2,4]) E11=matrix(C[2,5])
putexcel B12=matrix(C[3,2]) C12=matrix(C[3,3]) D12=matrix(C[3,4]) E12=matrix(C[3,5])
putexcel B13=matrix(C[4,2]) C13=matrix(C[4,3]) D13=matrix(C[4,4]) E13=matrix(C[4,5])
putexcel B14=matrix(C[5,2]) C14=matrix(C[5,3]) D14=matrix(C[5,4]) E14=matrix(C[5,5])

* Panel C
putexcel G3=("Panel C. Variation across years in α_k (by sign of mean g_k)")
putexcel G4=("Sum α_k (g_k<0)") H4=("Sum α_k (g_k>0)")
putexcel G5=(`C_sum_neg') H5=(`C_sum_pos')

* Panel D (TopK)
use "`TOPK'", clear
putexcel A16=("Panel D. Top `K' Rotemberg-weight groups")
putexcel A17=("Rank") B17=("Group") C17=("Wt(%)") D17=("Cum(%)") E17=("Mean s_1990") F17=("SD s_1990") G17=("Mean z") H17=("SD z") I17=("Corr(z,q)") J17=("R2(q~z)") K17=("# States")
local r = 18
quietly {
    sort rank
    forvalues i=1/`=_N' {
        putexcel A`r' = rank[`i']
        putexcel B`r' = group[`i']
        putexcel C`r' = weight_pct[`i']
        putexcel D`r' = cum_w_pct[`i']
        capture noisily putexcel E`r' = mean_s[`i']
        capture noisily putexcel F`r' = sd_s[`i']
        capture noisily putexcel G`r' = mean_z[`i']
        capture noisily putexcel H`r' = sd_z[`i']
        capture noisily putexcel I`r' = corr_z_q[`i']
        capture noisily putexcel J`r' = r2_q_on_z[`i']
        capture noisily putexcel K`r' = n_states_exposed[`i']
        local ++r
    }
}

* Panel E
putexcel G9=("Panel E. α-weighted β̂ by sign of α")
putexcel G10=("Negative α:") H10=(`E_neg')
putexcel G11=("Positive α:") H11=(`E_pos')

* First-stage summary
putexcel A35=("First-stage (fg on q, FE) — F-stat:") B35=(`F_fs') C35=("p-value:") D35=(`p_fs')
putexcel save

putexcel set "`xlsx'", sheet("TopK") modify
use "`TOPK'", clear
putexcel A1=("Rank") B1=("Group") C1=("Wt(%)") D1=("Cum(%)") E1=("Mean s_1990") F1=("SD s_1990") G1=("Mean z") H1=("SD z") I1=("Corr(z,q)") J1=("R2(q~z)") K1=("# States")
local r = 2
quietly {
    sort rank
    forvalues i=1/`=_N' {
        putexcel A`r' = rank[`i']
        putexcel B`r' = group[`i']
        putexcel C`r' = weight_pct[`i']
        putexcel D`r' = cum_w_pct[`i']
        capture noisily putexcel E`r' = mean_s[`i']
        capture noisily putexcel F`r' = sd_s[`i']
        capture noisily putexcel G`r' = mean_z[`i']
        capture noisily putexcel H`r' = sd_z[`i']
        capture noisily putexcel I`r' = corr_z_q[`i']
        capture noisily putexcel J`r' = r2_q_on_z[`i']
        capture noisily putexcel K`r' = n_states_exposed[`i']
        local ++r
    }
}
putexcel save



