clear all
do Globals
do Functions

use "${Data}/StateAnalysis.dta", clear

* Define the sample
loc samp inrange(year, 1994, 2021)

* Construct the Bartik instruments and left hand side variables - See Functions.do
qui PreRegProcessing
sort state year

loc dkraayband 9
loc depvarlags 4
loc ivlags     4
loc binnum     30

* Generate macro containing first-differenced variables
loc vars "Z Wage_Domestic Wage_Foreign L CapStock"
foreach v in `vars' {
    
    gen D0`v' = ln(`v' / L.`v')

}

/********
Visual IV - Full Bartik
*********/

* Residualize
qui reg fg L(1/`depvarlags').(fg) L(1/`ivlags').Bartik_1990 i.year i.state [pw = emp] if `samp'
predict double e_fg, residuals

qui reg Bartik_1990 L(1/`depvarlags').(fg) L(1/`ivlags').Bartik_1990 i.year i.state [pw = emp] if `samp'
predict double e_B, residuals

/******
Binscatter
******/
* Generate evenly size bins according to Bartik IV
xtile bin = e_B if `samp', nq(`binnum')

* Plot a binscatter
preserve
    collapse (mean) e_fg e_B [aw = emp], by(bin)
    tw (scatter e_fg e_B if bin != `binnum', mc(ebblue) ms(oh)) ///
       (scatter e_fg e_B if bin == `binnum', mc(ebblue)) ///
       (lfit    e_fg e_B, lc(orange) lp(dash)), ///
    xlab(, nogrid) ylab(, nogrid) legend(off) ytitle("{&Delta}{sup:0}ln(L{sup:F}), Residualized") ///
    xtitle("Bartik IV, 1990 Shares Residualized") name(Bartik1990_bs)
restore

* Look at which years are in the extreme bins
qui summ year if bin == `binnum'
hist year if bin == `binnum', d xlab(`r(min)'(1)`r(max)', nogrid) ylab(, nogrid) name(Bartik1990_olr_yr) xtitle("")

* Look at which states are in the extreme bins
encode statename, gen(stnum)
qui summ stnum if bin == `binnum'
hist stnum if bin == `binnum', d xlab(1(1)`r(max)', labsize(tiny) angle(45) valuelabel nogrid) ///
ylab(, nogrid) xtitle("") name(Bartik1990_olr_state)

* Export the binscatter
graph export "${Graphs}/Bartik1990_bs.pdf", replace name(Bartik1990_bs)

