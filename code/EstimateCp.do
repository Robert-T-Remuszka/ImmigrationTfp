clear all
do Globals

/*================================================================
Estimates the "Choice Probabilities" equation -- as an
individual-level probit on MakeIndividualAnalysis.do's ACS 2000-2022 sample:

    Pr(Foreign | tau, w) = Phi( b0 + b_w*w_tilde + b_tau*Phi^-1(tau) )

with  b_w = 1/xi_w,  b_tau = -xi_z/xi_w,  b0 = -mu_z/xi_w . This
identifies (mu_z, xi_w, xi_z) directly off the probit's three coefficients

w_tilde varies only at the state-year level (it's a market-clearing factor
price in the model, not an individual object) while tau varies only at the
occupation level -- both merged onto individual rows in
MakeIndividualAnalysis.do. Clustering by STATEFIP accounts for w_tilde's
coarser level of variation
================================================================*/

use "${Data}/IndividualCpAnalysis.dta", clear

count
tab Foreign
summ w_tilde phi_inv_tau [aw = PERWT]

/******************* PROBIT **************************************/
probit Foreign w_tilde phi_inv_tau [pweight = PERWT], cluster(STATEFIP)

loc N = e(N)

/******************* BACK OUT (mu_z, xi_w, xi_z) VIA THE DELTA METHOD ************/
nlcom (xi_w: 1/_b[w_tilde]) ///
      (xi_z: -_b[phi_inv_tau]/_b[w_tilde]) ///
      (mu_z: -_b[_cons]/_b[w_tilde]), post

loc xi_w  = _b[xi_w]
loc se_xi_w = _se[xi_w]
loc xi_z  = _b[xi_z]
loc se_xi_z = _se[xi_z]
loc mu_z  = _b[mu_z]
loc se_mu_z = _se[mu_z]

di as text "N (probit)      = " as result `N'
di as text "xi_w-hat (se)   = " as result `xi_w' as text "  (" as result `se_xi_w' as text ")"
di as text "xi_z-hat (se)   = " as result `xi_z' as text "  (" as result `se_xi_z' as text ")"
di as text "mu_z-hat (se)   = " as result `mu_z' as text "  (" as result `se_mu_z' as text ")"

/******************* SAVE FOR JULIA CONSUMPTION ***********************************/
clear
set obs 1
gen N          = `N'
gen xi_w_hat   = `xi_w'
gen xi_w_se    = `se_xi_w'
gen xi_z_hat   = `xi_z'
gen xi_z_se    = `se_xi_z'
gen mu_z_hat   = `mu_z'
gen mu_z_se    = `se_mu_z'

la var xi_w_hat "sigma_omega: dispersion of producer heterogeneity (tilde-omega ~ N(0,xi_w^2))"
la var xi_z_hat "sigma_z: dispersion of comparative advantage z(tau) across tasks"
la var mu_z_hat "location of z(tau)"

save "${Data}/CpEstimates.dta", replace
di as text "Wrote ${Data}/CpEstimates.dta"
