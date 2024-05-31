clear
set more off

insheet using "base_perturbed.csv"


gen lnsigma = log( 1 / (1 - rho) )
gen ndf = 21

gen obsnum = _n

rename id labnumber

merge m:m labnumber using "base_estimated_M2_1.dta", keep(match) keepusing(alphai deltai se_alpha se_delta ratei lnsigmai cparami se_rate se_lnsigma se_cparam cov_cparamdelta cil_delta ciu_delta neqone_delta) nogenerate

rename lag k
gen t0 = 0
replace t0 = 1 if time == 0
rename price_sooner pratio
rename cons_sooner cons1
rename income mnumer

replace cons1 = cons1 * pratio / mnumer

capture drop convergei
gen convergei = .

capture drop betai
gen betai = .

capture drop se_beta
gen se_beta = .

capture drop alphai_t0
capture drop deltai_t0
gen alphai_t0 = .
gen deltai_t0 = .

capture drop se_alpha_t0
capture drop se_delta_t0
gen se_alpha_t0 = .
gen se_delta_t0 = .

capture drop ratei_t0
capture drop lnsigmai_t0
capture drop cparami_t0
gen ratei_t0 = .
gen lnsigmai_t0 = .
gen cparami_t0 = .

capture drop se_rate_t0
capture drop se_lnsigma_t0
capture drop se_cparam_t0
gen se_rate_t0 = .
gen se_lnsigma_t0 = .
gen se_cparam_t0 = .

capture drop cov_cparamdelta_t0
gen cov_cparamdelta_t0 = .

capture drop ic
capture drop mss
capture drop mms
capture drop msr
capture drop rmse
capture drop ll
capture drop r2
capture drop r2_a
capture drop dev
capture drop rss
capture drop tss
gen ic = .
gen mss = .
gen mms = .
gen msr = .
gen rmse = .
gen ll = .
gen r2 = .
gen r2_a = .
gen dev = .
gen rss = .
gen tss = .

capture drop cil_beta
capture drop ciu_beta
capture drop neqone_beta
gen cil_beta = .
gen ciu_beta = .
gen neqone_beta = 0




levelsof labnumber, local(levels)

foreach i of local levels {

    disp "==========================="
    disp "subj:  " `i'
    disp "==========================="

    estimates clear
    return clear

    replace convergei = -100 if labnumber == `i'

    local TMPcparam = -0.07
    local TMPdelta = 1.0

    local NEWalpha = " ( 1 - exp( - 4 * tanh( {cparam=`TMPcparam'} ) - 1.5 ) ) "

    local DENOM1      = " pratio                                      ^ ( `NEWalpha' / ( `NEWalpha' - 1 ) ) "
    local DENOM2      = " (                    {delta=`TMPdelta'}^k ) ^ ( 1.0        / ( `NEWalpha' - 1 ) ) "

    local NUMER      = " (                    {delta=`TMPdelta'}^k * pratio ) ^ ( 1.0 / ( `NEWalpha' - 1 ) ) "

    capture nl ( cons1 = /*
*/                   ( /*
*/                       ( `NUMER'      * mnumer ) / ( 1 + `DENOM1' * `DENOM2'      ) /*
*/                   ) / ( mnumer / pratio ) /*
*/             ) if labnumber == `i' & t0 == 1, vce(jack) iter(200) eps(1e-5)


    if _rc == 0 {

        nl

        replace convergei = e(converge) if labnumber == `i'

        replace ic   = e(ic  ) if labnumber == `i'
        replace mss  = e(mss ) if labnumber == `i'
        replace mms  = e(mms ) if labnumber == `i'
        replace msr  = e(msr ) if labnumber == `i'
        replace rmse = e(rmse) if labnumber == `i'
        replace ll   = e(ll  ) if labnumber == `i'
        replace r2   = e(r2  ) if labnumber == `i'
        replace r2_a = e(r2_a) if labnumber == `i'
        replace dev  = e(dev ) if labnumber == `i'
        replace rss  = e(rss ) if labnumber == `i'
        replace tss  = e(tss ) if labnumber == `i'


        capture matrix covmat =  e(V)
        if _rc == 0 {
            matrix list covmat

            replace cov_cparamdelta_t0 = covmat[1,1] if labnumber == `i'
        }


        capture estimates
        if _rc == 0 {
            matrix rtable = r(table)
            matrix list rtable

            replace deltai_t0 = rtable[1,1] if labnumber == `i'
            replace cparami_t0 = rtable[1,2] if labnumber == `i'

            replace se_delta_t0 = rtable[2,1] if labnumber == `i'
            replace se_cparam_t0 = rtable[2,2] if labnumber == `i'
        }


        replace lnsigmai_t0 = 4 * tanh( cparami_t0 ) + 1.5 if labnumber == `i'
        replace se_lnsigma_t0 = 4 * ( 1 - tanh( cparami_t0 ) ^ 2 ) * abs( se_cparam_t0 ) if labnumber == `i'

        replace alphai_t0 = 1 - exp( - 4 * tanh( cparami_t0 ) - 1.5 ) if labnumber == `i'
        replace se_alpha_t0 = 4 * ( 1 - tanh( cparami_t0 ) ^ 2 ) * exp( - 4 * tanh( cparami_t0 ) - 1.5 ) * abs( se_cparam_t0 ) if labnumber == `i'

        replace ratei_t0 = max(1e-5, deltai_t0 ^ (-365)) - 1 if labnumber == `i'
        replace se_rate_t0 = abs( (-365) * max(1e-5, deltai_t0 ^ (-366)) ) * abs( se_delta_t0 ) if labnumber == `i'

        replace betai = ( deltai_t0 / deltai ) ^ 70
        replace se_beta = (( 70 * (deltai_t0 / deltai)^70 / deltai_t0 )^2 * ( se_delta_t0 )^2 + ( -70 * (deltai_t0 / deltai)^70 / deltai )^2 * ( se_delta )^2 )^(1/2)

    }

    drop if labnumber == `i' & mod(obsnum, 42) != 0

    replace cil_beta    = betai    - se_beta    * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace ciu_beta    = betai    + se_beta    * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace neqone_beta  = 1 if labnumber == `i' & beta  < 1 & ciu_beta  < 1 & se_beta != .
    replace neqone_beta  = 1 if labnumber == `i' & beta  > 1 & cil_beta  > 1 & se_beta != .
}

drop time k pratio price_later mnumer cons1 cons_later noise obsnum t0


save "base_estimated_M4.dta", replace
