clear
set more off

insheet using "base_perturbed.csv"


gen lnsigma = log( 1 / (1 - rho) )
gen ndf = 42

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

    local TMPbeta  = 1.0

    local DENOM1      = " pratio                                      ^ ( alphai / ( alphai - 1 ) ) "
    local DENOM2prsnt = " ( {beta=`TMPbeta'} * deltai ^k ) ^ ( 1.0        / ( alphai - 1 ) ) "
    local DENOM2      = " (                    deltai ^k ) ^ ( 1.0        / ( alphai - 1 ) ) "

    local NUMERprsnt = " ( {beta=`TMPbeta'} * deltai ^k * pratio ) ^ ( 1.0 / ( alphai - 1 ) ) "
    local NUMER      = " (                    deltai ^k * pratio ) ^ ( 1.0 / ( alphai - 1 ) ) "

    capture nl ( cons1 = /*
*/                   ( /*
*/                       t0         * (  ( `NUMERprsnt' * mnumer ) / ( 1 + `DENOM1' * `DENOM2prsnt' )  ) /*
*/                     + ( 1 - t0 ) * (  ( `NUMER'      * mnumer ) / ( 1 + `DENOM1' * `DENOM2'      )  ) /*
*/                   ) / ( mnumer / pratio ) /*
*/             ) if labnumber == `i', vce(jack) iter(200) eps(1e-5)


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


        capture estimates
        if _rc == 0 {
            matrix rtable = r(table)
            matrix list rtable

            replace betai  = rtable[1,1] if labnumber == `i'

            replace se_beta  = rtable[2,1] if labnumber == `i'
        }

    }

    drop if labnumber == `i' & mod(obsnum, ndf) != 0

    replace cil_beta    = betai    - se_beta    * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace ciu_beta    = betai    + se_beta    * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace neqone_beta  = 1 if labnumber == `i' & beta  < 1 & ciu_beta  < 1 & se_beta != .
    replace neqone_beta  = 1 if labnumber == `i' & beta  > 1 & cil_beta  > 1 & se_beta != .
}

drop time k pratio price_later mnumer cons1 cons_later noise obsnum t0


save "base_estimated_M2_2.dta", replace
