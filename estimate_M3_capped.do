clear
set more off

insheet using "base_perturbed.csv"


gen lnsigma = log( 1 / (1 - rho) )
gen ndf = 42

gen obsnum = _n

rename id labnumber
rename lag k
gen t0 = 0
replace t0 = 1 if time == 0
rename price_sooner pratio
rename cons_sooner cons1
rename income mnumer

replace cons1 = cons1 * pratio / mnumer

capture drop convergei
gen convergei = .

capture drop alphai
capture drop betai
capture drop deltai
gen alphai = .
gen betai = .
gen deltai = .

capture drop se_alpha
capture drop se_beta
capture drop se_delta
gen se_alpha = .
gen se_beta = .
gen se_delta = .

capture drop ratei
capture drop lnsigmai
capture drop cparami
gen ratei = .
gen lnsigmai = .
gen cparami = .

capture drop se_rate
capture drop se_lnsigma
capture drop se_cparam
gen se_rate = .
gen se_lnsigma = .
gen se_cparam = .

capture drop cov_betadelta
capture drop cov_cparambeta
capture drop cov_cparamdelta
gen cov_betadelta = .
gen cov_cparambeta = .
gen cov_cparamdelta = .

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

capture drop cil_delta
capture drop ciu_delta
capture drop cil_beta
capture drop ciu_beta
capture drop neqone_delta
capture drop neqone_beta
gen cil_delta = .
gen ciu_delta = .
gen cil_beta = .
gen ciu_beta = .
gen neqone_delta = 0
gen neqone_beta = 0


local lowerrho = 1 - exp( -2.5 )
gen capped_rho = rho
replace capped_rho = `lowerrho' if capped_rho > `lowerrho'

drop if capped_rho == rho


levelsof labnumber, local(levels)

foreach i of local levels {

    disp "==========================="
    disp "subj:  " `i'
    disp "==========================="

    estimates clear
    return clear

    replace convergei = -100 if labnumber == `i'

    local TMPbeta  = 1.0
    local TMPdelta = 1.0

    local DENOM1      = " pratio                                      ^ ( capped_rho / ( capped_rho - 1 ) ) "
    local DENOM2prsnt = " ( {beta=`TMPbeta'} * {delta=`TMPdelta'}^k ) ^ ( 1.0        / ( capped_rho - 1 ) ) "
    local DENOM2      = " (                    {delta=`TMPdelta'}^k ) ^ ( 1.0        / ( capped_rho - 1 ) ) "

    local NUMERprsnt = " ( {beta=`TMPbeta'} * {delta=`TMPdelta'}^k * pratio ) ^ ( 1.0 / ( capped_rho - 1 ) ) "
    local NUMER      = " (                    {delta=`TMPdelta'}^k * pratio ) ^ ( 1.0 / ( capped_rho - 1 ) ) "

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


        capture matrix covmat =  e(V)
        if _rc == 0 {
            matrix list covmat

            replace cov_betadelta   = covmat[2,1] if labnumber == `i'
        }


        capture estimates
        if _rc == 0 {
            matrix rtable = r(table)
            matrix list rtable

            replace betai  = rtable[1,1] if labnumber == `i'
            replace deltai = rtable[1,2] if labnumber == `i'

            replace se_beta  = rtable[2,1] if labnumber == `i'
            replace se_delta = rtable[2,2] if labnumber == `i'
        }


        replace ratei = max(1e-5, deltai ^ (-365)) - 1 if labnumber == `i'
        replace se_rate = abs( (-365) * max(1e-5, deltai ^ (-366)) ) * abs( se_delta ) if labnumber == `i'

    }

    drop if labnumber == `i' & mod(obsnum, ndf) != 0

    replace cil_delta   = deltai   - se_delta   * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace ciu_delta   = deltai   + se_delta   * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace cil_beta    = betai    - se_beta    * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace ciu_beta    = betai    + se_beta    * invttail(ndf-1, 0.05/2) if labnumber == `i'
    replace neqone_delta = 1 if labnumber == `i' & delta < 1 & ciu_delta < 1 & se_delta != .
    replace neqone_delta = 1 if labnumber == `i' & delta > 1 & cil_delta > 1 & se_delta != .
    replace neqone_beta  = 1 if labnumber == `i' & beta  < 1 & ciu_beta  < 1 & se_beta != .
    replace neqone_beta  = 1 if labnumber == `i' & beta  > 1 & cil_beta  > 1 & se_beta != .
}

drop time k pratio price_later mnumer cons1 cons_later noise obsnum t0


save "base_estimated_M3_capped.dta", replace
