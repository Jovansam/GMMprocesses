* The model we simulate and estimate here is an Permanent-transitory model where the permanent component is a random walk and the transitory component is a MA(1)
* Formally: y_it = u_{it} + v_{it} where:
	* u_{it} = u_{it-1} + w_{it}
	* v_{it} = e_{it-1} + \theta e_{it}

* I am using data in "long" format so I can make use of lag and difference operators
* Each moment relates to a specific period so is only non-missing for one value of t
* Need to specify "nocommonesample" in the call to gmm to make sure it still runs properly



* Simulate data
***************

clear all
set seed 337631
*set seed 37631

scalar N = 10000
scalar T = 5

scalar sigw = 0.014
scalar sigu1 = 0.866
scalar sige = 0.077
scalar theta = 0.498

set obs `= scalar(N) * scalar(T)'

gen long n = int((_n - 1) / scalar(T)) + 1
gen long t = mod(_n - 1, scalar(T)) + 1

sort n t

* Random walk component
gen double w = rnormal(0.0, scalar(sigw))
qui by n (t): gen double u = rnormal(0.0, scalar(sigu1)) if (_n == 1)
qui by n (t): replace u = u[_n-1] + w if (_n > 1)

* Moving average component
gen double e0 = rnormal(0.0, scalar(sige))
gen double e = rnormal(0.0, scalar(sige))
qui by n (t): gen double v = e + (scalar(theta) * e0) if (_n == 1)
qui by n (t): replace v = e + (scalar(theta) * e[_n-1]) if (_n > 1)

gen double y = u + v

xtset n t



* Moment evaluators
*******************

* Moment evaluator using first-differenced moments
capture program drop gmm_inc_process_first_diff
program gmm_inc_process_first_diff

    syntax varlist if, at(name) incvar(varname) tvar(varname)

	qui xtset
	local tmin = r(tmins)
	local tmax = r(tmaxs)
	local diff = `tmax' - `tmin' + 1
	if (`diff' < 4) {
		di as error "Too few periods (t needs to be at least 4)"
		exit 198
	}

    quietly {
        
		* Calculate differenced outcome variable
		tempvar D`incvar'
		*local D`incvar' "D`incvar'"
		gen double `D`incvar'' = D.`incvar'
		
		* De-mean differenced outcome variable
		tempvar D`incvar'_demean
		*local D`incvar'_demean "D`incvar'_demean"
		egen double `D`incvar'_demean' = mean(`D`incvar'') `if', by(`tvar')
		replace `D`incvar'_demean' = `D`incvar'' - `D`incvar'_demean' `if'

        * Parameters: sigw, sige, theta
		* (Note we exclude sigu1 because we can't estimate it based on first-differenced moments because none depend on it)

		* Calculate non-linear combinations of parameters
		scalar sig2_deltay = (`at'[1,1]^2) + (2.0 * ((`at'[1,3]^2) - `at'[1,3] + 1) * (`at'[1,2]^2))
		scalar cov_deltay_Ldeltay = -((`at'[1,3] - 1)^2) * (`at'[1,2]^2)
		scalar cov_deltay_L2deltay = -`at'[1,3] * (`at'[1,2]^2)

		local varlistcopy "`varlist'"
		
		* Loop across t
		forval t = `=`tmin'+1'/`tmax' {

			* Calculate variance
			gettoken outvar varlist : varlist
			replace `outvar' = ((`D`incvar'_demean')^2) - scalar(sig2_deltay) `if' & (`tvar' == `t')


			* If at least third period, calculate autocorrelation at lag 1
			if (`t' >= `tmin' + 2) {
				gettoken outvar varlist : varlist
				replace `outvar' = (`D`incvar'_demean' * L.`D`incvar'_demean') - scalar(cov_deltay_Ldeltay) `if' & (`tvar' == `t')
			}
			
			* If at least fourth period, calculate autocorrelation at lag 2
			if (`t' >= `tmin' + 3) {
				gettoken outvar varlist : varlist
				replace `outvar' = (`D`incvar'_demean' * L2.`D`incvar'_demean') - scalar(cov_deltay_L2deltay) `if' & (`tvar' == `t')
			}
		
		}
		*pause
		*noi su `varlistcopy'

    }


end 



* Moment evaluator using moments in levels
capture program drop gmm_inc_process_levels
program gmm_inc_process_levels

    syntax varlist if, at(name) incvar(varname) tvar(varname)

	qui xtset
	local tmin = r(tmins)
	local tmax = r(tmaxs)
	local maxdiff = `tmax' - `tmin' + 1
	if (`maxdiff' < 3) {
		di as error "Too few periods (t needs to be at least 3)"
		exit 198
	}

    quietly {
        
		* De-mean outcome variable
		tempvar `incvar'_demean
		*local `incvar'_demean "`incvar'_demean"
		egen double ``incvar'_demean' = mean(`incvar') `if', by(`tvar')
		replace ``incvar'_demean' = `incvar' - ``incvar'_demean' `if'

        * Parameters: sigw, sigu1, sige, theta

		* Calculate non-linear combinations of parameters
		scalar var_v = (1.0 + (`at'[1,4]^2)) * (`at'[1,3]^2)
		scalar cov_v_Lv = `at'[1,4] * (`at'[1,3]^2)
		
		local varlistcopy "`varlist'"
		
		* Loop across t
		forval t = `tmin'/`tmax' {
		
			forval s = `t'/`tmax' {

				local diff = `s' - `t'
				
				scalar cov_vt_vs = 0.0
				if (`diff' == 0) scalar cov_vt_vs = scalar(var_v)
				if (`diff' == 1) scalar cov_vt_vs = scalar(cov_v_Lv)
				
				scalar cov_ut_us = (`at'[1,2]^2) + ((`t' - `tmin') * (`at'[1,1]^2))
			
				gettoken outvar varlist : varlist
				replace `outvar' = (``incvar'_demean'  * F`diff'.``incvar'_demean') - scalar(cov_ut_us) - scalar(cov_vt_vs) `if' & (`tvar' == `t')
				
			}
		
		}

    }


end 



* Estimate models
*****************

local nequations = (scalar(T) - 1) + (scalar(T) - 2) + (scalar(T) - 3)
matrix define initval = (0.1, 0.1, 0.7)
gmm gmm_inc_process_first_diff, nequations(`nequations') parameters(sigw sige theta) winitial(unadjusted, independent) vce(unadjusted) from(initval) twostep nocommonesample incvar(y) tvar(t)

local nequations = 0.5 * scalar(T) * (scalar(T) + 1)
matrix define initval = (0.1, 0.5, 0.1, 0.7)
gmm gmm_inc_process_levels, nequations(`nequations') parameters(sigw sigu1 sige theta) winitial(unadjusted, independent) vce(unadjusted) from(initval) twostep nocommonesample incvar(y) tvar(t)



