

# Internal/testing function. DO NOT USE.
#TODO: This seems bizarre. What is the calibration that people are supposed to use?
testFunction_credibleInterval = function(prior, posterior, CIP = 0.95) {
	
	cips = c((1 - CIP) / 2, (1 + CIP) / 2)
	
	qs = as.numeric(stats::quantile(posterior, cips))
	
	list(lower=qs[1], upper=qs[2])
	
}

# Internal/testing function. DO NOT USE.
#I'm still not sure about this one working properly very much of the time.
#Also, choice of p is subjective and affects results.
testFunction_lowerTail = function(prior, posterior, p = 0.01) {
	
	qprior = stats::quantile(prior, p)
	qpost = stats::quantile(posterior, p)
	
	# Make sure that you have at least p in the lower tail of both distributions
	q = max(c(qprior, qpost))
	
	pbPrior = mean(prior < q)
	pbPost = mean(posterior < q)
	
	bfBelow = ((1 - pbPrior) / pbPrior) * (pbPost / (1 - pbPost))
	
	bfAbove = 1 / bfBelow
	
	list(bf01 = bfBelow, bf10 = bfAbove, p = p, q = q, success=TRUE)
}

# Sort of internal function.
# only appropriate if x has length 2, but this is not tested for.
devianceFunction_absDif = function(x) {
	abs(x[1] - x[2])
}

# Internal function
# To deal with the fact that logspline doesn't like NULL arguments for lbound and ubound.
logspline_null = function(x, lbound=NULL, ubound=NULL) {
	
	ls = NULL
	
	if (is.null(lbound)) {
		if (is.null(ubound)) {
			ls = polspline::logspline(x)
		} else {
			ls = polspline::logspline(x, ubound=ubound)
		}
	} else {
		if (is.null(ubound)) {
			ls = polspline::logspline(x, lbound=lbound)
		} else {
			ls = polspline::logspline(x, lbound=lbound, ubound=ubound)
		}
	}
	
	ls
}




#' Test Function: Savage-Dickey Density Ratio
#' 
#' Takes matrices of prior and posterior samples the effect parameters. 
#' Calculates a deviance measure, `D`, for the prior and posterior effect parameters. 
#' A deviance value of 0 happens if there is no effect. Thus, the tested hypotheses are 
#' H_0: `D = 0` and H_1: `D =/= 0`. 
#' Returns Bayes factors related to both the null and alternative hypotheses. 
#' See the "Deviance Measure" and "Bayes Factor Estimation with Savage-Dickey" 
#' sections of the manual for more information about this procedure.
#' 
#' Due to a limitation of the density estimation procedure, the matrices of samples fom the 
#' prior and posterior should have the same number of rows. You will get a warning about this 
#' if they are not the same length. Given that sampling from the prior is usually relatively 
#' easy, once you have samples from the posterior, sample the same amount from the prior.
#' 
#' It is assumed that the prior may be diffuse, such as a Cauchy prior.
#' In that case, it is possible to have prior values of `D` that are very far from the
#' largest posterior value (like 10^9 times farther). This hurts the ability of the
#' density estimation to estimate the density in the same way for the prior and
#' posterior, which could result in a bias. This function accounts for unusually large
#' prior `D` values by only using prior `D` at most 10 times larger than the largest 
#' posterior `D`.
#' 
#' Note that if this function is failing with errors from the density estimation procedure,
#' you should consider changing the defaults for `postMaxMult` and `min_pKept`. 
#' It is possible to have cases where the prior is about as diffuse as the posterior, except
#' for extreme values of `D`. In that case, using values from the prior substantialy larger than the
#' median of the prior can result in the left edge of the prior being very compacted and
#' the density estimation can't deal well with that. You can change the default values with 
#' a curried function that then gets passed to, e.g., [`testHypothesis`]. See the examples.
#' 
#' 
#' @param priorEffects Numeric matrix of effect parameters sampled from the priors. Each column is one parameter and each row is one iteration. It must have at least two columns.
#' @param postEffects Numeric matrix of effect parameters sampled from the posteriors. Must have the same number of columns as `priorEffects`. Mathematically, this procedure works even with different numbers of columns, but such a test would be bizarre and meaningless (some parameters didn't have priors?). If the number of columns is not the same, an error will be emitted.
#' @param devianceFunction A function used for calculating the deviation of the effect parameters, `D`. It takes a vector of effect parameters and calculates some measure of how dispersed they are. One example of such a function is the sample variance (see [`stats::var`]).
#' @param postMaxMult The prior is cut off above `postMaxMult * max(post_D)`, where `post_D` is the posterior deviance measure. See also `min_pKept`.
#' @param min_pKept The minimum proportion of the prior that will be kept. Overrides `postMaxMult` if the `postMaxMult` rule would keep less than `min_pKept` of the prior. This rule can be disabled if `min_pKept` is set to `NULL` or `1`.
#' @param truncatePosterior If `TRUE`, the posterior will be truncated to the same proportion as the prior. This is generally a good thing as it further equates the treatment of the prior and the posterior and keeps the number of prior and posterior samples the same. If `FALSE`, the posterior will not be altered.
#' @param warnOnLength If `TRUE` and the `prior` and `posterior` are not the same length, a warning will be emitted.
#'
#' @return A list with four elements: 
#' * `success`: A boolean indicating whether there was an exception during density estimation. If `success` is `FALSE`, all other values will be `NULL` or `NA`.
#' * `bf01`: The Bayes factor in favor of H0.
#' * `bf10`: The Bayes factor in favor of H1. 
#' * `pKept`: The proportion of the prior distribution that was used to estimate the density;
#' 
#' @md
#' @export
#' 
#' @examples
#' \dontrun{
#' curriedTestFunction = function(priorEffects, postEffects) { 
#'   testFunction_SDDR(priorEffects, postEffects, 
#'     devianceFunction = sd, postMaxMult = 1.5, min_pKept = 0.95)
#' }
#' 
#' testHypothesis(..., testFunction = curriedTestFunction)
#' }
testFunction_SDDR = function(priorEffects, postEffects, devianceFunction = NULL,
														 postMaxMult = 10, min_pKept = 0.90, truncatePosterior = TRUE,
														 warnOnLength = TRUE) 
{
	
	if (is.null(devianceFunction)) {
		if (ncol(priorEffects) == 2) {
			devianceFunction = devianceFunction_absDif
		} else {
			devianceFunction = stats::var
		}
	}
	
	prior_D = apply(priorEffects, 1, devianceFunction)
	post_D = apply(postEffects, 1, devianceFunction)
	
	if (warnOnLength && nrow(priorEffects) != nrow(postEffects)) {
		warning("The length of the prior and posterior should be the same for accurate density estimation.")
	}
	
	# If the prior extends to very large values, cut it off.
	# This helps with density estimation at smaller values.
	if (is.null(min_pKept)) {
		mpkq = 0
	} else {
		mpkq = stats::quantile(prior_D, min_pKept)
	}
	
	postmmq = postMaxMult * max(post_D)
	
	priorMax = max(mpkq, postmmq)
	
	# Only truncate if there will be truncation
	if (priorMax < max(prior_D)) {
		keptPriorValues = (prior_D <= priorMax)
		pKept = mean(keptPriorValues)
		prior_D = prior_D[ keptPriorValues ]
		
		if (truncatePosterior) {
			postMax = stats::quantile(post_D, pKept)
			post_D = post_D[ post_D < postMax ]
		}
		
	} else {
		pKept = 1
		priorMax = NULL
		postMax = NULL
	}
	
	
	#Do the prior and posterior estimation in a tryCatch because 
	#failures can happen during density estimation.
	success = tryCatch({
		priorLS = logspline_null(prior_D, lbound = 0, ubound = priorMax)
		postLS = logspline_null(post_D, lbound = 0, ubound = postMax)
		TRUE
	}, error = function(e) {
		print(e)
		return(FALSE)
	})
	
	
	if (success) {
		
		priorDens = polspline::dlogspline(0, priorLS)
		postDens = polspline::dlogspline(0, postLS)
		
		# Account for the fact that there is some density in the upper area 
		# that isn't accounted for by the logspline. Scale down the density
		# in the lower area because the upper area was not included and the 
		# lower area appears to have more data in it than it really does, 
		# which leads to an overestimate of the density.
		priorDens = priorDens * pKept
		if (truncatePosterior) {
			postDens = postDens * pKept
		}
		
		# Savage-Dickey step
		bf10 = priorDens / postDens
		
		rval = list(success = TRUE, bf01 = 1 / bf10, bf10 = bf10, pKept = pKept)
	} else {
		rval = list(success = FALSE, bf01 = NA, bf10 = NA, pKept = NA)
	}
	
	rval
}


#' Test Function: Encompassing Prior
#' 
#' See the "Encompassing Prior Approach" section of the package manual for an explanation of this approach to estimating Bayes factors. There is some encompassing model, `M_1`, and a constrained model, `M_0`, which is created by placing some constraint on some of the parameters of `M_1`. The function `I_M0` that is passed to this function determines whether the constraint is satisfied.
#'  
#' Note that you cannot pass this function directly as the `testFunction` argument of [`testHypothesis`] because there is no default value for `I_M0`. Thus, you must create a curried function to pass to [`testHypothesis`] or use [`create_EPA_intervalTF`]. See the examples for an example of currying.
#' 
#' @param priorEffects A numeric matrix of prior effect parameters. Rows are samples and columns are parameters.
#' @param postEffects A numeric matrix of posterior effect parameters. Rows are samples and columns are parameters.
#' @param I_M0 The indicator function for the constrained model, `M_0`. This is a function that takes a vector of effect parameters, checks some constraint on those parameters, and returns 1 (or `TRUE`) if the constraint is satisfied or 0 (or `FALSE`) if the constraint is not satisfied.
#' 
#' @return A list for the following elements:
#' * `success`: Logical. Whether or not the test was successful. This is `FALSE` if the Bayes factors are infinte, NaN, or NA, or `TRUE` otherwise.
#' * `bf10`: The Bayes factor in favor of the encompassing/alternative model, `M_1`.
#' * `bf01`: The Bayes factor in favor of the constrained/null model, `M_0`.
#' * `prior_satisfied`, `post_satisfied`: The proportions of the prior and posterior, respectively, that satisfy the constraint of `I_M0`. If these are both low and the number of samples satisfying the constraint is small, it indicates that your constraint is possibly too restrictive and that the estimated Bayes factor may be noisy. You should either use a less-restrictive constraint or take more prior and posterior samples.
#' 
#' @md
#' @export
#' 
#' @examples \dontrun{
#' curriedTestFun = function(priorEffects, postEffects) {
#'   #M_0: Constrain the effects to be less than 2 units from 0
#'   I_M0 = function(eff) {
#'     all(abs(eff) < 2)
#'   }
#'   testFunction_EPA(priorEffects, postEffects, I_M0)
#' }
#' 
#' testHypothesis(..., testFunction = curriedTestFun)
#' }
testFunction_EPA = function(priorEffects, postEffects, I_M0) {
	
	# Calculate the numerator: The average I_M0 for the posterior
	post_sat = apply(postEffects, 1, I_M0)
	numerator = mean(post_sat)
	
	# Calculate the denominator: The average I_M0 for the prior
	prior_sat = apply(priorEffects, 1, I_M0)
	denominator = mean(prior_sat)
	
	
	bf_01 = numerator / denominator #In favor of the constrained hypothesis/model
	bf_10 = 1 / bf_01 # In favor of the general hypothesis/model
	
	res = list(success = TRUE, bf01 = bf_01, bf10 = bf_10, prior_satisfied = denominator, post_satisfied = numerator)
	if (is.infinite(bf_01) || is.infinite(bf_10) || is.nan(bf_01) || is.na(bf_01)) {
		res$success = FALSE
	}
	res
}

#' Create Encompassing Prior Interval Test Function
#' 
#' Convenience function for using the encompassing prior test function in the case when an interval is to be tested. The null hypothesis is that all of the effect parameters are within the interval.
#' 
#' @param lower The lower end of the interval to be used. Single-sided intervals can be constructed by using `-Inf`.
#' @param upper The upper end of the interval to be used. Single-sided intervals can be constructed by using `Inf`.
#' 
#' @return A function that can be passed as the `testFunction` argument of, e.g., [`testHypothesis`].
#' 
#' @md
#' @export
#' @examples
#' \dontrun{
#' tf = create_EPA_intervalTF(-0.5, 0.5)
#' testHypothesis(prior, post, fact, "f", testFunction = tf)
#' }
create_EPA_intervalTF = function(lower, upper) {
	
	rval = local({
		force(lower)
		force(upper)
		function(prior, post) {
			I_M0 = function(eff) {
				all(eff >= lower & eff <= upper)
			}
			testFunction_EPA(prior, post, I_M0)
		}
	})
	
	rval
}

#' Test Specific Parameter Value
#' 
#' Performs a hypothesis test of whether a single parameter, `P`, has a given value, with the value chosen with `testVal`. Uses the Savage-Dickey density ratio.
#' 
#' This function assumes that the prior and posterior are both unbounded and potentially diffuse. To account for diffuse priors, it truncates both the prior and posterior in some way. By default, it keeps some proportion of the prior and posterior, which is set by the `pKept` argument. You use some combination of the `pKept` and `bounds` arguments to set the bounds.
# It assumes that the density at the truncation points is near 0.
#' 
#' @param prior A vector of samples from the prior of the parameter.
#' @param posterior A vector of samples from the posterior of the parameter.
#' @param testVal The value that will be tested. The null hypothesis is that `P == testVal` and the alternative hypothesis is that `P =/= testVal`.
#' @param pKept The proportion of the prior and posterior distributions that will be kept. The rest will be discarded from the tails.
#' @param bounds If `TRUE`, bounds are based on `pKept`. If `FALSE`, no bounds are used. If a length 2 numeric vector, those are the bounds that are used.
#' 
#' @return A list with several elements:
#' * `success`: Whether density estimation was successful. If `FALSE`, all of the other values will be `NA` or `NULL`.
#' * `bf10`: The Bayes factor in favor of the hypothesis that `P == testVal`.
#' * `bf01`: The Bayes factor in favor of the hypothesis that `P =/= testVal`.
#' * `prior_pKept`: The actual proportion of the prior that was kept. Should usually be equal to `pKept`.
#' * `post_pKept`: The actual proportion of the posterior that was kept. Should usually be equal to `pKept`.
#' 
#' @md
#' @export
valueTest_SDDR = function(prior, posterior, testVal, pKept = 0.96, bounds=TRUE) {
	
	if (length(prior) != length(posterior)) {
		warning("The length of the prior and posterior should be the same for accurate density estimation.")
	}
	
	cips = c((1 - pKept) / 2, (1 + pKept) / 2)
	
	prior_bounds = stats::quantile(prior, cips)
	post_bounds = stats::quantile(posterior, cips)
	
	if (is.numeric(bounds) && length(bounds) == 2) {
		post_bounds = prior_bounds = bounds
	} else if (is.logical(bounds)) {
		if (bounds == FALSE) {
			post_bounds = prior_bounds = c(-Inf, Inf)
		}
	} else {
		stop("Invalid value of the bounds argument.")
	}
	
	if (testVal < prior_bounds[1] || testVal > prior_bounds[2] || testVal < post_bounds[1] || testVal > post_bounds[2]) {
		warning("The testVal is outside of the bounds. The bounds have been adjusted so that the testVal is within bounds.")
		if (testVal < prior_bounds[1] && testVal < prior_bounds[2]) {
			prior_bounds[1] = testVal
		}
		if (testVal > prior_bounds[1] && testVal > prior_bounds[2]) {
			prior_bounds[2] = testVal
		}
		if (testVal < post_bounds[1] && testVal < post_bounds[2]) {
			post_bounds[1] = testVal
		}
		if (testVal > post_bounds[1] && testVal > post_bounds[2]) {
			post_bounds[2] = testVal
		}
	}
	
	priorKept = prior > prior_bounds[1] & prior < prior_bounds[2]
	postKept = posterior > post_bounds[1] & posterior < post_bounds[2]
	
	priorPKept = mean(priorKept)
	postPKept = mean(postKept)
	

	# If bounds == FALSE, set bounds to NULL. This is needed because non-null bounds are needed just above
	if (is.logical(bounds) && bounds == FALSE) {
		post_bounds = prior_bounds = NULL
	}
	
	success = tryCatch({
		priorLS = logspline_null(prior[ priorKept ], lbound=prior_bounds[1], ubound=prior_bounds[2])
		postLS = logspline_null(posterior[ postKept ], lbound=post_bounds[1], ubound=post_bounds[2])
		TRUE
	}, error = function(e) {
		print(e)
		return(FALSE)
	})
	
	if (success) {
		
		priorDens = polspline::dlogspline(testVal, priorLS)
		postDens = polspline::dlogspline(testVal, postLS)
		
		#Account for dropped data
		priorDens = priorDens * priorPKept
		postDens = postDens * postPKept
		
		bf10 = priorDens / postDens
		
		rval = list(success = TRUE, bf01 = 1 / bf10, bf10 = bf10, prior_pKept = priorPKept, post_pKept = postPKept)
	} else {
		rval = list(success = FALSE, bf01 = NA, bf10 = NA, prior_pKept = NA, post_pKept = NA)
	}
	
	rval
}
