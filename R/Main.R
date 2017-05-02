
splitFactorNames = function(n, checkForFormula = TRUE, convertToFormula = TRUE) {
	if (is.null(n)) {
		return(n)
	}
	
	if (convertToFormula && grepl("^\\s*~", n)) {
		n = stats::as.formula(n)
	}
	
	if (checkForFormula && class(n) == class(stats::formula())) {
		return(n)
	}
	if (length(n) == 1 && grepl(":", n, fixed=TRUE)) {
		n = strsplit(n, ":", fixed=TRUE)[[1]]
	}
	n
}

# Maybe external function?
# Renaming the columns doesn't really make sense, because the terms are the terms
# as named as used by the contrasts function. As such, there is no natural renaming,
# other than for contr.treatment (and contr.SAS). As such, no renaming happens 
# except for treatment and SAS contrasts.
makeDesignMatrix = function(factors, dmFactors, contrastType, renameCols=FALSE) {
	
	if (class(dmFactors) == class(stats::formula())) {
		form = dmFactors
	} else {
		dmFactors = splitFactorNames(dmFactors, checkForFormula=FALSE, convertToFormula=FALSE)
		form = paste("~", paste(dmFactors, collapse = " * "))
		form = stats::formula(form)
	}

	tt = stats::terms(form)
	
	o1terms = getFirstOrderTermLabels(tt)
	
	contrasts.arg = list()
	if (is.list(contrastType)) {
		contrasts.arg = contrastType
	} else {
		for (fn in o1terms) {
			contrasts.arg[[fn]] = contrastType
		}
	}
	
	# When using a factor, model.matrix seems to use levels(f) instead of unique(f).
	# When a factor has more levels than are present in the factor, this is bad.
	# It also doesn't like numbers that are not in a factor.
	for (n in names(factors)) {
		factors[ , n] = as.character(factors[ , n])
	}
	
	m = stats::model.matrix(form, factors, contrasts.arg)

	if (renameCols && is.character(contrastType) && contrastType %in% c("contr.treatment", "contr.SAS")) {
		colnames(m) = renameDesignMatrixColumns( colnames(m), factors, dmFactors, contrastType)
	}
	
	list(mat=m, terms=tt)
}

#This function is only used once. Consider inlining it.
getFirstOrderTermLabels = function(terms) {
	attr(terms, "term.labels")[ attr(terms, "order") == 1 ]
}

# Maybe external function? Should be external if makeDesignMatrix is external.
# If the design matrix is not full rank, this function removes columns 
# that contribute to it being not full rank.
stripExcessTermsFromDM = function(dm) {
	
	mat = dm$mat
	
	q = qr(mat)
	
	kept = q$pivot[ seq(q$rank) ]
	dropped = q$pivot[ !(seq(ncol(q$qr)) %in% seq(q$rank)) ]
	
	allAssignments = attr(mat, "assign")
	keptAssignments = unique(allAssignments[ kept ])
	fullyDroppedAssignments = unique(allAssignments[ !(allAssignments %in% keptAssignments) ])
	fullyDroppedTerms = attr(dm$terms, "term.labels")[fullyDroppedAssignments]
	
	newMat = subset(mat, select = kept)
	
	attr(newMat, "assign") = attr(mat, "assign")[ kept ]
	attr(newMat, "contrasts") = attr(mat, "contrasts")
	
	list(mat=newMat, kept = colnames(mat)[kept], dropped = colnames(mat)[dropped], fullyDroppedTerms = fullyDroppedTerms)
	
}


# Internal function?
# Gets the columns of the design matrix associated with particular effects.
# Those columns are the same as the elements of the effect parameters.
getEffectAssignmentColumns = function(dm, fNames) {
	
	if (length(fNames) == 1 && fNames == "(Intercept)") {
		return(which(attr(dm$mat, "assign") == 0))
	}
	
	labels = attr(dm$terms, "term.labels")
	
	#Select first order (main effect) terms.
	o1labels = labels[ attr(dm$terms, "order") == 1 ]
	
	if (!all(fNames %in% o1labels)) {
		stop("fNames and design matrix term labels do not match. Make sure that the design matrix has all of the terms that you want in it.")
	}
	
	fNames = o1labels[ o1labels %in% fNames ] # sort fNames by label names
	
	effect = paste0(fNames, collapse=":")
	
	assignments = attr(dm$mat, "assign")
	
	effectAssignment = which(labels == effect)

	mCols = which(assignments == effectAssignment)
	
	if (length(mCols) == 0) {
		stop( paste0("No design matrix columns associated with effect \"", effect, "\".") )
	}
	
	mCols
	
}



#' Calculate Effect Parameter Matrix
#' 
#' From a matrix of cell means, the factors in the experiment, and information about 
#' which factor to test, creates a matrix of ANOVA effect parameters.
#' 
#' The steps of the procedure:
#' 1. Calculate the design matrix, `X`.
#' 1b. If design is not fully crossed, strip excess terms from `X`.
#' 2. If `mu` are cell means and `beta` are effect parameters, then `mu = X * beta` and `beta = (X' * X)^-1 * X' * mu`. Calculate `S = (X' * X)^-1 * X'`.
#' 3. Calculate `beta = S * mu`.
#' 4. Select out only the `beta` needed for the `testedFactors`, `beta_s`.
#' 5. Select the corresponding columns (and rows) in `X`, `X_s`.
#' 6. Complete the required effect parameters with `X_s * beta_s`. This is not the same as `mu = X * beta` because only subsets of `X` and `beta` are used.
#' 
#' @seealso See also [`summarizeEffectParameters`] and [`groupEffectParameters`].
#'
#' @param cellMeans A matrix of cell means. Each column is a cell and each row is a sample from the prior or posterior. The columns must correspond to the rows of `factors`.
#' @param factors See [`testHypothesis`].
#' @param testedFactors See [`testHypothesis`].
#' @param dmFactors See [`testHypothesis`].
#' @param contrastType See [`testHypothesis`].
#' @param warnOnDrop Emit a warning if columns of the design matrix are dropped due to it not being full rank.
#'
#' @return A matrix of effect parameters where each column is one parameter and each row is a sample from the prior or posterior (depending on what `cellmeans` was). The columns are named with factor level names.
#'
#' @md
#' @export 
getEffectParameters = function(cellMeans, factors, testedFactors, dmFactors = NULL, 
															 contrastType = NULL, warnOnDrop=FALSE) {
	
	if (is.null(dmFactors)) {
		if (length(testedFactors) == 1 && testedFactors == "(Intercept)") {
			dmFactors = stats::formula(" ~ 1") #only use intercept
		} else {
			dmFactors = testedFactors
		}
	}

	testedFactors = splitFactorNames(testedFactors, convertToFormula = FALSE)
	dmFactors = splitFactorNames(dmFactors)
	
	fullyCrossed = isDesignFullyCrossed(factors, warnOnDuplicate=FALSE)
	if (is.null(contrastType)) {
		if (fullyCrossed) {
			contrastType = "contr.sum"
			if (getOption("verbose")) {
				cat("The design is fully crossed and sums-to-zero contrasts (contr.sum) have been chosen.\n")
			}
		} else {
			contrastType = "contr.treatment"
			if (getOption("verbose")) {
				cat("The design is not fully crossed and treatment contrasts (contr.treatment) have been chosen.\n")
			}
		}
	} else {
		if (!fullyCrossed && contrastType %in% c("contr.sum", "contr.helmert", "contr.poly")) {
			warning("The design is not fully crossed but you are using orthogonal contrasts. Contrasts cannot actually be orthogonal if the design is not fully crossed.")
		}
	}
	
	# 1. Calculate the design matrix, X.
	dm = makeDesignMatrix(factors, dmFactors, contrastType, renameCols=FALSE)

	# 1b. If design is not fully crossed, strip excess terms from X.
	strippedInfo = stripExcessTermsFromDM(dm)
	for (fdt in strippedInfo$fullyDroppedTerms) {
		parts = strsplit(fdt, ":", fixed=TRUE)[[1]]
		if (length(testedFactors) == length(parts) && all(testedFactors %in% parts)) {
			stop("Unable to get effect parameters for testedFactors \"", paste(testedFactors, collapse=":"), "\" because all terms related to that effect have been stripped from the design matrix. This error happens when your design is lacking the right cells to allow you to test this effect. You likely have some kind of unbalanced design. See the \"Non-Fully-Crossed/Unbalanced Designs\" section of the CMBBHT package manual.")
		}
	}
	
	if (warnOnDrop && length(strippedInfo$dropped) > 0) {
		warning( paste0("The design is not fully crossed (unbalanced). As a result, some terms were dropped from the design matrix: ", paste(strippedInfo$dropped, collapse=", "), ". The following terms were kept: ", paste(strippedInfo$kept, collapse=", "), ". The naming of these terms depends on the contrastType that was used.") )
	}
	
	X = dm$mat = strippedInfo$mat

	
	# 2. Calculate S = (X'X)^-1 X'
	S = solve(t(X) %*% X) %*% t(X)
	
	# 3. Calculate beta = S %*% mu
	# Some extra transposition because mu and beta are both transposed
	beta = t(S %*% t(cellMeans))
	
	# 4. Select out only the needed beta, beta_s
	mCols = getEffectAssignmentColumns(dm, testedFactors)
	
	beta_s = subset(beta, select = mCols)
	
	# 5. Select the corresponding columns (and rows) in X, X_s
	X_s = subset(X, select = mCols)
	
	#select only some rows of X_s
	
	if (length(testedFactors) == 1 && testedFactors == "(Intercept)") {
		X_s_rows = 1
	} else {
	
		usedFactorLevels = unique(subset(factors, select = testedFactors))
		X_s_rows = rep(NA, nrow(usedFactorLevels))
		
		for (i in 1:nrow(usedFactorLevels)) {
			
			rows = NULL
			for (j in 1:nrow(factors)) {
				if (all(factors[j, testedFactors] == usedFactorLevels[i,])) {
					rows = c(rows, j)
				}
			}
			# For duplicate rows, the X elements are the same, so just take the first.
			X_s_rows[i] = rows[1] 
		}
	}
	
	X_s = X_s[ X_s_rows, ]
	
	# 6. Complete the betas by X_s %*% beta_s (or a funny transposed version)
	fullEffects = beta_s %*% t(X_s)
	
	# Name the columns with cell names
	if (length(testedFactors) == 1 && testedFactors == "(Intercept)") {
		colnames(fullEffects) = "(Intercept)"
	} else {
		colnames(fullEffects) = makeCellName(usedFactorLevels)
	}
	
	fullEffects
}


#' Perform Hypothesis Test from Cell Means
#' 
#' This is the primary function in this package. It takes two matrices of prior and posterior cell means, a `data.frame` containing information about the factor levels in the design, and what effect you want to test. It computes a Bayes factor related to the hypothesis that the effect is present in the cell means.
#' 
#' @param priorCMs Numeric matrix. Cell means sampled from the priors. The columns must correspond to the rows of factors but do not need to be named.
#' @param postCMs Numeric matrix. Cell means sampled from the posterior distribution. The columns must correspond to the rows of factors but do not need to be named. 
#' @param factors A `data.frame` containing information about the experimental design. Each column is a factor of the design. Each row contains the levels of the factors that define a cell of the design. No additional columns may be included in factors. Factor names and factor levels must not include period (".") or colon (":").
#' @param testedFactors Character vector. The factors for which to perform the hypothesis test as a vector of factor names. A single factor name results in the test of the main effect of the factor. Multiple factor names result in the test of the interaction of all of those factors. You may provide either a vector with multiple elements, e.g. `c('A', 'B')`, or a vector with one element where factor names are separated by colon, e.g. `'A:B'`.
#' @param dmFactors Character vector or formula. The factors to use to construct the design matrix. Like `testedFactors`, you may separate factor names with colon. For a fully-crossed (balanced) design using orthogonal contrasts, this can always be equal to `testedFactors` (the default). For non-fully-crossed designs, you may sometimes want to create a design matrix using a set of factors, but perform a hypothesis test with only some of those factors (`testedFactors` must be a subset of `dmFactors`). This constraint is not tested for if a formula is provided. You may supply a formula like that taken by [`model.matrix`] which will be used to create the design matrix. The formula should be like ` ~ A * B`, where A and B are factor names with nothing on the left hand side of the "~".
#' @param contrastType Character, function, or list. The contrast to use to create the design matrix. If character, can be any of the function names on the documentation page for `contr.sum`. For a non-fully-crossed (unbalanced) design, you should use either `"contr.treatment"` or `"contr.SAS"`. For a balanced design, you can use anything, but psychologists are most used to `"contr.sum"`, which uses sums-to-zero constraints. If a function, it should produce contrasts. If a list, it should be able to be passed directly to the `contrasts.arg` argument of [`stats::model.matrix`].
#' @param testFunction A function that takes two matrices of prior and posterior effect parameters, in that order. For example, see [`testFunction_SDDR`]. You can probably leave this at the default value.
#' @param usedFactorLevels A `data.frame` with a column for each of the factors in `testedFactors`. Each row specifies factor levels that should be included in the test. This allows you to do things like pairwise comparisons of specific factor levels. The factor levels that are not in `usedFactorLevels` are dropped after calculation of the effect parameters, which means that the kept effect parameters are calculated in the context of any effects that are specified by `dmFactors`.
#' 
#' @return The return value depends on the choice of `testFunction`. See [`testFunction_SDDR`] for an example.
#' 
#' @md
#' @export
testHypothesis = function(priorCMs, postCMs, factors, testedFactors, dmFactors = testedFactors,
													contrastType = NULL, testFunction = testFunction_SDDR, usedFactorLevels = NULL) {
	
	#TODO: Do you want to do this? It's a minor burden on users to specify.
	#Due to lazy evaluation, dmFactors is set to testedFactors after this happens.
	if (missing(testedFactors) || is.null(testedFactors)) {
		if (ncol(factors) == 1) {
			testedFactors = names(factors)
		} else {
			stop("testedFactors is missing or NULL.")
		}
	}
	
	testedFactors = splitFactorNames(testedFactors, convertToFormula = FALSE)
	dmFactors = splitFactorNames(dmFactors)
	
	if (nrow(factors) != ncol(priorCMs) || nrow(factors) != ncol(postCMs)) {
		stop("The number of rows in factors does not correspond to the number of columns in the cell means arguments (priorCMs or postCMs).")
	}
	
	if (class(dmFactors) != class(stats::formula()) && !all(testedFactors %in% dmFactors)) {
		stop("testedFactors must be a subset of dmFactors.")
	}
	
	priorEffects = getEffectParameters(priorCMs, factors, testedFactors, 
																	 dmFactors=dmFactors, contrastType=contrastType)
	
	postEffects = getEffectParameters(postCMs, factors, testedFactors, 
																	dmFactors=dmFactors, contrastType=contrastType)
	
	
	if (!is.null(usedFactorLevels)) {
		cns = makeCellName(usedFactorLevels)
		
		priorEffects = priorEffects[ , cns ]
		postEffects = postEffects[ , cns ]
	}
	
	#TODO: Keep these errors?
	if (ncol(priorEffects) != ncol(postEffects)) {
		stop("The number of columns (parameters) in priorEffects does not equal the number of columns in postEffects. See the documentation for this function. This is probably the result of a conceptual error.")
	}
	
	if (ncol(priorEffects) < 2 || ncol(postEffects) < 2) {
		stop("priorEffects or postEffects has less than two columns. You must have at least two effect parameters to perform a hypothesis test of main effects or interactions.")
	}
	
	testFunction(priorEffects, postEffects)
}

#' Test Multiple Hypotheses
#' 
#' Convenience function for testing multiple hypotheses. E.g., for a two factor design, you could test the main effects of A and B plus the interaction of A and B. If using `testHypothesis` you would need to call that function three times, but this function requires only one function call.
#' 
#' @param prior See [`testHypothesis`].
#' @param post See [`testHypothesis`].
#' @param factors See [`testHypothesis`].
#' @param testedFactors Character vector (or list). Interactions should be indicated by putting colons between factor names. For example, the interaction of A and B is given by "A:B". The order of factor names does not matter. If a list, the elements should not be named.
#' @param dmFactors Character vector (or list) of the same length as `testedFactors`.
#' @param contrastType See [`testHypothesis`].
#' @param testFunction See [`testHypothesis`].
#' @param usedFactorLevels List of data frames. If provided, should be the same length as `testedFactors`.
#' @param testName Character vector (or list). An optional name for the tests. If provided, should be the same length as `testedFactors`.
#' 
#' @return A `data.frame` with one test on each row.
#' 
#' @md
#' @export 
testHypotheses = function(prior, post, factors, testedFactors, dmFactors = testedFactors, contrastType = NULL, testFunction = testFunction_SDDR, usedFactorLevels = NULL, testName = testedFactors) {
	
	#If testedFactors is not provided, use all factors.
	if (missing(testedFactors) || is.null(testedFactors)) {
		ns = names(factors)
		testedFactors = list()
		tfi = 1
		for (i in 1:length(ns)) {
			cb = utils::combn(ns, i)
			for (j in 1:ncol(cb)) {
				testedFactors[[tfi]] = cb[,j]
				tfi = tfi + 1
			}
		}
	}
	
	if (!is.list(testedFactors)) {
		testedFactors = as.list(testedFactors)
	}
	if (!is.list(dmFactors)) {
		dmFactors = as.list(dmFactors)
	}
	if (!is.list(testName)) {
		testName = as.list(testName)
	}
	
	if (length(testedFactors) != length(dmFactors) || (!is.null(usedFactorLevels) && length(testedFactors) != length(usedFactorLevels))) {
		stop("One of testedFactors, dmFactors, and usedFactorLevels has a different length than the others.")
	}
	
	allNames = NULL
	allRes = list()
	for (i in 1:length(testedFactors)) {

		temp = list(testName = paste(testName[[i]], collapse=":"))
		
		ht = tryCatch({
			testHypothesis(prior, post, factors, testedFactors[[i]], dmFactors = dmFactors[[i]], contrastType = contrastType, testFunction = testFunction, usedFactorLevels = usedFactorLevels[[i]])
		}, error = function(e) {
			print(e)
			list(success=FALSE)
		})
		
		allRes[[i]] = c(temp, ht)
		allNames = union(allNames, names(allRes[[i]]))
	}
	
	res = NULL
	for (i in 1:length(allRes)) {
		#Add NAs where needed
		for (n in allNames) {
			if (!(n %in% names(allRes[[i]]))) {
				allRes[[i]][[n]] = NA
			}
		}

		tr = as.data.frame(allRes[[i]], stringsAsFactors = FALSE)
		res = rbind(res, tr)
	}
	
	res
}

#' Test the Value of the Intercept
#' 
#' Tests whether the intercept/grand mean has some given value. For a general test of whether a parameter has a given value, see [`valueTest_SDDR`]. Note that in order for this test to be very meaningful (probably)
#' 
#' @param priorCMs See [`testHypothesis`].
#' @param postCMs See [`testHypothesis`].
#' @param factors See [`testHypothesis`].
#' @param testVal The value that will be tested. Depends on the choice of `testFunction`. See [`valueTest_SDDR`] for the interpretation of this argument by the default test function.
#' @param dmFactors See [`testHypothesis`].
#' @param contrastType See [`testHypothesis`].
#' @param testFunction A function of 3 positional arguments: the prior intercept, the posterior intercept, and the point at which the hypothesis is tested (i.e. `testVal` is passed as the third argument). A custom test function may choose to ignore `testVal` if the tested values are determined in some other way. The prior and posterior intercepts will be vectors or one-column matrices.
#' 
#' @return The result of the function depends on the choice of `testFunction`. See [`valueTest_SDDR`] for an example.
#' 
#' @md
#' @export
testIntercept = function(priorCMs, postCMs, factors, testVal, dmFactors = stats::formula(" ~ 1"), contrastType = NULL, testFunction=valueTest_SDDR) {

	priorInt = getEffectParameters(priorCMs, factors, 
									 testedFactors = "(Intercept)", dmFactors=dmFactors, 
									 contrastType=contrastType)
	
	postInt = getEffectParameters(postCMs, factors, 
				  				testedFactors = "(Intercept)", dmFactors=dmFactors, 
					  			contrastType=contrastType)
	

	testFunction(priorInt, postInt, testVal)
	
}


