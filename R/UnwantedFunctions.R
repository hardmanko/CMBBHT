

# Unwanted function, but is used by makeDesignMatrix.
#This function may not be needed (probably shouldn't be needed).
#fNames should be enough fNames to account for everything in oldNames
# contrastType may only be scalar character
renameDesignMatrixColumns = function(oldNames, factors, fNames, contrastType) {
	
	newNames = rep(NA, length(oldNames))
	
	for (i in 1:length(oldNames)) {
		
		if (oldNames[i] == "(Intercept)") {
			newNames[i] = oldNames[i]
			next
		}
		
		parts = strsplit(oldNames[i], split = ":", fixed=TRUE)[[1]]
		
		factorsLevels = list()
		
		for (j in 1:length(parts)) {
			dlm = ":"
			
			#Find the factor that starts the name of this part of this term
			for (fn in fNames) {
				if (substr(parts[j], 1, nchar(fn)) == fn) {
					break
				}
			}
			
			level_something = gsub(fn, "", parts[j])
			
			if (contrastType == "contr.sum" || contrastType == "contr.helmert") {
				levInd = as.integer(level_something)
				
				lev = unique(factors[, fn])[levInd]
			} else if (contrastType == "contr.treatment" || contrastType == "contr.SAS") {
				lev = level_something
			} else {
				stop("Invalid contrast type.")
			}
			
			factorsLevels[[fn]] = lev
		}
		
		newNames[i] = makeCellName(factorsLevels)
		
	}
	
	newNames
}

# Unwanted function
# assumes sums to 0 and maybe fully crossed design(?)
# unfilled values in fv should be NA.
# can be used with betas or weights
fillOutFullValues_sumsToZero = function(fv, uniqueFL) {
	for (col in colnames(fv)) {
		if (all(!is.na(fv[, col]))) {
			next #if this one is known, do nothing
		}
		
		fl = splitCellName(col)
		for (variedF in names(fl)) {
			variedLevels = unique(uniqueFL[,variedF])
			variedLevels = variedLevels[ variedLevels != fl[[variedF]] ]
			
			othersKnown = rep(FALSE, length(variedLevels))
			otherCols = rep("", length(variedLevels))
			
			for (vli in 1:length(variedLevels)) {
				flc = fl
				flc[[variedF]] = variedLevels[vli]
				variedCol = makeCellName(flc)
				
				othersKnown[vli] = all(!is.na(fv[, variedCol]))
				otherCols[vli] = variedCol
				
			}
			
			if (all(othersKnown)) {
				#we're in business
				
				otherColVals = fv[ , otherCols ]
				if (length(otherCols) == 1) {
					fv[ , col ] = -otherColVals
				} else {
					fv[ , col ] = apply(otherColVals, 1, function(x) { -sum(x) })
				}
				
			}
			
		}
		
	}
	
	colIsNA = apply(fv, 2, function(x) {any(is.na(x))})
	if (any(colIsNA)) {
		fv = fillOutFullValues_sumsToZero(fv, uniqueFL)
	}
	
	fv
}



# Unwanted function
# Assumes that partialWeights has proper names
# contrastType may only be a scalar character
fillInWeights = function(partialWeights, factors, fNames, contrastType, uniqueFL = NULL) {
	
	if (is.null(uniqueFL)) {
		uniqueFL = unique( subset(factors, select = fNames) )
	}
	
	# This does not assume fully crossed, but if the design is fully, crossed, it should work.
	fullWeights = matrix(NA, nrow=nrow(partialWeights), ncol=nrow(uniqueFL))
	
	# Set column names for fullWeights
	# This isn't really right, because the names of partialWeights may not be cell names
	allColNames = rep("", ncol(fullWeights))
	for (i in 1:nrow(uniqueFL)) {
		temp = subset(uniqueFL, subset = (i == 1:nrow(uniqueFL)) )
		allColNames[i] = makeCellName(as.list(temp))
	}
	colnames(fullWeights) = allColNames
	
	for (colname in colnames(partialWeights)) {
		fullWeights[ , colname ] = partialWeights[ , colname ]
	}
	
	if (contrastType == "contr.sum") {
		
		if (!isDesignFullyCrossed(factors, warnOnDuplicate=FALSE)) {
			stop("Design is not fully crossed: You can't use sums-to-zero contrasts for this purpose. Use contr.treatment instead.")
		}
		
		fullWeights = fillOutFullValues_sumsToZero(fullWeights, uniqueFL)
	} else if (contrastType == "contr.treatment" || contrastType == "contr.SAS") {
		fullWeights[ is.na(fullWeights) ] = 0
	} else {
		stop("Unsupported contrastType.")
	}
	
	fullWeights
	
}

# Unwanted function
# Gets part of S that has been filled in with information about the implicit effect parameters.
# contrastType may only be a scalar character
getPartialFilledS = function(factors, testedFactors, dmFactors = testedFactors, contrastType = NULL, warnOnDrop = FALSE) {
	
	################################################
	# This section is a C/P from getEffectParameters
	
	if (is.null(contrastType)) {
		if (isDesignFullyCrossed(factors, warnOnDuplicate=FALSE)) {
			contrastType = "contr.sum"
		} else {
			contrastType = "contr.treatment"
		}
	}
	
	# 1. Calculate the design matrix, X.
	dm = makeDesignMatrix(factors, dmFactors, contrastType, renameCols=FALSE)
	
	# 1b. If design is not fully crossed, strip excess terms from X.
	strippedInfo = stripExcessTermsFromDM(dm$mat)
	
	if (warnOnDrop && length(strippedInfo$dropped) > 0) {
		warning( paste0("The design is not fully crossed (unbalanced). As a result, some terms were dropped from the design matrix: ", paste(strippedInfo$dropped, collapse=", "), ". The following terms were kept: ", paste(strippedInfo$kept, collapse=", "), ". The naming of these terms depends on the contrastType that was used.") )
	}
	
	X = dm$mat = strippedInfo$mat
	
	
	# 2. Calculate S = (X'X)^-1 X'
	# If mu = X beta
	# then beta = (X'X)^-1 X'mu = S mu
	S = solve(t(X) %*% X) %*% t(X)
	
	# End C/P
	###############################
	
	
	S = t(S)
	
	mCols = getEffectAssignmentColumns(dm, testedFactors)
	
	S_s = subset(S, select=mCols)
	
	#S_s must have proper names before being passed to fillInWeights
	colnames(S_s) = renameDesignMatrixColumns(colnames(S_s), factors, testedFactors, contrastType = contrastType)
	
	full_s = fillInWeights(S_s, factors, testedFactors, contrastType = contrastType)
	
	full_s
}
