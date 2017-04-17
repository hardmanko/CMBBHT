
#' Determine Whether Factorial Design is Fully Crossed
#' 
#' Determine whether the factorial design in `factors` is a fully-crossed design.
#' 
#' @param factors See [`testHypothesis`].
#' @param warnOnDuplicate If `TRUE` and there are duplicate rows in `factors`, a warning that extent is emitted. Note that there are legitmate reasons to have duplicate rows in `factors`.
#' 
#' @return `TRUE` if the design is fully crossed, `FALSE` otherwise.
#' 
#' @md
#' @export
isDesignFullyCrossed = function(factors, warnOnDuplicate = TRUE) {
	
	totalCells = 1
	for (n in names(factors)) {
		totalCells = totalCells * length(unique(factors[,n]))
	}
	
	uf = unique(factors)
	if (warnOnDuplicate && nrow(uf) != nrow(factors)) {
		warning("There are duplicate rows in factors. This may be ok.")
	}
	
	factors = uf
	if (nrow(factors) == totalCells) {
		return(TRUE)
	}
	
	FALSE
}


# Internal function
# cnl can be a list or data.frame.
makeCellName = function(cnl) {
	n = ""
	fNames = names(cnl)
	for (i in 1:length(fNames)) {
		fact = fNames[i]
		lev = cnl[[fact]]
		n = paste(n, fact, ".", lev, sep="")
		if (i < length(fNames)) {
			n = paste(n, ":", sep="")
		}
	}
	n
}

# Internal function. Goes with makeCellName.
splitCellName = function(str) {
	parts = strsplit(str, split = ":", fixed=TRUE)[[1]]
	
	fl = list()
	
	for (i in 1:length(parts)) {
		
		parts2 = strsplit(parts[i], split=".", fixed=TRUE)[[1]]
		
		fact = parts2[1]
		lev = parts2[2]
		
		fl[[ fact ]] = lev
		
	}
	fl
}


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
		stop("Unsupported contrastType. The supported types are contr.sum, contr.treatment, and contr.SAS")
	}
	
	fullWeights
	
}

#' Get Part of S, Filled with Implicit Parameters
#' 
#' \code{S = solve(t(X) \%*\% X) \%*\% t(X)}.
#' Gets part of S, selected with `testedFactors`, that has been filled in with information about the implicit effect parameters.
#' 
#' 
#' @param factors See [`testHypothesis`].
#' @param testedFactors [`testHypothesis`].
#' @param dmFactors [`testHypothesis`].
#' @param contrastType May only be one of `"contr.sum"`, `"contr.treatment"`, or `"contr.SAS"`.
#' 
#' @return The selected part of the `S` matrix, with additional, implicit parameters included.
#' 
#' @md
#' @export
getPartialFilledS = function(factors, testedFactors, dmFactors, contrastType) {
	
	if (!is.character(contrastType)) {
		stop("For this function, contrastType may only be a string.")
	}
	
	testedFactors = splitFactorNames(testedFactors, convertToFormula = FALSE)
	dmFactors = splitFactorNames(dmFactors)
	
	dm = makeDesignMatrix(factors, dmFactors, contrastType, renameCols=FALSE)
	
	strippedInfo = stripExcessTermsFromDM(dm)
	X = dm$mat = strippedInfo$mat

	S = solve(t(X) %*% X) %*% t(X)

	# Transpose to make each column related to an effect
	S = t(S)
	
	mCols = getEffectAssignmentColumns(dm, testedFactors)
	
	S_s = subset(S, select=mCols)
	
	#S_s must have proper names before being passed to fillInWeights
	colnames(S_s) = renameDesignMatrixColumns(colnames(S_s), factors, testedFactors, contrastType = contrastType)
	
	fillInWeights(S_s, factors, testedFactors, contrastType = contrastType)
}


