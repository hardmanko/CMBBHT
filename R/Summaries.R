
#' Summary Statistics of Effect Parameters
#' 
#' Calculates summary statistics, like mean, median, and credible interval, of effect parameters. 
#' 
#' Note that equality/inequality of effect parameters should not be determined based on their individual credible intervals. You should instead use the credible interval of their difference, which is done by [`groupEffectParameters`] by default.
#' 
#' @param effects A matrix, such as from [`getEffectParameters`], where each column is one parameter and each row is an iteration. If the columns are named with the name of the effect, those names will be used.
#' @param fList A named list of functions of one vector argument that will be applied to the effects individually.
#' @param ps Proportions for which quantiles will be calculated. By default, the median and 95\% credible interval is calculated.
#' 
#' @return A `data.frame` with summary statistics for each column in `effects`.
#' 
#' @seealso For plotting the result of this function, [`plotEffectParameterSummary`]. For grouping indistinguishable effects, [`groupEffectParameters`].
#' 
#' @md
#' @export
#' 
summarizeEffectParameters = function(effects, fList = list(mean=mean), ps = c(0.025, 0.5, 0.975)) {
	
	if (is.null(colnames(effects))) {
		colnames(effects) = paste("effect", 1:ncol(effects), sep="_")
	}
	
	res = NULL
	for (i in 1:ncol(effects)) {
		
		ef = effects[,i]
		
		temp = data.frame(effect = colnames(effects)[i])
		
		for (n in names(fList)) {
			temp[ , n ] = fList[[n]](ef)
		}
		
		qs = stats::quantile(ef, ps )
		for (j in 1:length(qs)) {
			temp[,names(qs)[j]] = qs[j]
		}
		
		res = rbind(res, temp)
		
	}
	
	res
	
}

#' Plot Effects Parameter Summary
#' 
#' Takes the result of [`summarizeEffectParameters`] and plots a measure of central tendency (default mean) and variability (default credible interval).
#' 
#' @param summary A `data.frame` such as that returned by [`summarizeEffectParameters`].
#' @param ctName The name of the column in `summary` that contains the measure of central tendency to plot.
#' @param varNames A length 2 vector of names of the columns of `summary` that contain measures of variability in the form of endpoints of an error bar. If `NULL`, it selects the lowest and highest percentiles present in `summary`.
#' @param mar Passed to the `mar` argument of `par()`. If `NULL`, it is ignored.
#' 
#' @md
#' @export
plotEffectParameterSummary = function(summary, ctName = "mean", varNames = NULL, mar=c(3,6,1,1)) {
	
	if (is.null(varNames)) {
		n = names(summary)
		whichPercent = which(grepl("%", n, fixed=TRUE))
		n = n[ whichPercent ]
		n = gsub("%", "", n, fixed=TRUE)
		n = as.numeric(n)
		
		cimin = names(summary)[ whichPercent[ which.min(n) ] ]
		cimax = names(summary)[ whichPercent[ which.max(n) ] ]
		varNames = c(cimin, cimax)
	}
	
	xlim = c(min(summary[ , cimin ]), max(summary[ , cimax ]))
	ylim = c(0.5, nrow(summary) + 0.5)
	
	#Invert so that effects are listed top to bottom
	summary = summary[ seq.int(nrow(summary), 1), ]
	
	if (!is.null(mar)) {
		graphics::par(mar=mar)
	}
	
	graphics::plot(x=0, y=0, xlim=xlim, ylim=ylim, type='n', axes=FALSE, ylab="", xlab="")
	graphics::box()
	graphics::axis(1)
	graphics::axis(2, at=1:nrow(summary), labels=summary$effect, las=2)
	for (i in 1:nrow(summary)) {
		
		graphics::points(summary[ i, ctName ], i, pch=16)
		
		lower = summary[ i, varNames[1] ]
		upper = summary[ i, varNames[2] ]
		
		if (lower != upper) {
			graphics::lines(c(lower, upper), c(i,i))
			graphics::lines(c(lower, lower), i + c(-0.2, 0.2))
			graphics::lines(c(upper, upper), i + c(-0.2, 0.2))
		}
		
	}
	
}

#' Group Indistinguishable Effect Parameters
#' 
#' Perform pairwise comparisons of effect levels to determine which effects are indistinguishable from one another. By default, this is based on credible intervals. For people who do not want to, or cannot, use Bayes factors estimated with the Savage-Dickey density ratio, this function allows for another method of determining whether effects are present in the data. In addition, this allows for an analysis of which factor levels differ from one another.
#' 
#' @param postEffects A matrix of posterior effect parameters, as from [`getEffectParameters`].
#' @param method The method to use. `"credInt"`: Use credible intervals of the difference between effects to determine equality. `"BayesFactor"`: Use Bayes factors of a hypothesis test that the difference between effects is 0.
#' @param CIP Used if `method == "credInt"`. Credible interval proportion.
#' @param priorEffects Used if `method == "BayesFactor"`. A matrix of prior effects with columns ordered the same as \code{postEffects}.
#' @param equalBF Used if `method == "BayesFactor"`. If the Bayes factor in favor of the null is greater than equalBF, the effects are considered to be equivalent.
#' 
#' @return Invisibly, a list with the following elements:
#' * `grp`: A matrix containing grouping information. All effects that share a letter are Indistinguishable. This is what is printed.
#' * `eqm`: A boolean matrix contining information about which effects were equal (indistinguishable). The value `TRUE` indicates equality.
#' * `bf01`: If `method == "BayesFactor"`, a matrix containing Bayes factors in favor of the hypothesis that the two effects were equal. 
#' 
#' @md
#' @export
groupEffectParameters = function(postEffects, method="credInt", CIP = 0.95, priorEffects = NULL, equalBF = 3) {
	
	if (method == "credInt") {
		if (is.null(CIP)) {
			stop("CIP must be provided if method == \"credInt\".")
		}
	} else if (method == "BayesFactor") {
		if (is.null(priorEffects)) {
			stop("priorEffects must be provided if method == \"BayesFactor\".")
		}
		if (is.null(equalBF)) {
			stop("equalBF must be provided if method == \"BayesFactor\".")
		}
	} else {
		stop("Invalid method. The valid methods are \"credInt\" and \"BayesFactor\".")
	}
	
	
	cips = c((1 - CIP) / 2, (1 + CIP) / 2)
	
	
	ord = order( apply(postEffects, 2, mean) )
	postEffects = postEffects[ , ord ]
	if (!is.null(priorEffects)) {
		priorEffects = priorEffects[ , ord ]
	}
	
	
	n = ncol(postEffects)
	
	eqm = matrix(NA, nrow=n, ncol=n)
	diag(eqm) = TRUE
	rownames(eqm) = colnames(eqm) = colnames(postEffects)
	
	bfm = eqm
	diag(bfm) = Inf
	
	for (i in 1:n) {
		for (j in 1:n) {
			if (j <= i) {
				#For the lower tri, copy results from the upper tri.
				eqm[i,j] = eqm[j,i]
				bfm[i,j] = bfm[j,i]
			} else {
				
				if (method == "credInt") {
					
					difci = stats::quantile(postEffects[,i] - postEffects[,j], cips)
					eqm[i,j] = difci[1] <= 0 && difci[2] >= 0
					
				} else if (method == "BayesFactor") {
					
					postDif = postEffects[,i] - postEffects[,j]
					priorDif = priorEffects[,i] - priorEffects[,j]
					
					
					sdr = valueTest_SDDR(priorDif, postDif, testVal=0)
					eqm[i,j] = sdr$bf01 >= equalBF
					bfm[i,j] = sdr$bf01
					
				}
				#TODO: Maybe allow user to pass test function. It's not clear how to do this, 
				#given that priorEffects are not even passed by default. 
				#NULL[,i] == NULL, so maybe that is the solution?
				#else if (method == "userFun") {
				#	eqm[i,j] = userFun(priorEffects[,i], priorEffects[,j], postEffects[,i], postEffects[,j])
				#}
			}
			
		}
	}
	
	
	foundShared = NULL
	groupMat = NULL
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			if (eqm[i,j]) {
				
				shared = eqm[i,] & eqm[j,]
				
				
				sharedStr = paste(which(shared), collapse=",")
				
				if (!(sharedStr %in% foundShared)) {
					foundShared = c(foundShared, sharedStr)
					groupMat = rbind(groupMat, shared)
				}
				
			}
		}
	}
	
	#second pass: add groups for cells that are in no other group
	for (i in 1:n) {
		notInOtherGroup = !any(groupMat[,i])
		if (notInOtherGroup) {
			
			shared = i == 1:n
			
			groupMat = rbind(groupMat, shared)
			foundShared = c(foundShared, paste(which(shared), collapse=","))
		}
	}
	
	
	grpm = matrix("", nrow=n, ncol=length(foundShared))
	rownames(grpm) = colnames(postEffects)
	for (i in 1:length(foundShared)) {
		
		wh = as.numeric(strsplit(foundShared[i], split=",", fixed=TRUE)[[1]])
		
		grpm[wh,i] = LETTERS[i]
	}
	
	colnames(grpm) = rep("", ncol(grpm))
	
	print(grpm, quote=FALSE)
	
	res = list(grp = grpm, eq = eqm)
	if (method == "BayesFactor") {
		res$bf01 = bfm
	}
	
	invisible(res)
	
}
