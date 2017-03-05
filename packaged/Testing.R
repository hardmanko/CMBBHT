
library(CMBBHT)

?getEffectParameters
?testHypothesis


factors = data.frame(let = rep(letters[1:5], each=3),
										 num = factor(rep(1:3, 5)),
										 stringsAsFactors = FALSE)


factors = data.frame(let = c('a', 'a', 'b', 'b', 'c', 'c'),
										 num = c('1', '2', '1', '2', '1', '2'),
										 stringsAsFactors = FALSE)

mus = c( 0,0,   2,2,   4,4 )
mus = c( 1,0,   3,2,   3,4 )

factors = data.frame(let = rep(letters[1:3], each=3), 
										 num = factor(rep(1:3, 3)), 
										 stringsAsFactors = FALSE)

#not fully crossed (above minus something)
factors = factors[ -c(6,8,9), ]

mus = c( 0, 2, 4,  0, 2,  0 )



factors = data.frame(low = rep(c('a', 'b', 'c'), each=4),
										 num = rep(c('1', '2'), 3, each=2),
										 up = rep(c('X', 'Y'), 6),
										 stringsAsFactors = FALSE)

#above minus something
factors = factors[ c(-1,-5,-6), ]



#bad factor stuff
factors$let = factor(factors$let)
levels(factors$let) = c( levels(factors$let), 'd', 'e')
factors$num = as.numeric(factors$num)




#posteriorSD tunes how much data you have => more precise posterior.
getPriorAndPosteriorCMs = function(mus, posteriorSD, priorLocation, priorScale, samples=10000) {

	priorCMs = postCMs = matrix(nrow = samples, ncol=length(mus))
	for (i in 1:ncol(postCMs)) {
		priorCMs[,i] = rcauchy(nrow(priorCMs), priorLocation, priorScale)
		postCMs[,i] = rnorm(nrow(postCMs), mus[i], posteriorSD)
	}
	
	list(prior=priorCMs, post=postCMs)
	
}




citf = function(prior, posterior) {
	testFunction_credibleInterval(prior, posterior, 0.99)
}

posttf = function(prior, posterior) {
	hist(posterior)
}



cms = getPriorAndPosteriorCMs(mus, 0.5, 0, 4)

testedFactors = c("let")
dmFactors = formula(" ~ let * num")
contrType = "contr.sum"


testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, contrastType = contrType,
							 testFunction = testFunction_sd_diffuse)

testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, contrastType = contrType,
							 testFunction = testFunction_savageDickey_diffusePrior)

testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, contrastType = contrType,
							 testFunction = testFunction_lowerTail)

testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, contrastType = contrType,
							 testFunction = citf)

pd = testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, contrastType = contrType,
							 testFunction = posttf)


priorEffects = getEffectParameters(cms$prior, factors, testedFactors, dmFactors, contrastType = contrType)
postEffects = getEffectParameters(cms$post, factors, testedFactors, dmFactors, contrastType = contrType)

summary = summarizeEffectParameters(postEffects)
plotEffectParameterSummary(summary)

groupEffectParameters(postEffects, method="credInt", CIP = 0.95)

res = groupEffectParameters(postEffects, method="BayesFactor", priorEffects = priorEffects, equalBF = 3)

res$bf01



d = apply(priorEffects, 1, var)


pm = getEffectParameters(cms$post, factors, testedFactors = "(Intercept)", dmFactors="~ let", contrastType = "contr.treatment")
hist(pm)



tf = function(prior, posterior) {
	savageDickey_valTest(prior, posterior, 0)
}

testHypothesis(cms$prior, cms$post, factors, "(Intercept)", contrastType = "contr.treatment", testFunction = tf)

priorCMs = cms$prior
postCMs = cms$post





testIntercept(cms$prior, cms$post, factors, 0)

summary = summarizeEffectParameters(postEffects, fList = list(mean=mean, median=median, sd=sd))
summary$lower = summary$median - summary$sd
summary$upper = summary$median + summary$sd
plotEffectParameterSummary(summary, ctName = "median", varNames=c("lower", "upper"))


#The first three are orthogonal (for balanced designs); the last two are not orthogonal.
#For balanced designs, you should get the same results from the first three. (Actually, for all designs.)
#For all designs, you should get the same results from the last two.
contrastTypes = c("contr.sum", "contr.helmert", "contr.poly", "contr.treatment", "contr.SAS")

res = NULL

for (ct in contrastTypes) {
	
	ht = testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, contrastType = ct)
	
	temp = data.frame(ct = ct, 
										testedFactors = paste(testedFactors, collapse=","), 
										dmFactors = paste(dmFactors, collapse=" "),
										bf10 = ht$bf10, bf01 = ht$bf01)
	res = rbind(res, temp)
	
}
res



#######################
# Try different test functions

testedFactors = c("let")
dmFactors = c("let", "num")

tf1 = function(prior, posterior) {
	testFunction_sd_diffuse(prior, posterior, postMaxMult = 10)
}

tf2 = function(prior, posterior) {
	testFunction_sd_diffuse(prior, posterior, postMaxMult = 10, truncatePosterior = FALSE)
}

tf3 = function(prior, posterior) {
	testFunction_sd_diffuse(prior, posterior, postMaxMult = 10000)
}

testFuns = list(tf1=tf1, tf2=tf2)
res = NULL

for (i in 1:50) {
	for (tfn in names(testFuns)) {
		
		cms = getPriorAndPosteriorCMs(mus, 0.5, 0, 4)
		
		tf = testFuns[[ tfn ]]
		
		ht = testHypothesis(cms$prior, cms$post, factors, testedFactors, dmFactors=dmFactors, testFunction = tf)
		
		temp = data.frame(tfn = tfn, i = i,
											testedFactors = paste(testedFactors, collapse=","), 
											dmFactors = paste(dmFactors, collapse=" "),
											bf10 = ht$bf10, bf01 = ht$bf01, pKept = ht$pKept)
		res = rbind(res, temp)
		
	}
}
#res

aggregate(bf10 ~ tfn, res, mean)
aggregate(bf10 ~ tfn, res, median)

aggregate(pKept ~ tfn, res, mean)


# Testing truncated normal density estimate.
# How to pick bw?

library(msm)

makeSmoothedDensity = local({
	tndens = function(value, bw) {
		force(value); force(bw)
		function(z) dtnorm(z, mean = value, sd = bw, lower=0)
	}
	function(x, bw) {
		
		
		flist = lapply(x, tndens, bw = bw)
		
		f = function(z) {
			
			if(length(z) <= 1) {
				dens = mean(sapply(flist, function(fun) fun(z)))
			}	else {
				dens = rowMeans(sapply(flist, function(fun) fun(z)))
			}
			dens
		}
		
		f
	}
})

calcDens = function(new, old, bw) {
	
	
	
}

crossValidate = function(x, densFunMaker, startingBW) {
	
	allBW = startingBW ^ seq(-6, 6, length.out = 20)
	
	meanDens = NULL
	
	for (bw in allBW) {
		
		samp = sample(1:length(x), floor(length(x) / 2))
		
		h1 = x[ samp ]
		h2 = x[ -samp ]
		
		h1f = densFunMaker(h1, bw)
		h2f = densFunMaker(h2, bw)
		
		allDens = c(h1f(h2), h2f(h1))
		meanDens = c(meanDens, mean(allDens))
		
	}
	
	plot(meanDens)
	
}


bw = 5

prior_f = makeSmoothedDensity(prior_D, bw)
post_f = makeSmoothedDensity(post_D, bw)

prior_f(0)
post_f(0)

post_f(0) / prior_f(0)


xs = seq(0, 50, 0.2)
plot(xs, prior_f(xs), type='l')

hist(prior_D[ prior_D < 50 ], prob=TRUE)
lines(xs, prior_f(xs))

dens = prior_f(prior_D)
mean(dens)




library(polspline)

g = rgamma(1000, 2, 2)
hist(g)

ls1 = logspline(g)
ls2 = logspline(g, lbound=0)
ls3 = logspline(g, lbound=0, ubound=max(g))
ls4 = logspline(g, lbound=0, ubound=max(g) * 10)

lsl = list(ls1, ls2, ls3, ls4)

dens = NULL
for (ls in 1:length(lsl)) {
	d = dlogspline(0, lsl[[ls]])
	dens = c(dens, d)
}
dens / max(dens)

logspline(g, lbound=NULL)





de = d[ d < 100 ]

gamfun = function(par, d) {
	
	alpha = max(par[1], 0.001)
	beta = max(par[2], 0.001)
	
	ld = dgamma(d, alpha, beta, log=TRUE)
	-sum(ld)
}

opt = optim(c(1,1), gamfun, d=d)

expfun = function(par, d) {
	ld = dexp(d, par, log=TRUE)
	-sum(ld)
}

opt = optimize(expfun, interval=c(0.001, 50), d=de)

alpha = max(opt$par[1], 0.001)
beta = max(opt$par[2], 0.001)

hist(d[ d < 100 ], breaks=100, prob=TRUE)
xs = seq(0, 100, 1)
lines(xs, dgamma(xs, alpha, beta))
lines(xs, dexp(xs, 0.035))



ls2 = logspline(de, lbound=0)
dlogspline(0, ls2)



contained = function(x, limits) {
	all(limits[1] < x & x < limits[2])
}

rangeI = function(eff) {
	contained(eff, limits=c(-1, 1))
}

#Based on Wetzels Et Al 2010
goodIntervalTest = function(prior_eff, post_eff, I_fun) {
	
	prior_I = apply(prior_eff, 1, I_fun)
	post_I = apply(post_eff, 1, I_fun)
	
	bf_01 = mean(post_I) / mean(prior_I) #In favor of the constrained hypothesis/model
	bf_10 = 1 / bf_01 # In favor of the general hypothesis/model
	
	list(bf_10 = bf_10, bf_01 = bf_01)
}

goodIntervalTest(prior_eff, post_eff, rangeI)



