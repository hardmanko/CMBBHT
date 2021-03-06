% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{testIntercept}
\alias{testIntercept}
\title{Test the Value of the Intercept}
\usage{
testIntercept(priorCMs, postCMs, factors, testVal,
  dmFactors = stats::formula(" ~ 1"), contrastType = NULL,
  testFunction = valueTest_SDDR)
}
\arguments{
\item{priorCMs}{See \code{\link{testHypothesis}}.}

\item{postCMs}{See \code{\link{testHypothesis}}.}

\item{factors}{See \code{\link{testHypothesis}}.}

\item{testVal}{The value that will be tested. Depends on the choice of \code{testFunction}. See \code{\link{valueTest_SDDR}} for the interpretation of this argument by the default test function.}

\item{dmFactors}{See \code{\link{testHypothesis}}.}

\item{contrastType}{See \code{\link{testHypothesis}}.}

\item{testFunction}{A function of 3 positional arguments: the prior intercept, the posterior intercept, and the point at which the hypothesis is tested (i.e. \code{testVal} is passed as the third argument). A custom test function may choose to ignore \code{testVal} if the tested values are determined in some other way. The prior and posterior intercepts will be vectors or one-column matrices.}
}
\value{
The result of the function depends on the choice of \code{testFunction}. See \code{\link{valueTest_SDDR}} for an example.
}
\description{
Tests whether the intercept/grand mean has some given value. For a general test of whether a parameter has a given value, see \code{\link{valueTest_SDDR}}. Note that in order for this test to be very meaningful (probably)
}
