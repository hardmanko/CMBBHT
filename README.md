CMBBHT: Cell-Means-Based Bayesian Hypothesis Tests
=================

This package performs Bayesian hypothesis tests of designs where cell means are estimated with single parameters. It also works for designs where it is not cell means themselves that are estimated, but differences from a grand mean or intercept. It includes functions for testing main effects and interactions. It can used with nested and unbalanced designs, with a little effort.

Installation
------------

This package can be installed directly from this repository with the following commands in R:
```{r}
#You don't need to do this if you already have devtools installed.
install.packages("devtools")

#Install the package, including the manual
devtools::install_github("hardmanko/CMBBHT", build_vignettes = TRUE)
```

Once the package is installed, open the manual with
```{r}
vignette("Manual", "CMBBHT")
```
