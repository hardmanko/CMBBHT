CMBBHT: Cell-Means-Based Bayesian Hypothesis Tests
=================

This package performs Bayesian hypothesis tests of designs where cell means are estimated with single parameters. It also works for designs where it is not cell means themselves that are estimated, but differences from a grand mean or intercept. It includes functions for testing main effects and interactions. It can used with nested and unbalanced designs, with a little effort.

Installation
------------

This package can be installed directly from this repository with the following commands in R:
```{r}
#You can skip this if you already have these packages installed and up-to-date
install.packages(c("devtools", "R.rsp"))

#Install the package
devtools::install_github("hardmanko/CMBBHT", build_vignettes = TRUE)
```

Once the package is installed, open the manual with
```{r}
vignette("Manual", "CMBBHT")
```
You can also find a pdf of the manual in the `vignettes` directory of this repository.
