CMBBHT: Cell-Means-Based Bayesian Hypothesis Tests
=================

This package performs Bayesian hypothesis tests of designs where cell means are estimated with single parameters. It also works for designs where it is not cell means themselves that are estimated, but differences from a grand mean. It includes functions for testing main effects and interactions. It can used with nested and unbalanced designs, with a little effort.

Installation
------------

This package can be installed directly from this repository with the following commands in R:
```{r}
#You don't need to do this if you already have devtools installed.
install.packages("devtools")

#Install the package
devtools::install_github("hardmanko/CMBBHT")
```

Once the package is installed, open the manual with
```{r}
vignette("Manual", "CMBBHT")
```

It is slightly better to install the package from a package file, which you can get by downloading this repository and getting the file from the "packaged" subdirectory. Once you have that file on you computer somewhere, you can install the package using:
```{r}
install.packages(type = "source", repos = NULL, pkgs = "path_to_package_file/CMBBHT_0.1.tar.gz")
```

To just test the package without installing it, you can simply source `./R/CMBBHT-functions.R`, which also contains all of the functions and the documentation as comments.
