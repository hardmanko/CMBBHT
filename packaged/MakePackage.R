library(devtools)
library(roxygen2)

#install.packages("R.rsp")

setwd("~/../Programming/R/CMBBHT/")

devtools::use_package("polspline")

file.remove("NAMESPACE")
devtools::document()

devtools::install()




devtools::check()



devtools::build(path="packaged")

# File name depends on the version number.
install.packages("packaged/CMBBHT_0.1.tar.gz", repos=NULL)

