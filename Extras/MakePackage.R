library(devtools)
library(roxygen2)

#install.packages("R.rsp")



#increment pkgVer!!!
pkgVer = "0.1.3"

setwd("~/../Programming/R/CMBBHT/")

devtools::use_package("polspline")

file.remove("NAMESPACE")
devtools::document()

desc = read.dcf("DESCRIPTION")
desc[,"Version"] = pkgVer
desc[,"Date"] = format(Sys.Date(), "%Y")
write.dcf(desc, "DESCRIPTION")

lic = read.dcf("LICENSE")
lic[,"YEAR"] = format(Sys.Date(), "%Y")
write.dcf(lic, "LICENSE")

devtools::install(build_vignettes = TRUE)







devtools::check()





remove.packages("CMBBHT")

devtools::install_github("hardmanko/CMBBHT", build_vignettes = TRUE)

devtools::install_github("hardmanko/CMBBHT@v0.1.3", build_vignettes = TRUE)


vignette("Manual", "CMBBHT")
vignette("BinomialTutorial", "CMBBHT")



#devtools::build(path="packaged")

# File name depends on the version number.
install.packages( paste0("packaged/CMBBHT_", pkgVer, ".tar.gz"), repos=NULL)
