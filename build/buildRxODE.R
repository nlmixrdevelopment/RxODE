setwd(devtools::package_file("build"))
source("../R/RxODE-win-setup.R")
rxPhysicalDrives <- function(...){"C:\\"}
.normalizePath <- normalizePath
rxWinSetup()
devtools::install()
