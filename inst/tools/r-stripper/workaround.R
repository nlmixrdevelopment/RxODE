## This is only for RxODE
for (f in c("inst/include/RxODE_RcppExports.h", "src/RcppExports.cpp")) {
  l <- readLines(f)
  w <- which(regexpr("^[#]include <RcppArmadillo.h>", l) != -1)
  if (length(w) == 1) {
    l <- l[-w]
    message("Excluding RcppArmadillo from", f)
    writeLines(l, f)
  }
}

if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")){
  writeLines(gsub("@ISYSTEM@", "I",
                  gsub("@CXX14STD@", "CXX14STD = -std=c++1y",
                       suppressWarnings(readLines("src/Makevars.in")))),
             "src/Makevars.win")
}


unlink("R/RxODE_md5.R")

cpp <- list.files("src", pattern = ".(c|h|cpp|f)$")
include <- list.files("inst/include")
Rfiles <- list.files("R/", pattern = ".R")
md5 <- digest::digest(lapply(c(paste0("src/", cpp),
                               paste0("inst/include/", include),
                               paste0("R/", Rfiles)), digest::digest, file = TRUE))
writeLines(sprintf("RxODE.md5 <- \"%s\"\n", md5), "R/RxODE_md5.R")

l <- readLines("DESCRIPTION")
w <- which(regexpr("Version[:] *(.*)$", l) != -1)
v <- gsub("Version[:] *(.*)$", "\\1", l[w])

unlink("src/ode.h")
writeLines(c(sprintf("#define __VER_md5__ \"%s\"", md5),
             "#define __VER_repo__ \"https://github.com/nlmixrdevelopment/RxODE\"",
             sprintf("#define __VER_ver__ \"%s\"", v)),
           "src/ode.h")
