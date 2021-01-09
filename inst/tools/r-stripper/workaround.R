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

if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")) {
  writeLines(gsub("@ISYSTEM@", "I",
                  gsub("@CXX14STD@", "CXX14STD = -std=c++1y",
                       suppressWarnings(readLines("src/Makevars.in")))),
             "src/Makevars.win")
} else {
  writeLines(gsub("@ISYSTEM@", "isystem",
                  gsub("@CXX14STD@", "CXX14STD = -std=gnu++14",
                       suppressWarnings(readLines("src/Makevars.in")))),
             "src/Makevars")
}


unlink("R/RxODE_md5.R")

cpp <- list.files("src", pattern = ".(c|h|cpp|f)$")
include <- list.files("inst/include")
Rfiles <- list.files("R/", pattern = ".R")
md5 <- digest::digest(lapply(c(paste0("src/", cpp),
                               paste0("inst/include/", include),
                               paste0("R/", Rfiles)), digest::digest, file = TRUE))
unlink("R/RxODE_md5.R")
md5file <- file("R/RxODE_md5.R", "wb")
writeLines(sprintf("RxODE.md5 <- \"%s\"\n", md5), md5file)
close(md5file)

l <- readLines("DESCRIPTION")
w <- which(regexpr("Version[:] *(.*)$", l) != -1)
v <- gsub("Version[:] *(.*)$", "\\1", l[w])

unlink("src/ode.h")
ode.h <- file("src/ode.h", "wb")
writeLines(c(sprintf("#define __VER_md5__ \"%s\"", md5),
             "#define __VER_repo__ \"https://github.com/nlmixrdevelopment/RxODE\"",
             sprintf("#define __VER_ver__ \"%s\"", v)),
           ode.h)
close(ode.h)
