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

if (file.exists("src/Makevars.in.r-stripper.bak")) {
  l <- readLines("src/Makevars.in.r-stripper.bak")
  writeLines(l, "src/Makevars.in")
  unlink("src/Makevars.in.r-stripper.bak")
}

if (file.exists("man/reexports.Rd")) {
  l <- readLines("man/reexports.Rd")
  if (!any(regexpr("[\\]value", l) != -1)) {
    l <- c(l, "\\value{ Inherited from parent routine }")
    writeLines(l, "man/reexports.Rd")
  }
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

unlink("src/codegen2.h")
l <- readLines("inst/include/RxODE_model_shared.c")

l <- l[l != ""]
l <- gsub(" *= *NULL;", "=NULL;", l)

def <- l
w <- which(regexpr("double _prod", def) != -1) - 1
def <- def[1:w]
def <- gsub("=NULL", "", def)
def <- gsub("[^ ]* *[*]?([^;]*);", "\\1", def)

def <- c(def, c("_sum", "_sign", "_prod", "_max", "_min", "_transit4P", "_transit3P", "_assignFuns0", "_assignFuns", "_getRxSolve_"))

## deparse1 came from R 4.0, use deparse2
deparse2 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
  paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

final <- c("#include <time.h>",
           "void writeHeader() {",
           "time_t timeId;",
           "timeId = time(NULL);",
           paste0("sAppend(&sbOut, \"#define ", def, " ", def, "%ld\\n\", timeId);"),
           "}",
           "void writeBody() {",
           paste0("sAppendN(&sbOut, ", vapply(paste0(l, "\n"), deparse2, character(1)), ", ", nchar(l) + 1, ");"),
           "}"
           )

codegen2.h <- file("src/codegen2.h", "wb")
writeLines(final,
           codegen2.h)
close(codegen2.h)
