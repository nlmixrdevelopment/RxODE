## This is only for RxODE
for (f in c("inst/include/RxODE_RcppExports.h", "src/RcppExports.cpp")) {
  l <- readLines(f)
  w <- which(regexpr("^[#]include <RcppArmadillo.h>", l) != -1)
  if (length(w) == 1) {
    l <- l[-w]
    message("Excluding RcppArmadillo from", f)
    unlink(f)
    file.out <- file(f, "wb")
    writeLines(l, file.out)
    close(file.out)
  }
}

.in <- suppressWarnings(readLines("src/Makevars.in"))
.in <- gsub("@ARMA@", file.path(find.package("RcppArmadillo"),"include"), .in)
.in <- gsub("@BH@", file.path(find.package("BH"),"include"), .in)
.in <- gsub("@RCPP@", file.path(find.package("Rcpp"),"include"), .in)
.in <- gsub("@EG@", file.path(find.package("RcppEigen"),"include"), .in)

.badStan <- ""
if (R.version$major < 4 && isTRUE(.Platform$OS.type == "windows")) {
  # See https://github.com/nlmixrdevelopment/nlmixr/issues/429#issuecomment-729993005
  if (compareVersion(as.character(packageVersion("BH")), "1.66.0-1") > 0) {
    .badStan <- " -D__Rx_noSTAN__ "
  }
}
.in <- gsub("@SH@", gsub("-I", "-@ISYSTEM@",
                         paste(## capture.output(StanHeaders:::CxxFlags()),
                               ## capture.output(RcppParallel:::CxxFlags()),
                               paste0("-@ISYSTEM@'", system.file('include', 'src', package = 'StanHeaders', mustWork = TRUE), "'"),
                               .badStan)),
            .in)

.in <- gsub("@SL@", "", ##paste(capture.output(StanHeaders:::LdFlags()), capture.output(RcppParallel:::RcppParallelLibs())),
            .in)

if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")) {
  .in <- gsub("@CXX14STD@", "-std=c++1y", .in)
  file.out <- file("src/Makevars.win", "wb")
  writeLines(gsub("@ISYSTEM@", "I", .in),
             file.out)
  close(file.out)
} else {
  .in <- gsub("@CXX14STD@", "-std=gnu++14", .in)
  file.out <- file("src/Makevars", "wb")
  writeLines(gsub("@ISYSTEM@", "isystem", .in),
             file.out)
  close(file.out)
}

if (file.exists("src/Makevars.in.r-stripper.bak")) {
  ## Reset to old Makevars.in
  l <- readLines("src/Makevars.in.r-stripper.bak")
  file.out <- file("src/Makevars.in", "wb")
  writeLines(l, file.out)
  close(file.out)
  unlink("src/Makevars.in.r-stripper.bak")
}

if (file.exists("man/reexports.Rd")) {
  l <- readLines("man/reexports.Rd")
  if (!any(regexpr("[\\]value", l) != -1)) {
    l <- c(l, "\\value{ Inherited from parent routine }")
    file.out <- file("man/reexports.Rd", "wb")
    writeLines(l, file.out)
    close(file.out)
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

def <- unique(c(def, c("_sum", "_sign", "_prod", "_max", "_min", "_transit4P", "_transit3P", "_assignFuns0", "_assignFuns", "_getRxSolve_", "_solveData")))

## deparse1 came from R 4.0, use deparse2
deparse2 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
  paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}

final <- c("#include <time.h>",
           "#include <stdlib.h>",
           "void writeHeader(const char *md5, const char *extra) {",
           "unsigned long int timeId=0;",
           paste0("sAppend(&sbOut, \"#define ", def, " _rx%s%s%ld\\n\", extra, md5, timeId++);"),
           "}",
           "void writeBody() {",
           paste0("sAppendN(&sbOut, ", vapply(paste0(l, "\n"), deparse2, character(1)), ", ", nchar(l) + 1, ");"),
           "}",
           "void writeFooter() {",
           paste0("sAppendN(&sbOut, \"#undef ", def, "\\n\", ", nchar(def) + 8, ");"),
           "}"
           )

codegen2.h <- file("src/codegen2.h", "wb")
writeLines(final,
           codegen2.h)
close(codegen2.h)
