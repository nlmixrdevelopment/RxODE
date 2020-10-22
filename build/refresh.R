library(devtools)

if (!file.exists(devtools::package_file("src/liblsoda"))){
    owd <- getwd();
    setwd(devtools::package_file("src"));
    system("git clone git@github.com:nlmixrdevelopment/liblsoda")
    setwd(owd);
}

for (file in list.files(devtools::package_file("src/liblsoda/src"), pattern="[.]([hc]|inc)$")){
    if (file == "cfode_static.c"){
    } else{
        message(sprintf("\tcopy %s", file))
        lines <- suppressWarnings(readLines(file.path(devtools::package_file("src/liblsoda/src"), file)));
        lines <- gsub("#include \"cfode_static.inc\"", "#include \"cfode_static.h\"", lines, fixed=TRUE);
        lines <- gsub("fprintf[(] *stderr *, *", "REprintf(", lines);
        lines <- gsub("([ \t])printf[(] *", "\\1REprintf(", lines);
        if (any(regexpr("REprintf", lines) != -1)){
            lines <- c("#include <R.h>", "#include <Rinternals.h>", lines);
        }
        if (any(regexpr("#define ERROR", lines, fixed=TRUE) != -1)){
            lines <- gsub("#define ERROR", "#ifdef ERROR\n#undef ERROR\n#endif\n#define ERROR", lines, fixed=TRUE);
        }
        if (file == "solsy.c"){
            lines <- c("#include <R.h>", "#include <Rinternals.h>", gsub("abort[(] *[)]", "error(\"liblsoda does not implement this. (solsy)\")", lines));
        }
        if (file == "lsoda.c"){
            lines <- gsub("int[ \t]+i(hit)?[ \t]*([;,])", "int i\\1=0\\2", lines, perl=TRUE);
            lines <- gsub("i=0, *ihit;", "i=0, ihit=0;", lines);
        }
        writeLines(lines, file.path(devtools::package_file("src"), gsub("[.]inc$", ".h", file)))
    }
}


## library(tidyr)
## library(dplyr)
cat("Copy header to inst directory")

file.copy(devtools::package_file("src/RxODE_types.h"),
          devtools::package_file("inst/include/RxODE_types.h"),
          overwrite=TRUE);

ode.h();

cat("Update Parser c file\n");
dparser::mkdparse(devtools::package_file("inst/tran.g"),
         devtools::package_file("src/"),
         grammar_ident="RxODE");
file <- gsub("^([#]line [0-9]+ )\".*(src)/+(.*)\"","\\1\"\\2/\\3\"",
             readLines(devtools::package_file("src/tran.g.d_parser.c")))
sink(devtools::package_file("src/tran.g.d_parser.c"))
cat(paste(file,collapse="\n"));
sink();
## sink(devtools::package_file("R/version.R"))
## cat("##\' Version and repository for this dparser package.
## ##\'
## ##\' @return A character vector with the version and repository.
## ##\' @author Matthew L. Fidler
## ##\' @keywords internal
## ##\' @export
## rxVersion <- function(){return(c(version=\"");
## ver <- readLines(devtools::package_file("DESCRIPTION"));
## ver <- ver[regexpr("^Version:", ver) != -1]
## ver <- ver[1];
## cat(gsub("Version: +", "", ver))
## cat("\",repo=\"");
## cat("https://github.com/")
## tmp <- readLines(devtools::package_file(".git/config"))
## cat(gsub("\\.git$", "", gsub(".*git@github.com:", "", tmp[which(tmp == '[remote "origin"]')[1]+1])))
## cat("\"))}\n");
## sink();
## ## devtools::load_all();

cat("Update README\n");
owd <- getwd();
on.exit({setwd(owd)});
setwd(devtools::package_file());
knitr::knit(devtools::package_file("README.Rmd"))

gen.ome <- function(mx){
    ret <- paste0(sprintf("//Generated from refresh.R for %s dimensions\n#include <R.h>\n#include <Rdefines.h>\n#include <R_ext/Error.h>\n#include <Rmath.h>\nSEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn){\nint dm=INTEGER(dms)[0];\nif (dm == 0){\n  SEXP ret=  PROTECT(allocVector(INTSXP,1));\n  INTEGER(ret)[0] = %s;\n  UNPROTECT(1);\n  return(ret);\n}", mx, mx),
                  paste(sapply(1:mx, function(x){
                      tmp <- rxSymInvC2(matrix(rep(1, x * x), x), allow.cache=FALSE)[[1]]
                      sprintf("else if (dm == %s){\n%s\n}\n", x, tmp);
                  }), collapse=""), "\n  return R_NilValue;\n}");
    sink(devtools::package_file("src/omegaChol.c"));
    cat(ret)
    sink();
}

## create.syminv.cache();

## create.syminv.cache(2);

## create.syminv.cache(3);

## create.syminv.cache(4);

## cpp code

if (Sys.getenv("RxODE_derivs") == "TRUE"){

    gen.ome(12);

}

document();

