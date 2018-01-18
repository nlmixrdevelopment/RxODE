library(devtools)

if (!file.exists(devtools::package_file("src/liblsoda"))){
    owd <- getwd();
    setwd(devtools::package_file("src"));
    system("git clone git@github.com:sdwfrost/liblsoda")
    setwd(owd);
}

for (file in list.files(devtools::package_file("src/liblsoda/src"), pattern="[.]([hc]|inc)$")){
    message(sprintf("\tcopy %s", file))
    lines <- suppressWarnings(readLines(file.path(devtools::package_file("src/liblsoda/src"), file)));
    lines <- gsub("#include \"cfode_static.inc\"", "#include \"cfode_static.h\"", lines, fixed=TRUE);
    lines <- gsub("fprintf[(] *stderr *, *", "Rprintf(", lines);
    lines <- gsub("([ \t])printf[(] *", "\\1Rprintf(", lines);
    if (any(regexpr("Rprintf", lines) != -1)){
        lines <- c("#include <R.h>", "#include <Rinternals.h>", lines);
    }
    if (any(regexpr("#define ERROR", lines, fixed=TRUE) != -1)){
        lines <- gsub("#define ERROR", "#ifdef ERROR\n#undef ERROR\n#endif\n#define ERROR", lines, fixed=TRUE);
    }
    if (file == "solsy.c"){
        lines <- c("#include <R.h>", "#include <Rinternals.h>", gsub("abort[(] *[)]", "error(\"liblsoda does not implement this. (solsy)\")", lines));
    }
    if (file == "lsoda.c"){
        lines <- gsub("int[ \t]+i[ \t]*([;,])", "int i=0\\1", lines);
    }
    writeLines(lines, file.path(devtools::package_file("src"), gsub("[.]inc$", ".h", file)))
}


## library(tidyr)
## library(dplyr)
cat("Copy header to inst directory")

file.copy(devtools::package_file("src/RxODE_types.h"),
          devtools::package_file("inst/include/RxODE_types.h"),
          overwrite=TRUE);

cat("Generate header string.\n");
odec <- readLines(devtools::package_file("inst/ode.c"));
solvec <- readLines(devtools::package_file("src/solve.h"));
w <- which(regexpr("#define R_pow_di Rx_pow_di", odec, fixed=TRUE) != -1)[1];
odec <- c(odec[1:w], solvec, odec[-(1:w)])
w <- which(regexpr("// CODE HERE", odec) != -1)[1];
ode <- odec[seq(1, w - 1)];
solve <- odec[seq(w + 1, length(odec))];
solve <- paste(gsub("%", "%%", gsub("\"", "\\\\\"", solve)), collapse="\\n")
if (nchar(solve) > 4095){
    solve1 <- substr(solve, 1, 4094);
    solve2 <- substr(solve, 4095, nchar(solve))
} else {
    solve1 <- solve;
    solve2 <- "";
}

found <- FALSE
hd <- sapply(strsplit(sprintf("#define __HD_ODE__ \"%s\\n\"\n#define __HD_SOLVE1__ \"%s\"\n#define __HD_SOLVE2__ \"%s\"",
              paste(gsub("%", "%%", gsub("\"", "\\\\\"", ode)), collapse="\\n"),
              solve1, solve2), "\n")[[1]],
             function(s){
    if (found){
        s <- gsub("#define __HD_SOLVE2__ \"n", "#define __HD_SOLVE2__ \"\\n", s, fixed=TRUE)
        found <<- FALSE
    }
    r1 <- substr(s, 0, nchar(s) - 2)
    r2 <- substr(s, nchar(s) - 1, nchar(s));
    if (r2 == "\\\""){
        found <<- TRUE
        return(paste0(r1, "\""))
    } else {
        return(paste0(r1, r2))
    }
});

writeLines(hd, devtools::package_file("src/ode.h"))

cat("Update Parser c file\n");
dparser::mkdparse(devtools::package_file("inst/tran.g"),
         devtools::package_file("src/"),
         grammar_ident="RxODE");
file <- gsub("^([#]line [0-9]+ )\".*(src)/+(.*)\"","\\1\"\\2/\\3\"",
             readLines(devtools::package_file("src/tran.g.d_parser.c")))
sink(devtools::package_file("src/tran.g.d_parser.c"))
cat(paste(file,collapse="\n"));
sink();
unlink(devtools::package_file("src/tran.o"))
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
                      sprintf("else if (dm == %s){\n%s\n}\n", x, rxSymInvC2(matrix(rep(1, x * x), x), allow.cache=FALSE));
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

    lin.diff <- function(logify=TRUE){
        ## These are the derivatives in infusion
        reg <- rex::rex(any_spaces, or(";", ""), or("\n", ""), any_spaces, end);
        rxLogifyModel <- function(x){
            ret <- RxODE(strsplit(gsub(" +", "", rxSumProdModel(x, sum=TRUE, prod=TRUE)), "\n")[[1]])
            ret <- gsub(rex::rex("(__0__)"), "", gsub(reg, ";", rxModelVars(ret)$model["parseModel"]));
            cat(".")
            return(ret)
        }

        tmp <- RxODE({
            ret = rate*A / alpha * (1 - exp(-alpha * t1)) * exp(-alpha * t2);
        })

        vars <- c("A", "alpha")

        one <- sapply(vars, function(x){
            rxSymPySetup(tmp)
            diff <- sprintf("diff(ret, %s)", x);
            diff <- rxSymPy(diff);
            ret <- rxLogifyModel(sprintf("ret=%s", diff));
            diff <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
            ret <- sprintf("  } else if ((diff1 == %s && diff2 == 0) || (diff1 == 0 && diff2 == %s)){\n    %s", diff, diff, ret);
            return(ret)
        })

        two <- sapply(vars, function(x){
            sapply(vars, function(y) {
                rxSymPySetup(tmp)
                diff <- sprintf("diff(diff(ret, %s),%s)", x, y);
                diff <- rxSymPy(diff);
                ret <- rxLogifyModel(sprintf("ret=%s", diff));
                diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                ret <- sprintf("  } else if (diff1 == %s && diff2 == %s){\n    %s", diff.x, diff.y, ret);
                return(ret)
            })
        })

        ## Now tlag derivatives during infusion  (i.e. t2=0)
        tmp <- RxODE({
            ret = rate*A / alpha * (1 - exp(-alpha * (tT- tlag)));
        })

        rxSymPySetup(tmp);
        diff <- rxSymPy("diff(ret,tlag)");
        ret <- rxLogifyModel(sprintf("ret=%s", diff));
        one <- c(one, sprintf("  } else if (t2 <= 0.0 && ((diff1 == dTlag && diff2 == 0) || (diff1 == 0 && diff2 == dTlag))){\n    %s", ret))

        two <- c(two, sapply(c(vars, "tlag"), function(y) {
                          rxSymPySetup(tmp)
                          diff <- sprintf("diff(diff(ret, tlag),%s)", y);
                          diff <- rxSymPy(diff);
                          ret <- rxLogifyModel(sprintf("ret=%s", diff));
                          ## diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                          diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                          ret.1 <- sprintf("  } else if (diff1 == dTlag && diff2 == %s && t2 <= 0.0){\n    %s", diff.y, ret);
                          rxSymPySetup(tmp)
                          diff <- sprintf("diff(diff(ret, %s), tlag)", y);
                          diff <- rxSymPy(diff);
                          ret <- rxLogifyModel(sprintf("ret=%s", diff));
                          ## diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                          diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                          ret.1 <- sprintf("%s\n  } else if (t2 <= 0.0 && diff1 == %s && diff2 == dTlag){\n    %s", ret, diff.y, ret);
                          return(ret)
                      }))

        ## Now tlag derivatives during after infusion  (i.e. t2!=0)
        tmp <- RxODE({
            ret = rate*A / alpha * (1 - exp(-alpha * (tinfA - tlag))) * exp(-alpha * (tT- tlag));
        })

        rxSymPySetup(tmp);
        diff <- rxSymPy("diff(ret,tlag)");
        ret <- rxLogifyModel(sprintf("ret=%s", diff));
        one <- c(one, sprintf("  } else if (t2 > 0.0 && ((diff1 == dTlag && diff2 == 0) || (diff1 == 0 && diff2 == dTlag))){\n    %s", ret))


        two <- c(two, sapply(c(vars, "tlag"), function(y) {
                          rxSymPySetup(tmp)
                          diff <- sprintf("diff(diff(ret, tlag),%s)", y);
                          diff <- rxSymPy(diff);
                          ret <- rxLogifyModel(sprintf("ret=%s", diff));
                          ## diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                          diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                          ret.1 <- sprintf("  } else if (diff1 == dTlag && diff2 == %s && t2 > 0.0){\n    %s", diff.y, ret);
                          rxSymPySetup(tmp)
                          diff <- sprintf("diff(diff(ret, %s), tlag)", y);
                          diff <- rxSymPy(diff);
                          ret <- rxLogifyModel(sprintf("ret=%s", diff));
                          ## diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                          diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                          ret.1 <- sprintf("%s\n  } else if (t2 > 0.0 && diff1 == %s && diff2 == dTlag){\n    %s", ret, diff.y, ret);
                          return(ret)
                      }))

        infusion <- c(one, two)
        infusion[1] <- sub(rex::rex(start, any_spaces, "} else "), "  ", infusion[1])

        ## Now oral derivatives
        tmp <- RxODE({
            ret = dose * A *(exp(-alpha * (tT - tlag)) - exp(-ka * (tT - tlag)))
        })

        vars <- c("A", "alpha", "ka", "tlag");

        one <- sapply(vars, function(x){
            rxSymPySetup(tmp)
            diff <- sprintf("diff(ret, %s)", x);
            diff <- rxSymPy(diff);
            ret <- rxLogifyModel(sprintf("ret=%s", diff));
            diff <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
            ret <- sprintf("  } else if (ka > 0 && ((diff1 == %s && diff2 == 0) || (diff1 == 0 && diff2 == %s))){\n    %s", diff, diff, ret);
            return(ret)
        })

        two <- sapply(vars, function(x){
            sapply(vars, function(y) {
                rxSymPySetup(tmp)
                diff <- sprintf("diff(diff(ret, %s),%s)", x, y);
                diff <- rxSymPy(diff);
                ret <- rxLogifyModel(sprintf("ret=%s", diff));
                diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                ret <- sprintf("  } else if (ka > 0 && diff1 == %s && diff2 == %s){\n    %s", diff.x, diff.y, ret);
                return(ret)
            })
        })

        ## Bolus model
        tmp <- RxODE({
            ret = dose * A *exp(-alpha * (tT - tlag))
        })

        vars <- c("A", "alpha", "tlag");

        one <- c(one, sapply(vars, function(x){
                          rxSymPySetup(tmp)
                          diff <- sprintf("diff(ret, %s)", x);
                          diff <- rxSymPy(diff);
                          ret <- rxLogifyModel(sprintf("ret=%s", diff));
                          diff <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                          ret <- sprintf("  } else if (ka <= 0 && ((diff1 == %s && diff2 == 0) || (diff1 == 0 && diff2 == %s))){\n    %s", diff, diff, ret);
                          return(ret)
                      }))

        two <- c(two, sapply(vars, function(x){
                          sapply(vars, function(y) {
                              rxSymPySetup(tmp)
                              diff <- sprintf("diff(diff(ret, %s),%s)", x, y);
                              diff <- rxSymPy(diff);
                              ret <- rxLogifyModel(sprintf("ret=%s", diff));
                              diff.x <- paste0("d",toupper(substr(x,0,1)),substr(x,2,nchar(x)));
                              diff.y <- paste0("d",toupper(substr(y,0,1)),substr(y,2,nchar(y)));
                              ret <- sprintf("  } else if (ka <= 0 && diff1 == %s && diff2 == %s){\n    %s", diff.x, diff.y, ret);
                              return(ret)
                          })
                      }))

        oral.bolus <- c(one, two);
        oral.bolus[1] <- sub(rex::rex(start, any_spaces, "} else "), "  ", oral.bolus[1])

        writeLines(c("// Generated by refresh.R;Can be recreated by refresh(derivs=TRUE) in source directory after loading RxODE by library(devtools);load_all();",
                     "#include <R.h>",
                     "#include <Rinternals.h>",
                     "#include <Rmath.h>","",
                     "#define dTlag 8",
                     "#define dKa 7","",
                     "#define _prod RxODE_prodV",
                     "#define _sum  RxODE_sumV",
                     "#define _sign RxODE_signV",
                     "#define Rx_pow RxODE_pow",
                     "#define Rx_pow_di RxODE_pow_di",
                     "#define R_pow RxODE_pow",
                     "#define R_pow_di RxODE_pow_di",
                     "#define safe_zero RxODE_safe_zero","",
                     "extern double RxODE_signV(int n, ...);",
                     "extern double RxODE_prodV(int n, ...);",
                     "extern double RxODE_sumV(int n, ...);",
                     "extern double RxODE_pow(double x, double y);",
                     "extern double RxODE_pow_di(double x, int i);",
                     "extern double RxODE_safe_zero(double x);","",
                     ##
                     "extern double rxSolveLinBdInf(int diff1, int diff2, int dA, int dAlpha, double rate, double tT, double t1, double t2, double tinf, double A, double alpha, double tlag){",
                     "double ret = 0, tinfA=tinf+tlag;",
                     infusion,
                     "  } else {\n    return 0;\n  }\n",
                     "  return ret;\n}", "",
                     ##
                     "extern double rxSolveLinBDiff(int diff1, int diff2, int dA, int dAlpha, double dose, double tT, double A, double alpha, double ka, double tlag){",
                     "  double ret = 0;",
                     oral.bolus,
                     "  } else {\n    return 0;\n  }\n",
                     "  return ret;\n}"
                     ), devtools::package_file("src/lincmtDiff.c"))

        rxClean();
        ## load_all();
        cat("done.\n");
    }

    lin.diff();

}

document();

