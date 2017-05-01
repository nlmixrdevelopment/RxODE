library(devtools)
## library(tidyr)
## library(dplyr)
cat("Copy header to inst directory")

file.copy(devtools::package_file("src/RxODE_types.h"),
          devtools::package_file("inst/include/RxODE_types.h"),
          overwrite=TRUE);

cat("Generate header string.\n");
odec <- readLines(devtools::package_file("inst/ode.c"));
w <- which(regexpr("__ODE_SOLVER__", odec) != -1)[1];
ode <- odec[seq(1, w - 1)];
solve <- odec[seq(w, length(odec))];
hd <- sprintf("#define __HD_ODE__ \"%s\\n\"\n#define __HD_SOLVE__ \"%s\\n\"\n",
              paste(gsub("%", "%%", gsub("\"", "\\\\\"", ode)), collapse="\\n"),
              paste(gsub("%", "%%", gsub("\"", "\\\\\"", solve)), collapse="\\n"));
sink(devtools::package_file("src/ode.h"))
cat(hd);
sink();
cat("Generate dplyr and tidyr compatability functions.\n")
tidyr.fns <- c("spread_", "unite_", "separate_", "gather_");
dplyr.fns <- c("sample_frac", "sample_n", "group_by_", "rename_",
               "summarise_", "transmute_", "mutate_", "distinct_",
               "select_", "arrange_", "slice_", "filter_");
## require(tidyr)
## require(dplyr)
sink(package_file("R/rxsolve-gen.R"))
cat("## Generated code from build/refresh.R\n\n");
one <- "data"
for (f in tidyr.fns){
    ## fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
    ## fn <- paste0(fn[-length(fn)], collapse="\n");
    ## one <- eval(parse(text=sprintf("attr(formals(%s)[1],\"names\")", f)));
    cat(sprintf("##' @name %s
##' @export %s.solveRxDll
##'
##' @method %s solveRxDll
##'
##' @title %s for \\code{solveRxDll} object
##' @description compatability function for tidyr
##' @param data Solved ODE, an \\code{solveRxDll} object.
##' @param ... Additional arguments
##'
%s.solveRxDll <- function(%s, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- dplyr::as.tbl(%s)
  return(do.call(getFromNamespace(\"%s\",\"tidyr\"), call, envir = parent.frame(1)));
}\n\n",f, f, f, f, f, one, one, one, f));
}
as.tbls <- c("sample_n", "sample_frac");
## merge?
for (f in dplyr.fns){
    ## fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
    ## one <- eval(parse(text=sprintf("attr(formals(%s)[1],\"names\")", f)))
    if (any(f == as.tbls)){
        one <- "tbl";
    } else {
        one <- ".data";
    }
    cat(sprintf("##' @name %s
##' @export %s.solveRxDll
##'
##' @method %s solveRxDll
##' @description compatability function for dplyr
##' @title %s for \\code{solveRxDll} object
##' @param %s Solved equation, an \\code{solveRxDll} object.
##' @param ... Additional arguments
##'
%s.solveRxDll <- function(%s, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- %s(%s)
  return(do.call(getFromNamespace(\"%s\",\"dplyr\"), call, envir = parent.frame(1)));
}\n\n",f, f, f, f, one, f, one, one, ifelse(one == "tbl", "dplyr::as.tbl",
                                ifelse(any(f == as.tbls),
                                       "dplyr::as.tbl", "asTbl")),  one, f))
}

for (f in c("row.names", "by", "aggregate", "anyDuplicated", "droplevels",
            "duplicated", "edit", "is.na", "Math", "rowsum",
            "split", "subset", "stack", "unstack", "unique",
            "within", "with")){
    if (f == "Math"){
        fn <- "function(x, ...)";
    } else {
        fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
        fn <- paste0(fn[-length(fn)], collapse="\n");
    }
    one <- eval(parse(text=sprintf("attr(formals(%s)[1],\"names\")", f)))
    ns <- gsub("package:", "", find(f));
    if (!any(ns == c("base"))){
        cat(sprintf("##' @importFrom %s %s\n", ns, f));
    }
    cat(sprintf("##' @export
%s.solveRxDll <- %s{
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- as.data.frame(%s)
  return(do.call(\"%s.data.frame\", call, envir = parent.frame(1)));
}\n\n",f, fn, one, one, f));
}

for (f in c("dim", "dimnames", "t")){
    fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
    fn <- paste0(fn[-length(fn)], collapse="\n");
    one <- "x"
    cat(sprintf("##' @export
%s.solveRxDll <- %s{
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- as.matrix(%s)
  return(do.call(\"%s\", call, envir = parent.frame(1)));
}\n\n",f, fn, one, one, f));
}
sink();

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
sink(devtools::package_file("R/version.R"))
cat("##\' Version and repository for this dparser package.
##\'
##\' @return A character vector with the version and repository.
##\' @author Matthew L. Fidler
##\' @keywords internal
##\' @export
rxVersion <- function(){return(c(version=\"");
ver <- readLines(devtools::package_file("DESCRIPTION"));
ver <- ver[regexpr("^Version:", ver) != -1]
ver <- ver[1];
cat(gsub("Version: +", "", ver))
cat("\",repo=\"");
cat("https://github.com/")
tmp <- readLines(devtools::package_file(".git/config"))
cat(gsub("\\.git$", "", gsub(".*git@github.com:", "", tmp[which(tmp == '[remote "origin"]')[1]+1])))
cat("\"))}\n");
sink();
devtools::load_all();

document();

cat("Update README\n");
owd <- getwd();
on.exit({setwd(owd)});
setwd(devtools::package_file());
knitr::knit(devtools::package_file("README.Rmd"))

create.syminv.cache <- function(n=1, xforms=c("sqrt", "log", "identity")){
    owd <- getwd();
    on.exit(setwd(owd))
    setwd(devtools::package_file("inst/inv"));
    i <- 0;
    if (n == 1){
        m <- matrix(1, 1);
        for (xform in c("sqrt", "log", "identity")){
            rxSymInvC(m, xform);
        }
    } else {
        m <- matrix(rep(1, n * n), n);
        m[lower.tri(m)] <- 0;
        m[upper.tri(m)] <- 0;
        n2 <- sum(lower.tri(m) * 1);
        grd <- eval(parse(text=sprintf("expand.grid(%s)", paste(rep("c(0,1)", n2), collapse=", "))));
        for (i in 1:length(grd[, 1])){
            v <- as.vector(m);
            vw <- which(as.vector(lower.tri(m) * 1) == 1);
            v[vw] <- grd[i, ];
            mc <- matrix(v, n);
            mc[upper.tri(mc)] <- t(mc)[upper.tri(mc)];
            dmat <- dim(mc)[1] -1;
            block <- list();
            last <- 1;
            if (dmat != 0){
                for (i in 1:dmat){
                    if (all(mc[rxBlockZeros(mc,i)] == 0)){
                        s <- seq(last, i);
                        cur <-matrix(as.double(mc[s, s]), length(s));
                        last <- i + 1;
                        block[[length(block) + 1]] <- cur;
                    }
                }
            }
            if (length(block) == 0){
                for (xform in xforms){
                    rxSymInvC(mc, xform);
                    rxSymInvC(mc, xform, chol=TRUE);
                    cat(".");
                    i <- i + 1;
                    if (i %% 5 == 0){
                        cat(i);
                    }
                    if (i %% 50 == 0){
                        cat("\n");
                    }
                }
            }
            ## print(mc);
            ## mc[upper.tri(mc)] <- grd[i, ]
            ## mc[lower.tri(mc)] <- t(mc)[lower.tri(mc)];
            ## print(mc);
        }
        ## print(n2)
        ## print(m);
        ## print(grd);
    }
    cat("\n");
}

## create.syminv.cache();

## create.syminv.cache(2);

## create.syminv.cache(3);

## create.syminv.cache(4);

## cpp code

if (Sys.getenv("RxODE_derivs") == "TRUE"){

    cat("Generating deriatives for solved linear models:");
    i <- 0;
    cur.par <- 1;
    create.cpp.code <- function(model){
        rxSymPySetup(model);
        p <- rxParams(model);
        lhs <- rxLhs(model);
        lhs <- lhs[regexpr(rex::rex(or("alpha", "beta", "gamma", "A", "B", "C")), lhs) != -1]
        ncmt <- length(lhs) / 2;
        ks <- c("alpha", "beta", "gamma");
        cfs <- c("A", "B", "C");
        up <- gsub("K21P", "K21", gsub("BETAP", "BETA", gsub("ALPHAP", "ALPHA", toupper(p))));
        lines <- c();
        lines[length(lines) + 1] <- sprintf("if (par == %s && ncmt == %s && oral == %s){\n", cur.par, ncmt, ifelse(any(up == "KA"), "1", "0"));
        lines[length(lines) + 1] <- sprintf("  if (%s){\n", paste(sprintf("e.exists(\"%s\")", up), collapse=" && "));
        lines[length(lines) + 1] <- paste0(paste(sprintf("    double %s = as<double>(e[\"%s\"]);", p, up), collapse="\n"), "\n");
        lines[length(lines) + 1] <- paste0(paste(sprintf("    mat d%s = mat(%s, 2);", up, ncmt), collapse="\n"), "\n");
        for (curP in p){
            for (curL in lhs){
                cat(".")
                i <<- i + 1;
                if (i %% 5 == 0){
                    cat(i)
                }
                if (i %% 50 == 0){
                    cat("\n");
                }
                tmp <- rxSymPy(sprintf("diff(%s, %s)",
                                       rxToSymPy(curL),
                                       rxToSymPy(curP)))
                if (any(curL == ks)){
                    par <- sprintf("(%s, 0)", which(curL == ks) - 1);
                } else {
                    par <- sprintf("(%s, 1)", which(curL == cfs) - 1)
                }
                tmp <- sprintf("    d%s%s = %s;", gsub("K21P", "K21", gsub("BETAP", "BETA", gsub("ALPHAP", "ALPHA", toupper(curP)))), par, sympyC(tmp));
                lines[length(lines) + 1] <- paste0(tmp, "\n");
            }
        }
        lines[length(lines) + 1] <- paste0(paste(sprintf("    e[\"d%s\"] = d%s;", up, up), collapse="\n"), "\n");
        lines[length(lines) + 1] <- paste0(sprintf("  return;\n  } else {\n    stop(\"Some required parameters not in environment (need %s)\");\n  }\n}\n",
                                                   paste(up, collapse=", ")));
        rxSymPyClean();
        return(lines);
    }

    cpp <- c()

    lin.1cmt.p1.oral <- RxODE({
        volume = V
        k = CL / volume
        alpha = k;
        A = ka / (ka - alpha)/volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.1cmt.p1.oral))

    lin.1cmt.p1 <- RxODE({
        volume = V
        k = CL / volume
        alpha = k;
        A = 1/volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.1cmt.p1))

    lin.2cmt.p1 <- RxODE({
        volume = V
        k = CL / volume
        k12 = Q / volume
        k21 = Q / V2
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = (alpha - k21) / (alpha - beta) / volume;
        B     = (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p1));

    lin.2cmt.p1.oral <- RxODE({
        volume = V
        k = CL / volume
        k12 = Q / volume
        k21 = Q / V2
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = ka / (ka - alpha) * (alpha - k21) / (alpha - beta) / volume;
        B     = ka / (ka - beta) * (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p1.oral));


    ## Fixme 3 cmt returns zoo :(
    lin.3cmt.p1 <- RxODE({
        volume = V
        k = CL / volume
        k12 = Q / volume
        k21 = Q / V2
        k13     = Q2 / volume
        k31     = Q2 / V3
        a0      = k * k21 * k31;
        a1      = k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2      = k + k12 + k13 + k21 + k31;
        p       = a1 - a2 * a2 / 3;
        q       = 2 * a2 * a2 * a2 / 27 - a1 * a2 /3 + a0;
        r1      = sqrt(-p * p * p / 27);
        r2      = 2 * r1 ^ (1 / 3);
        theta   = acos(-q / (2 * r1)) / 3
        alpha   = -(cos(theta) * r2 - a2 / 3);
        beta    = -(cos(theta + 2 / 3 * pi) * r2 - a2 / 3);
        gamma   = -(cos(theta + 4 / 3 * pi) * r2 - a2 / 3);
        A       = (k21 - alpha) * (k31 - alpha) / (alpha - beta) / (alpha - gamma) / volume;
        B       = (k21 - beta) * (k31 - beta) / (beta - alpha) / (beta - gamma) / volume;
        C       = (k21 - gamma) * (k31 - gamma) / (gamma - alpha) / (gamma - beta) / volume;
    })

    ## p and theta do not seem to work with sympy.
    ## also cos(th + 2/3*pi) does not give the same results as as cos(th+2*pi/3)
    ## This is due to the sympy numbers :(
    ## See http://docs.sympy.org/dev/gotchas.html#python-numbers-vs-sympy-numbers
    cpp <- c(cpp, create.cpp.code(lin.3cmt.p1));

    lin.3cmt.p1.oral <- RxODE({
        volume = V
        k = CL / volume
        k12 = Q / volume
        k21 = Q / V2
        k13     = Q2 / volume
        k31     = Q2 / V3
        a0      = k * k21 * k31;
        a1      = k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2      = k + k12 + k13 + k21 + k31;
        pp       = a1 - a2 * a2 / 3;
        q       = 2 * a2 * a2 * a2 / 27 - a1 * a2 /3 + a0;
        r1      = sqrt(-pp * pp * pp / 27);
        r2      = 2 * r1 ^ (1 / 3);
        th   = acos(-q / (2 * r1)) / 3
        alpha   = -(cos(th) * r2 - a2 / 3);
        beta    = -(cos(th + 2 * pi / 3) * r2 - a2 / 3);
        gamma   = -(cos(th + 4 * pi / 3) * r2 - a2 / 3);
        A       = ka / (ka - alpha) * (k21 - alpha) * (k31 - alpha) / (alpha - beta) / (alpha - gamma) / volume;
        B       = ka / (ka - beta) * (k21 - beta) * (k31 - beta) / (beta - alpha) / (beta - gamma) / volume;
        C       = ka / (ka - gamma) * (k21 - gamma) * (k31 - gamma) / (gamma - alpha) / (gamma - beta) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.3cmt.p1.oral));

    ## parameterization #2
    cur.par <- 2;

    lin.1cmt.p2.oral <- RxODE({
        volume = V
        alpha = k;
        A = ka / (ka - alpha)/volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.1cmt.p2.oral))

    lin.1cmt.p2 <- RxODE({
        volume = V
        alpha = k;
        A = 1/volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.1cmt.p2))

    lin.2cmt.p2 <- RxODE({
        volume = V
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = (alpha - k21) / (alpha - beta) / volume;
        B     = (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p2));

    lin.2cmt.p2.oral <- RxODE({
        volume = V
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = ka / (ka - alpha) * (alpha - k21) / (alpha - beta) / volume;
        B     = ka / (ka - beta) * (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p2.oral));

    lin.3cmt.p2 <- RxODE({
        volume = V
        a0      = k * k21 * k31;
        a1      = k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2      = k + k12 + k13 + k21 + k31;
        pp       = a1 - a2 * a2 / 3;
        q       = 2 * a2 * a2 * a2 / 27 - a1 * a2 /3 + a0;
        r1      = sqrt(-pp * pp * pp / 27);
        r2      = 2 * r1 ^ (1 / 3);
        th   = acos(-q / (2 * r1)) / 3
        alpha   = -(cos(th) * r2 - a2 / 3);
        beta    = -(cos(th + 2 / 3 * pi) * r2 - a2 / 3);
        gamma   = -(cos(th + 4 / 3 * pi) * r2 - a2 / 3);
        A       = (k21 - alpha) * (k31 - alpha) / (alpha - beta) / (alpha - gamma) / volume;
        B       = (k21 - beta) * (k31 - beta) / (beta - alpha) / (beta - gamma) / volume;
        C       = (k21 - gamma) * (k31 - gamma) / (gamma - alpha) / (gamma - beta) / volume;
    })

    ## p and theta do not seem to work with sympy.
    cpp <- c(cpp, create.cpp.code(lin.3cmt.p2));

    lin.3cmt.p2.oral <- RxODE({
        volume = V
        a0      = k * k21 * k31;
        a1      = k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2      = k + k12 + k13 + k21 + k31;
        pp       = a1 - a2 * a2 / 3;
        q       = 2 * a2 * a2 * a2 / 27 - a1 * a2 /3 + a0;
        r1      = sqrt(-pp * pp * pp / 27);
        r2      = 2 * r1 ^ (1 / 3);
        theta   = acos(-q / (2 * r1)) / 3
        alpha   = -(cos(theta) * r2 - a2 / 3);
        beta    = -(cos(theta + 2 / 3 * pi) * r2 - a2 / 3);
        gamma   = -(cos(theta + 4 / 3 * pi) * r2 - a2 / 3);
        A       = ka / (ka - alpha) * (k21 - alpha) * (k31 - alpha) / (alpha - beta) / (alpha - gamma) / volume;
        B       = ka / (ka - beta) * (k21 - beta) * (k31 - beta) / (beta - alpha) / (beta - gamma) / volume;
        C       = ka / (ka - gamma) * (k21 - gamma) * (k31 - gamma) / (gamma - alpha) / (gamma - beta) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.3cmt.p2.oral));


    cur.par <- 3;

    lin.2cmt.p3 <- RxODE({
        volume = V
        k = CL / volume
        k12 = Q / volume
        k21 = Q / (VSS - volume)
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = (alpha - k21) / (alpha - beta) / volume;
        B     = (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p3));

    lin.2cmt.p3.oral <- RxODE({
        volume = V
        k = CL / volume
        k12 = Q / volume
        k21 = Q / (VSS - volume)
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = ka / (ka - alpha) * (alpha - k21) / (alpha - beta) / volume;
        B     = ka / (ka - beta) * (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p3.oral));

    cur.par <- 4;

    lin.2cmt.p4 <- RxODE({
        volume = V
        k21 = (aob*betaP+alphaP)/(aob+1);
        k   = (alphaP*betaP)/k21;
        k12 = alphaP+betaP-k21-k;
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = (alpha - k21) / (alpha - beta) / volume;
        B     = (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p4));

    lin.2cmt.p4.oral <- RxODE({
        volume = V
        k21 = (aob*betaP+alphaP)/(aob+1);
        k   = (alphaP*betaP)/k21;
        k12 = alphaP+betaP-k21-k;
        beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4 * k21 * k));
        alpha = k21 * k / beta;
        A     = ka / (ka - alpha) * (alpha - k21) / (alpha - beta) / volume;
        B     = ka / (ka - beta) * (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p4.oral));

    cur.par <- 5;

    lin.2cmt.p5 <- RxODE({
        volume = V
        beta  = betap
        alpha = alphap
        k21 = k21p
        A     = (alpha - k21) / (alpha - beta) / volume;
        B     = (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p5));

    lin.2cmt.p5.oral <- RxODE({
        volume = V
        beta  = betap
        alpha = alphap
        k21 = k21p
        A     = ka / (ka - alpha) * (alpha - k21) / (alpha - beta) / volume;
        B     = ka / (ka - beta) * (beta - k21) / (beta - alpha) / volume;
    })

    cpp <- c(cpp, create.cpp.code(lin.2cmt.p5.oral));

    sink(devtools::package_file("src/lincmtDiff.cpp"));
    cat("// [[Rcpp::depends(RcppArmadillo)]]
// Generated by refresh.R;Can be recreated by refresh(derivs=TRUE) in source directory after loading RxODE by library(devtools);load_all();
#include <RcppArmadillo.h>
#include <R.h>
using namespace Rcpp;
using namespace R;
using namespace arma;
// [[Rcpp::export]]
void getLinDerivs(SEXP rho){
  Environment e = as<Environment>(rho);
  int par = as<int>(e[\"parameterization\"]);
  int ncmt = as<int>(e[\"ncmt\"]);
  int oral = as<int>(e[\"oral\"]);
  // double zoo = R_PosInf; // Zoo is in dV :(
");
    cat(paste(cpp, collapse="\n"));
    cat(" stop(\"environment not setup properly for this function.\");\n");
    cat("}\n");
    sink()
    rxClean();
    load_all();
    cat("done.\n");
}
