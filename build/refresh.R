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

    lin.diff <- function(){
        ## These are the derivatives in infusion
        tmp <- RxODE({
            ret = rate*A / alpha * (1 - exp(-alpha * t1)) * exp(-alpha * t2);
        })

        rxSymPySetup(tmp)
        dA <- rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, A)")))

        rxSymPySetup(tmp)
        dA.dA <- sprintf("ret=%s", rxSymPy("diff(diff(ret, A),A)"));

        rxSymPySetup(tmp)
        dA.dAlpha <- rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A),alpha)")));

        rxSymPySetup(tmp)
        dAlpha = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, alpha)")))

        rxSymPySetup(tmp)
        dAlpha.dA = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha),A)")));

        rxSymPySetup(tmp)
        dAlpha.dAlpha = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha),alpha)")));

        tmp <- RxODE({
            ret = rate*A / alpha * (1 - exp(-alpha * (tT- tlag)));
        })

        ## Now tlag derivatives during infusion  (i.e. t2=0)

        rxSymPySetup(tmp)
        dTlag.Inf = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, tlag)")));

        rxSymPySetup(tmp)
        dTlag.dA.Inf = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag),A)")));

        rxSymPySetup(tmp)
        dTlag.dAlpha.Inf = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag),alpha)")));

        rxSymPySetup(tmp)
        dTlag.dTlag.Inf = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag),tlag)")));

        rxSymPySetup(tmp) ## These are likely exactly the same...
        dAlpha.dTlag.Inf = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha),tlag)")));

        rxSymPySetup(tmp) ## These are likely exactly the same...
        dA.dTlag.Inf = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A),tlag)")));

        ## Now tlag derivatives during after infusion  (i.e. t2!=0)
        tmp <- RxODE({
            ret = rate*A / alpha * (1 - exp(-alpha * (tinfA - tlag))) * exp(-alpha * (tT- tlag));
        })

        rxSymPySetup(tmp)
        dTlag = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, tlag)")));

        rxSymPySetup(tmp)
        dTlag.dA = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag),A)")));

        rxSymPySetup(tmp)
        dTlag.dAlpha = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag),alpha)")));

        rxSymPySetup(tmp)
        dTlag.dTlag = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), tlag)")));

        rxSymPySetup(tmp) ## These are likely exactly the same...
        dAlpha.dTlag = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha),tlag)")));

        rxSymPySetup(tmp) ## These are likely exactly the same...
        dA.dTlag = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A),tlag)")));

        ## Now oral derivatives
        tmp <- RxODE({
            ret = dose * A *(exp(-alpha * (tT - tlag)) - exp(-ka * (tT - tlag)))
        })

        rxSymPySetup(tmp)
        dA.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, A)")));

        rxSymPySetup(tmp)
        dA.dA.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), A)")));

        rxSymPySetup(tmp)
        dA.dAlpha.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), alpha)")));

        rxSymPySetup(tmp)
        dA.dKa.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), ka)")));

        rxSymPySetup(tmp)
        dA.dTlag.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), tlag)")));

        rxSymPySetup(tmp)
        dAlpha.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, alpha)")));

        rxSymPySetup(tmp)
        dAlpha.dA.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), A)")));

        rxSymPySetup(tmp)
        dAlpha.dAlpha.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), alpha)")));

        rxSymPySetup(tmp)
        dAlpha.dKa.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), ka)")));

        rxSymPySetup(tmp)
        dAlpha.dTlag.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), tlag)")));

        rxSymPySetup(tmp)
        dTlag.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, tlag)")));

        rxSymPySetup(tmp)
        dTlag.dA.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), A)")));

        rxSymPySetup(tmp)
        dTlag.dAlpha.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), alpha)")));

        rxSymPySetup(tmp)
        dTlag.dKa.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), ka)")));

        rxSymPySetup(tmp)
        dTlag.dTlag.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), tlag)")));

        rxSymPySetup(tmp)
        dKa.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, ka)")));

        rxSymPySetup(tmp)
        dKa.dA.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, ka),A)")));

        rxSymPySetup(tmp)
        dKa.dAlpha.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, ka),alpha)")));

        rxSymPySetup(tmp)
        dKa.dKa.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, ka),ka)")));

        rxSymPySetup(tmp)
        dKa.dTlag.Oral = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, ka), tlag)")));

        ## Bolus model
        tmp <- RxODE({
            ret = dose * A *exp(-alpha * (tT - tlag))
        })

        rxSymPySetup(tmp)
        dA.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, A)")));

        rxSymPySetup(tmp)
        dA.dA.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), A)")));

        rxSymPySetup(tmp)
        dA.dAlpha.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), alpha)")));

        rxSymPySetup(tmp)
        dA.dTlag.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, A), tlag)")));

        rxSymPySetup(tmp)
        dAlpha.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, alpha)")));

        rxSymPySetup(tmp)
        dAlpha.dA.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), A)")));

        rxSymPySetup(tmp)
        dAlpha.dAlpha.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), alpha)")));

        rxSymPySetup(tmp)
        dAlpha.dTlag.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, alpha), tlag)")));

        rxSymPySetup(tmp)
        dTlag.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(ret, tlag)")));

        rxSymPySetup(tmp)
        dTlag.dA.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), A)")));

        rxSymPySetup(tmp)
        dTlag.dAlpha.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), alpha)")));

        rxSymPySetup(tmp)
        dTlag.dTlag.Bolus = rxLogifyModel(sprintf("ret=%s", rxSymPy("diff(diff(ret, tlag), tlag)")));

        reg <- rex::rex(any_spaces, or("\n", ""), any_spaces, end);
        regRx <- rex::rex(capture(or("sign_exp(", "abs_log(", "safe_zero(",
                                     "abs_log1p(", "podo(", "tlast(", "factorial(")))

        writeLines(c("// Generated by refresh.R;Can be recreated by refresh(derivs=TRUE) in source directory after loading RxODE by library(devtools);load_all();",
                     "#include <R.h>",
                     "#include <Rinternals.h>",
                     "#define dTlag 8",
                     "#define dKa 7","",
                     "extern double RxODE_safe_zero(double x);",
                     "extern double RxODE_sign_exp(double sgn, double x);",
                     "extern double RxODE_abs_log(double x);",
                     "extern double RxODE_abs_log1p(double x);",
                     "extern double RxODE_podo();",
                     "extern double RxODE_tlast();",
                     "extern double RxODE_factorial(double x);","",
                     ##
                     "extern double rxSolveLinBdInf(int diff1, int diff2, int dA, int dAlpha, double rate, double tT, double t1, double t2, double tinf, double A, double alpha, double tlag){",
                     "double ret = 0, tinfA=tinf+tlag;",
                     ##
                     "if ((diff1 == dA && diff2 == 0) || (diff1 == 0 && diff2 == dA)) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA)),
                     "} else if (diff1 == dA && diff2 == dA){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dA)),
                     "} else if (diff1 == dA && diff2 == dAlpha){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dAlpha)),
                     "} else if (diff1 == dA && diff2 == dTlag && t2 == 0){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dTlag.Inf)),
                     "} else if (diff1 == dA && diff2 == dTlag && t2 != 0){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dTlag)),
                     ##
                     "} else if ((diff1 == dAlpha && diff2 == 0) || (diff1 == 0 && diff2 == dAlpha)) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha)),
                     "} else if (diff1 == dAlpha && diff2 == dA){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dA)),
                     "} else if (diff1 == dAlpha && diff2 == dAlpha){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dAlpha)),
                     "} else if (diff1 == dAlpha && diff2 == dTlag && t2 == 0){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dTlag.Inf)),
                     "} else if (diff1 == dAlpha && diff2 == dTlag && t2 != 0){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dTlag)),
                     ##
                     "} else if (t2 == 0 && ((diff1 == dTlag && diff2 == 0) || (diff1 == 0 && diff2 == dTlag))){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.Inf)),
                     "} else if (t2 != 0 && ((diff1 == dTlag && diff2 == 0) || (diff1 == 0 && diff2 == dTlag))){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag)),
                     "} else if (t2 == 0 && diff1 == dTlag && diff2 == dA){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.dA.Inf)),
                     "} else if (t2 != 0 && diff1 == dTlag && diff2 == dA){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.dA)),
                     "} else if (t2 == 0 && diff1 == dTlag && diff2 == dAlpha){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.dAlpha.Inf)),
                     "} else if (t2 != 0 && diff1 == dTlag && diff2 == dAlpha){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.dAlpha)),
                     "} else if (t2 == 0 && diff1 == dTlag && diff2 == dTlag){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.dTlag.Inf)),
                     "} else if (t2 != 0 && diff1 == dTlag && diff2 == dTlag){",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.dTlag)),
                     "} else { return 0 ;}",
                     "return ret;}", "",
                     ##
                     "extern double rxSolveLinBDiff(int diff1, int diff2, int dA, int dAlpha, double dose, double tT, double A, double alpha, double ka, double tlag){",
                     "double ret = 0;",
                     ##
                     "if (ka > 0 && ((diff1 == dA && diff2 == 0) || (diff1 == 0 && diff2 == dA))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.Oral)),
                     "} else if (ka > 0 && ((diff1 == dAlpha && diff2 == 0) || (diff1 == 0 && diff2 == dAlpha))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.Oral)),
                     "} else if (ka > 0 && ((diff1 == dKa && diff2 == 0) || (diff1 == 0 && diff2 == dKa))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dKa.Oral)),
                     "} else if (ka > 0 && ((diff1 == dTlag && diff2 == 0) || (diff1 == 0 && diff2 == dTlag))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.Oral)),
                     ##
                     "} else if (ka > 0 && diff1 == dA && diff2 == dA) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dA.Oral)),
                     "} else if (ka > 0 && diff1 == dA && diff2 == dAlpha) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dAlpha.Oral)),
                     "} else if (ka > 0 && diff1 == dA && diff2 == dKa) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dKa.Oral)),
                     "} else if (ka > 0 && diff1 == dA && diff2 == dTlag) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dTlag.Oral)),
                     ##
                     "} else if (ka > 0 && diff1 == dAlpha && diff2 == dA) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dA.Oral)),
                     "} else if (ka > 0 && diff1 == dAlpha && diff2 == dAlpha) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dAlpha.Oral)),
                     "} else if (ka > 0 && diff1 == dAlpha && diff2 == dKa) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dKa.Oral)),
                     "} else if (ka > 0 && diff1 == dAlpha && diff2 == dTlag) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dTlag.Oral)),
                     ##
                     "} else if (ka <= 0 && ((diff1 == dA && diff2 == 0) || (diff1 == 0 && diff2 == dA))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.Bolus)),
                     "} else if (ka <= 0 && ((diff1 == dAlpha && diff2 == 0) || (diff1 == 0 && diff2 == dAlpha))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.Bolus)),
                     "} else if (ka <= 0 && ((diff1 == dTlag && diff2 == 0) || (diff1 == 0 && diff2 == dTlag))) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dTlag.Bolus)),
                     ##
                     "} else if (ka <= 0 && diff1 == dA && diff2 == dA) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dA.Bolus)),
                     "} else if (ka <= 0 && diff1 == dA && diff2 == dAlpha) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dAlpha.Bolus)),
                     "} else if (ka <= 0 && diff1 == dA && diff2 == dTlag) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dA.dTlag.Bolus)),
                     ##
                     "} else if (ka <= 0 && diff1 == dAlpha && diff2 == dA) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dA.Bolus)),
                     "} else if (ka <= 0 && diff1 == dAlpha && diff2 == dAlpha) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dAlpha.Bolus)),
                     "} else if (ka <= 0 && diff1 == dAlpha && diff2 == dTlag) {",
                     gsub(regRx, "RxODE_\\1", gsub(reg, ";", dAlpha.dTlag.Bolus)),
                     ##
                     "} else { return 0 ;}",
                     "return ret;}"
                     ), devtools::package_file("src/lincmtDiff.c"))

        rxClean();
        ## load_all();
        cat("done.\n");
    }
}

document();

