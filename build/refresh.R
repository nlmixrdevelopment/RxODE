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
## devtools::load_all();

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

    template.diff <- "if(rate>0){
 if(T<=tinf){
<%=CpIin%>
 }else{
<%=CpIout%>
 }
}else{
<%=Cp%>
}";


    i <- 0;
    oral.tlag <- "exp(-ka * (tT - tlag))"
    oral <- "exp(-ka * T)";

    fn <- function(A, alpha, tlag=F){
        if (tlag){
            c(sprintf("Dose*%s*(exp(-%s*(tT - tlag)))", A, alpha), ## Cp =
              sprintf("Dose*%s*(exp(-%s*(tT - tlag)) - %s)", A, alpha, oral.tlag), ## CpO =
              ## t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
              ## t2 = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion
              sprintf("rate * %s / %s * (1 - exp(-%s* tT))", A, alpha, alpha), ## CpIin
              sprintf("rate * %s / %s * (1 - exp(-%s * tinf))*(exp(-%s*(tT-tinf)))", A, alpha, alpha, alpha) ## CpIout
              )
        } else {
            c(sprintf("Dose*%s*(exp(-%s*T))", A, alpha), ## Cp =
              sprintf("Dose*%s*(exp(-%s*T) - %s)", A, alpha, oral), ## CpO =
              ##t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
              ##t2 = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion
              sprintf("rate * %s / %s * (1 - exp(-%s* T))", A, alpha, alpha), ## CpIin
              sprintf("rate * %s / %s * (1 - exp(-%s * tinf))*(exp(-%s*(T-tinf)))", A, alpha, alpha, alpha) ## CpIout
              )
        }
    }

    mod.1  <- paste(c("Cp=", "CpO=", "CpIin=", "CpIout="), fn("A", "alpha"));
    mod.1.tlag  <- paste(c("Cp=", "CpO=", "CpIin=", "CpIout="), fn("A", "alpha", T));

    mod.2  <- paste(mod.1, "+", fn("B", "beta"));
    mod.2.tlag  <- paste(mod.1.tlag, "+", fn("B", "beta", T));

    mod.3  <- paste(mod.2, "+", fn("C", "gamma"));
    mod.3.tlag  <- paste(mod.2.tlag, "+", fn("C", "gamma", T));

    ## Create RxODE models
    mod.1 <- RxODE(paste(mod.1, collapse="\n"))
    mod.1.tlag <- RxODE(paste(mod.1.tlag, collapse="\n"))

    mod.2 <- RxODE(paste(mod.2, collapse="\n"))
    mod.2.tlag <- RxODE(paste(mod.2.tlag, collapse="\n"))

    mod.3 <- RxODE(paste(mod.3, collapse="\n"))
    mod.3.tlag <- RxODE(paste(mod.3.tlag, collapse="\n"))

    i <- 0;

    cur.p <- c("A", "alpha", "B", "beta", "C", "gamma", "ka", "tlag");
    cur.p.rep <- c("par(0,1)", "par(0,0)", "par(1,1)", "par(1,0)", "par(2,1)", "par(2,0)", "ka", "tlag");

    ipp <- function(){
        cat(".")
        i <<- i + 1;
        if (i %% 5 == 0){
            cat(i)
        }
        if (i %% 50 == 0){
            cat("\n");
        }
    }

    create.cpp.code <- function(totcmt=3){
        sink(devtools::package_file("src/lincmtDiff.cpp"));
        cat("// [[Rcpp::depends(RcppArmadillo)]]
// Generated by refresh.R;Can be recreated by refresh(derivs=TRUE) in source directory after loading RxODE by library(devtools);load_all();
#include <RcppArmadillo.h>
#include <R.h>
using namespace Rcpp;
using namespace R;
using namespace arma;
extern \"C\" double getLinDeriv(int ncmt, int diff1, int diff2, double rate, double tinf, double Dose, double ka, double tlag, double T, double tT, mat par);
double getLinDeriv(int ncmt, int diff1, int diff2, double rate, double tinf, double Dose, double ka, double tlag, double T, double tT, mat par){
double ret = 0;
")
        sink();
        for (ncmt in seq(1, totcmt)){
            for (oral in c(1, 0)){
                lines <- c();
                i <<- 0;
                cat(sprintf("Calculating base %s cmt %s derivatives\n", ncmt, ifelse(oral == 1, "oral", "bolus/infusion")))
                pars <- cur.p[seq(1, ncmt * 2)];
                if (oral) pars <- c(pars, "ka", "tlag");
                for (curP in pars){
                    if (curP == "tlag"){
                        if (ncmt == 1){
                            rxSymPySetup(mod.1.tlag);
                        } else if (ncmt == 2){
                            rxSymPySetup(mod.2.tlag);
                        } else if (ncmt == 3){
                            rxSymPySetup(mod.3.tlag);
                        }
                    } else {
                        if (ncmt == 1){
                            rxSymPySetup(mod.1);
                        } else if (ncmt == 2){
                            rxSymPySetup(mod.2);
                        } else if (ncmt == 3){
                            rxSymPySetup(mod.3);
                        }
                    }
                    cpn <- NULL;
                    if (is.null(cpn)){
                        w <- which(cur.p == curP);
                        if (length(w) == 1){
                            cpn <- w;
                        }
                    }
                    if (oral){
                        CpO <- rxSymPy(sprintf("diff(%s, %s)",
                                               rxToSymPy("CpO"),
                                               rxToSymPy(curP)));
                        CpO <- paste(sub("+=-", "-=", paste0(" ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(CpO))), ";"), fixed=TRUE), collapse="\n");
                        ipp();
                    } else {
                        Cp <- rxSymPy(sprintf("diff(%s, %s)",
                                              rxToSymPy("Cp"),
                                              rxToSymPy(curP)));
                        Cp <- paste(sub("+=-", "-=", paste0("  ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(Cp))), ";"), fixed=TRUE), collapse="\n");
                        ipp();
                        CpIin <- rxSymPy(sprintf("diff(%s, %s)",
                                                 rxToSymPy("CpIin" ),
                                                 rxToSymPy(curP)));
                        CpIin <- paste(sub("+=-", "-=", paste0("  ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(CpIin))), ";"), fixed=TRUE), collapse="\n");
                        ipp();
                        CpIout <- rxSymPy(sprintf("diff(%s, %s)",
                                                  rxToSymPy("CpIout" ),
                                                  rxToSymPy(curP)));
                        CpIout <- paste(sub("+=-", "-=", paste0("  ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(CpIout))), ";"), fixed=TRUE), collapse="\n");
                        ipp();
                    }
                    saveName <- sprintf("d%s", toupper(curP));
                    tmf <- tempfile()
                    if (oral){
                        tmp <- CpO
                    } else {
                        brew::brew(text=template.diff, output=tmf)
                        tmp <- readLines(tmf);
                        unlink(tmf);
                    }
                    lines[length(lines) + 1] <- sprintf("if(ncmt==%s&&%s&&diff1==%s&&diff2==0){//%s",ncmt,
                                                        ifelse(oral == 1, "ka>0", "ka<=0"), cpn, saveName);
                    first <- FALSE;
                    ii <- 0;
                    for (i in seq_along(pars[seq(1, ncmt * 2)])){
                        tmp <- gsub(rex::rex(boundary, cur.p[i], boundary), cur.p.rep[i], tmp, perl=TRUE)
                    }
                    lines <- c(lines, tmp);
                    lines[length(lines) + 1] <- paste0(" return ret;\n}//", saveName);
                    for (curP2 in pars){
                        cpn2 <- NULL;
                        if (is.null(cpn2)){
                            w <- which(cur.p == curP2);
                            if (length(w) == 1){
                                cpn2 <- w;
                            }
                        }
                        ## Schwarz's theorem says this is communicative since
                        ## this has continuous second order partial
                        ## derivatives
                        if (cpn2 >= cpn){
                            if (any(c(curP, curP2) == "tlag")){
                                if (ncmt == 1){
                                    rxSymPySetup(mod.1.tlag);
                                } else if (ncmt == 2){
                                    rxSymPySetup(mod.2.tlag);
                                } else if (ncmt == 3){
                                    rxSymPySetup(mod.3.tlag);
                                }
                            } else {
                                if (ncmt == 1){
                                    rxSymPySetup(mod.1);
                                } else if (ncmt == 2){
                                    rxSymPySetup(mod.2);
                                } else if (ncmt == 3){
                                    rxSymPySetup(mod.3);
                                }
                            }
                            ipp();
                            if (oral){
                                CpO <- rxSymPy(sprintf("diff(diff(%s, %s), %s)",
                                                       rxToSymPy("CpO"),
                                                       rxToSymPy(curP),
                                                       rxToSymPy(curP2)));
                                CpO <- paste(sub("+=-", "-=", paste0(" ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(CpO))), ";"), fixed=TRUE), collapse="\n");
                                ipp();
                            } else {
                                Cp <- rxSymPy(sprintf("diff(diff(%s, %s), %s)",
                                                      rxToSymPy("Cp"),
                                                      rxToSymPy(curP),
                                                      rxToSymPy(curP2)));
                                Cp <- paste(sub("+=-", "-=", paste0("  ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(Cp))), ";"), fixed=TRUE), collapse="\n");
                                ipp();
                                CpIin <- rxSymPy(sprintf("diff(diff(%s, %s), %s)",
                                                         rxToSymPy("CpIin" ),
                                                         rxToSymPy(curP),
                                                         rxToSymPy(curP2)));
                                CpIin <- paste(sub("+=-", "-=", paste0("  ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(CpIin))), ";"), fixed=TRUE), collapse="\n");
                                ipp();
                                CpIout <- rxSymPy(sprintf("diff(diff(%s, %s), %s)",
                                                          rxToSymPy("CpIout"),
                                                          rxToSymPy(curP),
                                                          rxToSymPy(curP2)));
                                CpIout <- paste(sub("+=-", "-=", paste0("  ret+=", gsub(" +", "", rxSplitPlusQ(sympyC(CpIout))), ";"), fixed=TRUE), collapse="\n");

                           }
                            saveName <- sprintf("d%s.d%s", toupper(curP), toupper(curP2));
                            if (oral){
                                tmp <- CpO;
                            } else {
                                tmf <- tempfile()
                                brew::brew(text=template.diff, output=tmf)
                                tmp <- readLines(tmf);
                                unlink(tmf);
                            }

                            lines[length(lines) + 1] <- sprintf("if(ncmt==%s&&%s&&diff1==%s&&diff2==%s){//%s", ncmt,
                                                                ifelse(oral == 1, "ka>0", "ka<=0"), cpn, cpn2, saveName);
                            first <- FALSE;
                            ii <- 0;
                            for (i in seq_along(pars[seq(1, ncmt * 2)])){
                                tmp <- gsub(rex::rex(boundary, cur.p[i], boundary), cur.p.rep[i], tmp, perl=TRUE)
                            }
                            lines <- c(lines, tmp);
                            lines[length(lines) + 1] <- paste0(" return ret;\n}//", saveName);
                        }
                    }
                }
                rxSymPyClean();
                sink(devtools::package_file("src/lincmtDiff.cpp"), append=TRUE);
                cat("\n", paste(lines, collapse="\n"));
                sink();
                cat("done\n");
            }
        }
        sink(devtools::package_file("src/lincmtDiff.cpp"), append=TRUE);
        cat("\n stop(\"Linear derivatives not calculated; Somethings wrong.\");\n");
        cat(" return ret;\n")
        cat("}\n");
        sink();
    }

    create.cpp.code();

    rxClean();
    ## load_all();
    cat("done.\n");
}

document();

