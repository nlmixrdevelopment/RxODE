.cliRule <- function(...){
    cat(utils::capture.output(cli::rule(...)), "\n", sep="")
}
.normalizePath <- function(path, ...){
    ifelse(.Platform$OS.type=="windows",
           suppressWarnings(utils::shortPathName(normalizePath(path, ...))),
    ifelse(regexpr("^[/~]", path) != -1,
           suppressWarnings(normalizePath(path, ...)),
           suppressWarnings(normalizePath(file.path(getwd(), path), ...))))
}

##' Require namespace, otherwise throw error.
##'
##' @param pkg Package required for function to work.
##' @return Nothing
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxReq <- function(pkg){
    ## nocov start
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package \"%s\" needed for this function to work. Please install it.", pkg), call. = FALSE);
    }
    ## nocov end
}
##' Use cat when RxODE.verbose is TRUE
##'
##' @param ... Parameters sent to cat
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxCat <- function(a, ...){
    ## nocov start
    if (RxODE.verbose){
        if (is(a, "RxODE")){
            message(RxODE::rxNorm(a), appendLF=FALSE);
        } else {
            message(a, ..., appendLF=FALSE);
        }
    }
    ## nocov end
}



##' Cleanup anonymous DLLs
##'
##' This cleans up any DLLs created by text files
##'
##' @param wd What directory should be cleand
##'
##' This cleans up all files named rx-*.dll and associated files as
##' well as call_dvode.o and associated files
##'
##' @return TRUE if successful
##'
##' @author Matthew L. Fidler
##' @export
rxClean <- function(wd){
    rxTempDir()
    if (missing(wd)){
        ret <- rxClean(getwd()) && rxClean(rxTempDir());
        if (getFromNamespace("RxODE.cache.directory", "RxODE") != "."){
            ret <- ret && rxClean(getFromNamespace("RxODE.cache.directory", "RxODE"))
        }
        return(ret);
    } else if (dir.exists(wd)){
        owd <- getwd();
        setwd(wd);
        on.exit(setwd(owd));
        pat <- "^(Makevars|.*[.]lock|(rx|ui|saem)(.*)[.](o|dll|s[ol]|c|rx|prd|inv|dvdx|rxd|saemd|uid|bad))$"
        files <- list.files(pattern = pat, full.names=TRUE);
        for (f in files){
            if (f == "Makevars"){
                l1 <- readLines("Makevars", n=1L);
                if (l1 == "#RxODE Makevars"){
                    unlink(f);
                }
            } else {
                try(dyn.unload(f), silent = TRUE);
                unlink(f, recursive = TRUE, force=TRUE);
            }
        }
        if (.normalizePath(wd) != .normalizePath(getFromNamespace("RxODE.cache.directory", "RxODE"))){
            ## rxCat("Cleaning cache directory as well.\n");
            rxClean(getFromNamespace("RxODE.cache.directory", "RxODE"));
        }
        .unloadRx()
        return(length(list.files(pattern = pat)) == 0);
    } else {
        return(TRUE)
    }
}

refresh <- function(derivs=FALSE){
    ## nocov start
    cat("Dparser Version\n");
    print(dparser::dpVersion())
    if (derivs){
        Sys.setenv(RxODE_derivs=TRUE)
    } else {
        Sys.setenv(RxODE_derivs=FALSE)
    }
    source(devtools::package_file("build/refresh.R"))
    ## nocov end
}

ode.h <- function(){
    ## nocov start
    cat("Generate header string.\n");
    r.files <- list.files(devtools::package_file("R"), "[.]R$", full.names=TRUE);
    r.files <- r.files[regexpr("RxODE_md5.R", r.files, fixed=TRUE) == -1]
    md5 <- digest::digest(c(sapply(list.files(devtools::package_file("src"),
                                              pattern="[.](c|cpp|h|hpp|f|R|in)$", full.names=TRUE),
                                   function(x){digest::digest(x, file=TRUE)}),
                            sapply(r.files, function(x){digest::digest(x, file=TRUE)}),
                            sapply(list.files(devtools::package_file("vignettes"), pattern="[.](Rmd)$",
                                              full.names=TRUE),
                                   function(x){digest::digest(x, file=TRUE)}),
                            sapply(list.files(devtools::package_file("demo"), pattern="(R|index)$",
                                              full.names=TRUE),
                                   function(x){digest::digest(x, file=TRUE)}),
                            sapply(list.files(devtools::package_file(),
                                              pattern="(cleanup.*|configure.*|DESCRIPTION|NAMESPACE)",
                                              full.names=TRUE),
                                   function(x){digest::digest(x, file=TRUE)})))
    hd <- c(sprintf("#define __VER_2__ \"    SET_STRING_ELT(version,2,mkChar(\\\"%s\\\"));\\n\"", md5),
            sprintf("#define __VER_1__ \"    SET_STRING_ELT(version,1,mkChar(\\\"%s\\\"));\\n\"",
                    as.vector(RxODE::rxVersion()["repo"])),
            sprintf("#define __VER_0__ \"    SET_STRING_ELT(version,0,mkChar(\\\"%s\\\"));\\n\"",
                    sessionInfo()$otherPkgs$RxODE$Version),
            sprintf("#define __VER_md5__ \"%s\"", md5),
            sprintf("#define __VER_repo__ \"%s\"", as.vector(RxODE::rxVersion()["repo"])),
            sprintf("#define __VER_ver__ \"%s\"", sessionInfo()$otherPkgs$RxODE$Version))
    writeLines(hd, devtools::package_file("src/ode.h"))
    writeLines(sprintf("RxODE.md5 <- \"%s\"", md5),
               devtools::package_file("R/RxODE_md5.R"));
    ## nocov end
}



##' Choose the type of sums to use for RxODE.
##'
##' Choose the types of sums to use in RxODE.  These are used in the
##' RxODE \code{sum} blocks and the \code{rxSum} function
##'
##' @param type Sum type to use for \code{rxSum} and \code{sum()} in
##'     RxODE code blocks.
##'
##' \code{pairwise} uses the pairwise sum (fast, default)
##'
##' \code{fsum} uses Python's fsum function (most accurate)
##'
##' \code{kahan} uses kahan correction
##'
##' \code{neumaier} uses Neumaier correction
##'
##' \code{c} uses no correction, bud default/native summing
##'
##' @return nothing
##' @author Matthew L. Fidler
##' @export
rxSetSum <- function(type=c("pairwise", "fsum", "kahan", "neumaier", "c")){
    PreciseSums::psSetSum(type);
}

##' Choose the type of product to use in RxODE.  These are used in the
##' RxODE \code{prod} blocks
##'
##' @param type  Product to use for \code{prod()} in RxODE blocks
##'
##' \code{long double} converts to long double, performs the
##' multiplication and then converts back.
##'
##' \code{double} uses the standard double scale for multiplication.
##'
##' @return nothing
##' @author Matthew L. Fidler
##' @export
rxSetProd <- function(type=c("long double", "double", "logify")){
    PreciseSums::psSetProd(type);
}

##' Setup C++14 support in windows (required for nlmixr)
##'
##' @return nothing
##'
##' @export
rxC14 <- function(){
    ## nocov start
    .dotR <- file.path(Sys.getenv("HOME"), ".R")
    if (!file.exists(.dotR)) dir.create(.dotR)
    .M <- file.path(.dotR, ifelse(.Platform$OS.type!="windows", "Makevars.win", "Makevars"))
    if (!file.exists(.M)) file.create(.M)
    .lines <-suppressWarnings(readLines(.M));
    .write <- FALSE
    .w <- which(regexpr(rex::rex(any_spaces, "CXX14", any_spaces, "="), .lines) != -1)
    if (length(.w)==0L){
        .write <- TRUE
        .lines[length(.lines)+1] <- ifelse(.Platform$OS.type!="windows","CXX14=$(BINPREF)g++ $(M_ARCH)", "CXX14=g++")
    }
    if (.Platform$OS.type!="windows"){
        .w <- which(regexpr(rex::rex(any_spaces, "CXX14STD", any_spaces, "="), .lines) != -1)
        if (length(.w)==0L){
            .write <- TRUE
            .lines[length(.lines)+1] <- "CXX14STD=-std=c++1y"
        }
    }
    .w <- which(regexpr(rex::rex(any_spaces, "CXX14FLAGS", any_spaces, "="), .lines) != -1)
    if (length(.w)==0L){
        .write <- TRUE
        .lines[length(.lines)+1] <- "CXX14FLAGS=-O2 -Wall"
    }
    if (.write) {
        writeLines(.lines, .M);
        message("C++14 setup")
    } else {
        message("C++14 was already setup")
    }
    return(invisible(""));
    ## nocov end
}
