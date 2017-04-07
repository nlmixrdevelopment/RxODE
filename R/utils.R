rxReq <- function(pkg){
    ## nocov start
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package \"%s\" needed for this function to work. Please install it.", pkg), call. = FALSE);
    }
    ## nocov end
}
asTbl <- function(obj){
    if (getOption("RxODE.prefer.tbl", TRUE) && requireNamespace("dplyr", quietly = TRUE)){
        return(dplyr::as.tbl(as.data.frame(obj)));
    } else {
        return(as.data.frame(obj))
    }
}

rxTbl <- function(x, msg){
    if (getOption("RxODE.prefer.tbl", TRUE) && class(x) == "data.frame" && requireNamespace("dplyr", quietly = TRUE)){
        if (!missing(msg)){
            rxCat(sprintf("Change solved object to dplyr's tbl for %s\n", msg));
        }
        return(dplyr::as.tbl(x))
    } else {
        if (!missing(msg)){
            rxCat(sprintf("Change solved object to data.frame for %s\n", msg))
        }
        return(x)
    }
}
##' Use cat when RxODE.verbose is TRUE
##'
##' @param ... Parameters sent to cat
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxCat <- function(a, ...){
    ## nocov start
    if (getOption("RxODE.verbose", TRUE)){
        if (class(a) == "RxODE"){
            cat(rxNorm(a));
        } else {
            cat(a, ...);
        }
    }
    ## nocov end
}
##' Cleanup anonymous dlls
##'
##' This cleans up any dlls created by text files
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
rxClean <- function(wd = getwd()){
    owd <- getwd();
    setwd(wd);
    on.exit(setwd(owd));
    pat <- "^(Makevars|(rx.*)[.](o|dll|s[ol]|c|rx|prd))$"
    files <- list.files(pattern = pat);
    for (f in files){
        if (f == "Makevars" && file.exists("tran.c")){
            ## nocov start
            warning("Ignoring Makevars since 'tran.c' is in the same directory.")
            ## nocov end
        } else {
            try(dyn.unload(f), silent = TRUE);
            unlink(f);
        }
    }
    return(length(list.files(pattern = pat)) == 0);
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

