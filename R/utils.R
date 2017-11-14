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
asTbl <- function(obj){
    if (RxODE.prefer.tbl && requireNamespace("dplyr", quietly = TRUE)){
        return(dplyr::as.tbl(as.data.frame(obj)));
    } else {
        return(as.data.frame(obj))
    }
}

rxTbl <- function(x, msg){
    if (RxODE.prefer.tbl && is(x,"data.frame") && requireNamespace("dplyr", quietly = TRUE)){
        if (!missing(msg)){
            rxCat(sprintf("Change solved object to dplyr's tbl for %s.\n", msg));
        }
        return(dplyr::as.tbl(x))
    } else {
        if (!missing(msg)){
            rxCat(sprintf("Change solved object to data.frame for %s.\n", msg))
        }
        attr(x, ".env") <- NULL;
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
    if (RxODE.verbose){
        if (is(a,"RxODE")){
            message(rxNorm(a), appendLF=FALSE);
        } else {
            message(a, ..., appendLF=FALSE);
        }
    }
    ## nocov end
}

##' Print x using the message facility
##'
##' This allows the suppressMessages to work on print functions.  This
##' captures the output via R.Util's captureOutput function and then
##' sends it through the message routine.
##'
##' catpureOutput was used since it is much faster than the internal
##' capture.output see https://www.r-bloggers.com/performance-captureoutput-is-much-faster-than-capture-output/
##' @param x object to print
##' @param ... Other things output
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxPrint <- function(x, ...){
    message(invisible(paste(R.utils::captureOutput(x <<- print(x, ...)), collapse="\n")), appendLF=TRUE);
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
rxClean <- function(wd = getwd()){
    owd <- getwd();
    setwd(wd);
    on.exit(setwd(owd));
    pat <- "^(Makevars|(rx.*)[.](o|dll|s[ol]|c|rx|prd|inv))$"
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
    if (normalizePath(wd) != normalizePath(RxODE.cache.directory)){
        ## rxCat("Cleaning cache directory as well.\n");
        rxClean(RxODE.cache.directory);
    }
    return(length(list.files(pattern = pat)) == 0);
}

##' Give the cout string for a number
##'
##' @param number Number for converting
##' @return String of cout << number;
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxCout <- function(number){
    tmp <- tempfile();
    orig.sink.number <- sink.number()
    on.exit({while(sink.number() > orig.sink.number){
        sink();
    }
        unlink(tmp)}
    );
    sink(tmp);
    rxCoutEcho(number);
    sink();
    lines <- suppressWarnings({readLines(tmp)});
    return(lines)
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
