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
    if (RxODE.prefer.tbl && class(x) == "data.frame" && requireNamespace("dplyr", quietly = TRUE)){
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
        if (class(a) == "RxODE"){
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
##' do.call while setting up windows path for compiling
##'
##' @param ... Arguments sent to do.call
##' @param envir Environment to do the call in
##' @return do.call return
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rx.do.call <- function(..., envir=parent.frame()){
    rtools.base <- rxRtoolsBaseWin();
    if (file.exists(rtools.base)){
        old.path <- Sys.getenv("PATH");
        old.binpref <- Sys.getenv("BINPREF");
        on.exit(Sys.setenv(PATH=old.path, BINPREF=old.binpref));
        path <- unique(sapply(gsub("/", "\\\\", strsplit(Sys.getenv("PATH"), ";")[[1]]), function(x){
            if (file.exists(x)){
                return(normalizePath(x));
            } else {
                return("");
            }
        }))
        path <- path[path != ""];
        path <- path[regexpr(rex::rex(or("Rtools", "RTOOLS", "rtools")), path) == -1];
        gcc <- list.files(rtools.base, "gcc",full.names=TRUE)[1]
        if (is.na(gcc)){
            gcc <- "";
        }
        for (x in rev(c(file.path(rtools.base, "bin"),
                        file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                        file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/opt/bin", "mingw_64/opt/bin")),
                        ifelse(gcc == "", "", file.path(gcc, "bin")),
                        ifelse(gcc == "", "", ifelse(.Platform$r_arch == "i386",file.path(gcc, "bin32"), file.path(gcc, "bin64")))))){
            if (file.exists(x)){
                path <- c(normalizePath(x), path);
            }
        }
        path <- paste(unique(path), collapse=";");
        Sys.setenv(PATH=path);
        x <- file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin"));
        if (file.exists(x)){
            Sys.setenv(BINPREF=gsub("([^/])$", "\\1/", gsub("\\\\", "/", normalizePath(x))));
        }
    }
    do.call(..., envir=envir);
}



#' Using the Kahan method, take a more accurate sum
#'
#' @param numbers A vector of numbers to sum.
#' @return Sum of numbers
#' @references
#' \url{https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
#' @export
rxKahanSum <- function(numbers) {
    .Call(`_rxKahanSum`, as.double(numbers))
}

#' Using the Neumaier method, take a more accurate sum
#'
#' @inheritParams rxKahanSum
#' @return Sum of numbers, a bit more accurate than rxKahanSum
#' @references
#' \url{https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
#' @export
rxNeumaierSum <- function(numbers) {
    .Call(`_rxNeumaierSum`, as.double(numbers))
}

#' Return an accurate floating point sum of values
#'
#' This method avoids loss of precision by tracking multiple
#' intermediate partial sums. Based on python's math.fsum
#'
#' @inheritParams rxKahanSum
#'
#' @return Sum of numbers without loss of precision
#'
#' The algorithm's accuracy depends on IEEE-754 arithmetic guarantees
#' and the typical case where the rounding mode is half-even. On some
#' non-Windows builds, the underlying C library uses extended
#' precision addition and may occasionally double-round an
#' intermediate sum causing it to be off in its least significant bit.
#'
#' @export
#' @author Matthew Fidler (R implementation), Raymond Hettinger,
#'     Jonathan Shewchuk, Python Team
#' @references
#'
#' \url{https://docs.python.org/2/library/math.html}
#'
#' \url{https://code.activestate.com/recipes/393090/}
#'
#' \url{https://github.com/python/cpython/blob/a0ce375e10b50f7606cb86b072fed7d8cd574fe7/Modules/mathmodule.c}
#'
#' Shewchuk, JR. (1996)
#' \emph{Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.}
#' \url{http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps}
#'
rxPythonFsum <- function(numbers) {
    .Call(`_rxPythonSum`, as.double(numbers))
}

#' Return an accurate floating point sum of values
#'
#' This was taken by NumPy and adapted for use here.
#'
#' @inheritParams rxKahanSum
#'
#' @return A sum of numbers with a rounding error of O(lg n) instead
#'     of O(n).
#' @author Matthew Fidler (R implementation), Julian Taylor, Nathaniel
#'     J Smith, and others in NumPy team.
#' @references
#' \url{https://github.com/juliantaylor/numpy/blob/b0bc01275cac04483e6df021211c1fa2ba65eaa3/numpy/core/src/umath/loops.c.src}
#'
#' \url{https://github.com/numpy/numpy/pull/3685}
#'
rxPairwiseSum <- function(numbers) {
    .Call(`_rxPairwiseSum`, as.double(numbers))
}

#' Using RxODE's default method, take a sum
#'
#' @inheritParams rxKahanSum
#' @return Sum of numbers
#' @export
rxSum <- function(numbers) {
    .Call(`_rxSum`, as.double(numbers))
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
    i <- which(type == c("pairwise", "fsum", "kahan", "neumaier", "c"))
    invisible(.Call(`_rxSetSum`, as.integer(i)))
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
    i <- which(type == c("long double", "double", "logify"));
    invisible(.Call(`_rxSetProd`, as.integer(i)))
}


#' Using RxODE's default method, take a product
#'
#' @inheritParams rxKahanSum
#' @return Product of numbers
#' @export
rxProd <- function(numbers) {
    .Call(`_rxProd`, as.double(numbers))
}
