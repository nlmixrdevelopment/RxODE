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
        stop(sprintf("package \"%s\" needed for this function to work", pkg), call. = FALSE);
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



##' Cleanup anonymous DLLs by unloading them
##'
##' This cleans up any RxODE loaded DLLs
##'
##' @param wd What directory should be cleaned; (DEPRECIATED), this no
##'     longer does anything.
##'
##' This unloads all RxODE anonymous dlls.
##'
##' @return TRUE if successful
##'
##' @author Matthew L. Fidler
##' @export
rxClean <- function(wd){
    if (!missing(wd)) warning("'wd' is depreciated")
    rxUnloadAll();
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

##' Set timing for progress bar
##'
##' @param seconds This sets the number of seconds that need to elapse
##'     before drawing the next segment of the progress bar.  When
##'     this is zero or below this turns off the progress bar.
##'
##' @export
##' @author Matthew Fidler
rxSetProgressBar <- function(seconds=1.0) {
    invisible(.Call(`_rxParProgress`, as.double(seconds)))
}

##' Sample a covariance Matrix from the Posterior Inverse Wishart
##' distribution.
##'
##' Note this Inverse wishart rescaled to match the original scale of
##' the covariance matrix.
##'
##' If your covariance matrix is a 1x1 matrix, this uses an scaled
##' inverse chi-squared which is equivalent to the Inverse Wishart
##' distribution in the uni-directional case.
##'
##' @param nu Degrees of Freedom (Number of Observations) for
##'        covariance matrix simulation.
##'
##' @param omega Either the estimate of covariance matrix or the
##'     estimated standard deviations in matrix form each row forming
##'     the standard deviation simulated values
##'
##' @param n Number of Matrices to sample.  By default this is 1.
##'     This is only useful when \code{omega} is a matrix.  Otherwise
##'     it is determined by the number of rows in the input
##'     \code{omega} matrix of standard deviations
##'
##' @param omegaIsChol is an indicator of if the omega matrix is in
##'   the Cholesky decomposition. This is only used when code{type="invWishart"}
##'
##' @param returnChol Return the Cholesky decomposition of the
##'   covariance matrix sample. This is only used when code{type="invWishart"}
##'
##' @return a matrix (n=1) or a list of matrices  (n > 1)
##'
##' @author Matthew L.Fidler & Wenping Wang
##'
##' @examples
##'
##' ## Sample a single covariance.
##' draw1 <- cvPost(3, matrix(c(1,.3,.3,1),2,2))
##'
##' ## Sample 3 covariances
##' set.seed(42)
##' draw3 <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3)
##'
##' ## Sample 3 covariances, but return the cholesky decomposition
##' set.seed(42)
##' draw3c <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3, returnChol=TRUE)
##'
##' ## Sample 3 covariances with lognormal standard deviations via LKJ correlation sample
##' cvPost(3,sapply(1:3,function(...){rnorm(10)}), type="lkj")
##'
##' ## Sample 3 covariances with lognormal standard deviations via separation
##' ## strategy using inverse Wishart correlation sample
##' cvPost(3,sapply(1:3,function(...){rnorm(10)}), type="separation")
##'
##' @export
cvPost <- function(nu, omega, n = 1L, omegaIsChol = FALSE, returnChol = FALSE,
                   type = c("invWishart", "lkj", "separation"),
                   diagXformType = c("log", "identity", "nlmixrSqrt", "nlmixrLog", "nlmixrIdentity")){
    if (inherits(type, "numeric") || inherits(type, "integer")){
        .type <- as.integer(type)
    } else {
        .type <- as.vector(c("invWishart"=1L, "lkj"=2L, "separation"=3L)[match.arg(type)]);
    }
    if (.type == 1L){
        .xform <- 1L
    }  else if (inherits(diagXformType, "numeric") || inherits(diagXformType, "integer")){
        .xform <- as.integer(diagXformType)
    } else {
        .xform <- as.vector(c("log"=5, "identity"=4, "nlmixrSqrt"=1, "nlmixrLog"=2, "nlmixrIdentity"=3)[match.arg(diagXformType)])
    }
    return(.Call(`_RxODE_cvPost_`, nu, omega, n, omegaIsChol, returnChol, .type, .xform))
}
