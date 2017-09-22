.onAttach <- function(libname, pkgname){ ## nocov start
    ## Setup RxODE.prefer.tbl
    rxPermissive(respect=TRUE); ## need to call respect on the first time
    ## memoise needs to be called at load to use the right package.
    ## See https://github.com/hadley/r-pkgs/issues/203
    ## They suggest an environment, but I used the current namespace.
    rxSetupMemoize()
    ## Setup the path
    tmp <- try(rxWinRtoolsPath(), silent=TRUE);
    if (!rxWinRtoolsPath()){
        packageStartupMessage("Rtools is not set up correctly!\n\nYou need a working Rtools installation for RxODE to work.\nYou can set up Rtools using the command 'rxWinSetup()'.\nThis will also set up Python and SymPy to run a bit faster than rSymPy.\n");
    }
} ## nocov end

##' This setups the memoized functions.
##'
##' To easily create a memozied function by adding a \code{.slow <- NULL}
##' to the end of a function.
##'
##' For example, to memozie the function in the namespace
##' \code{rxModelVars.character} you would add a line:
##' \code{rxModelVars.character.slow <- NULL}
##'
##' @author Matthew L. Fidler
rxSetupMemoize <- function(){
    reSlow <- rex::rex(".slow",end)
    f <- sys.function(-1)
    ns <- environment(f)
    .slow <- ls(pattern=reSlow,envir=ns);
    for (slow in .slow){
        fast <- sub(reSlow, "", slow);
        if (!memoise::is.memoised(get(fast, envir=ns)) && is.null(get(slow, envir=ns))){
            utils::assignInMyNamespace(slow, get(fast, envir=ns))
            utils::assignInMyNamespace(fast, memoise::memoise(get(slow, envir=ns)))
        }

    }
}

##' Clear memoise cache for RxODE
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxForget <- function(){
    reSlow <- rex::rex(".slow",end)
    f <- sys.function(-1)
    ns <- environment(f)
    .slow <- ls(pattern=reSlow,envir=ns);
    for (slow in .slow){
        fast <- sub(reSlow, "", slow);
        memoise::forget(get(fast, envir=ns));
    }
}

## strict/permissive
rxOpt <- list(RxODE.prefer.tbl               =c(FALSE, FALSE),
              RxODE.display.tbl              =c(TRUE, TRUE),
              RxODE.echo.compile             =c(FALSE, FALSE),
              RxODE.warn.on.assign           =c(TRUE, TRUE),
              RxODE.compile.on.load          =c(TRUE, TRUE),
              RxODE.syntax.assign            =c(FALSE, TRUE),
              RxODE.syntax.star.pow          =c(FALSE, TRUE),
              RxODE.syntax.require.semicolon =c(TRUE, FALSE),
              RxODE.syntax.allow.dots        =c(FALSE, TRUE),
              RxODE.syntax.allow.ini0        =c(FALSE, TRUE),
              RxODE.syntax.allow.ini         =c(FALSE, TRUE),
              RxODE.calculate.jacobian       =c(FALSE, FALSE),
              RxODE.calculate.sensitivity    =c(FALSE, FALSE),
              RxODE.verbose                  =c(TRUE, TRUE),
              RxODE.suppress.syntax.info     =c(FALSE, FALSE),
              RxODE.sympy.engine             =c("", ""),
              RxODE.cache.directory          =c(".", "."),
              RxODE.delete.unnamed           =c(FALSE, FALSE)
              );

RxODE.prefer.tbl <- NULL
RxODE.display.tbl <- NULL
RxODE.echo.compile <- NULL
RxODE.warn.on.assign <- NULL
RxODE.compile.on.load <- NULL
RxODE.syntax.assign <- NULL
RxODE.syntax.star.pow <- NULL
RxODE.syntax.require.semicolon <- NULL
RxODE.syntax.allow.dots <- NULL
RxODE.syntax.allow.ini0 <- NULL
RxODE.syntax.allow.ini <- NULL
RxODE.calculate.jacobian <- NULL
RxODE.calculate.sensitivity <- NULL
RxODE.verbose <- NULL
RxODE.suppress.syntax.info <- NULL
RxODE.sympy.engine <- NULL
RxODE.cache.directory <- NULL
RxODE.delete.unnamed <- NULL


##' Permissive or Strict RxODE sytax options
##'
##' This sets the RxODE syntax to be permissive or strict
##'
##' @param expr Expression to evaluate in the permissive/strict
##'     environment.  If unspecified, set the options for the current
##'     environment.
##' @param respect when TRUE, respect any options that are specified.
##'     This is called at startup, but really should not be called
##'     elsewhere, otherwise the options are not changed.
##' @param cran When specifyed and true, run on CRAN. Otherwise it is skipped on cran.
##' @param on.validate When TRUE run only when validating.
##' @param silent when true, also silence the syntax errors and
##'     interactive output (useful in testing).
##' @param rxclean when TRUE, call rxClean before and after the expr
##'     is called.
##' @author Matthew L. Fidler
##' @export
rxPermissive <- function(expr, silent=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != -1),
                         respect=FALSE,
                         rxclean=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != -1),
                         cran=FALSE, on.validate=FALSE){
    args  <- as.list(match.call())[-1];
    args$op.rx <- 2;
    do.call(getFromNamespace("rxOptions", "RxODE"), args, envir=parent.frame(1));
}
##' @rdname rxPermissive
##' @export
rxStrict <- function(expr, silent=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != -1), respect=FALSE,
                     rxclean=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != -1),
                     cran=FALSE, on.validate=FALSE){
    args  <- as.list(match.call())[-1];
    args$op.rx <- 1;
    do.call(getFromNamespace("rxOptions", "RxODE"), args, envir=parent.frame(1));
}
##' Options for RxODE
##'
##' This is a backend for \code{rxPermissive} (with
##' \code{op.rx} = \code{2}) and \code{rxStrict} (with
##' \code{op.rx} =\code{1})
##'
##' When \code{expr} is missing and \code{op.rx} is NULL, this
##' desplays the current RxODE options.
##'
##' @inheritParams rxPermissive
##' @param op.rx A numeric for strict (1) or permissive (2) syntax.
##' @author Matthew L. Fidler
##' @export
rxOptions <- function(expr, op.rx=NULL, silent=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != 1), respect=FALSE,
                      rxclean=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != 1),
                      cran=FALSE, on.validate=FALSE){
    do.it <- TRUE
    if (!identical(Sys.getenv("NOT_CRAN"), "true") && !cran){
        ## on Cran, but only tested when not on cran, skip.
        do.it <- FALSE
    }
    if (on.validate && !identical(Sys.getenv("RxODE_VALIDATION_FULL"), "true")){
        do.it <- FALSE
    }
    if (!on.validate && identical(Sys.getenv("RxODE_VALIDATION_FULL"), "true")){
        do.it <- FALSE
    }
    if (do.it){
        if (missing(expr) && is.null(op.rx)){
            op <- options()
            op <- op[regexpr(rex::rex("RxODE."), names(op)) != -1];
            op <- op[order(names(op))];
            sapply(names(op), function(n){rxCat(sprintf("%s: %s\n", n, op[[n]]))});
            return(invisible(op));
        } else {
            if (class(op.rx) == "character"){
                if (op.rx == "strict"){
                    op.rx  <- 1;
                } else {
                    op.rx <- 2;
                }
            }
            if (class(op.rx) == "numeric"){
                if (op.rx <= 2){
                    x  <- op.rx;
                    op.rx  <- list()
                    for (v in names(rxOpt)){
                        op.rx[[v]] <- rxOpt[[v]][x];
                    }
                }
            }
            if (!missing(silent)){
                op.rx$RxODE.verbose=!silent;
                op.rx$RxODE.suppress.syntax.info=silent;
            }
            if (!missing(expr)){
                if (rxclean){
                    rxClean();
                }
                opOld <- options();
                on.exit({options(opOld); rxSyncOptions(); if (rxclean){rxClean();}});
            }
            if (respect){
                op <- options();
                w <- !(names(op.rx) %in% names(op))
                if (any(w)) options(op.rx[w]);
                rxSyncOptions()
            } else if (length(op.rx) > 0){
                options(op.rx);
                rxSyncOptions()
            }
            if (class(substitute(expr)) == "{"){
                if (silent){
                    return(suppressMessages(eval(substitute(expr), envir=parent.frame(1))));
                } else {
                    return(eval(substitute(expr), envir=parent.frame(1)));
                }
            }
        }
    }
}

##' Sync options with RxODE varaibles
##'
##' Accessing RxODE options via getOption slows down solving.  This
##' allows the options to be synced with variables.
##'
##' @author Matthew L. Fidler
##' @export
rxSyncOptions <- function(){
    for (var in names(rxOpt)){
        assignInMyNamespace(var, getOption(var, rxOpt[[var]][1]));
    }
}

rxSkipValidate <- function(){
    if(!identical(Sys.getenv("RxODE_VALIDATION_FULL"), "true")){
        testthat::skip("Only run on full validation")
    }
}

rxSyncOptions();
