.onLoad <- function(libname, pkgname){ ## nocov start
    ## Setup RxODE.prefer.tbl
    rxPermissive(respect=TRUE); ## need to call respect on the first time
    ## memoise needs to be called at load to use the right package.
    ## See https://github.com/hadley/r-pkgs/issues/203
    ## They suggest an environment, but I used the current namespace.
    utils::assignInMyNamespace("rxModelVars.character", memoise::memoise(rxModelVars.character.slow));
    utils::assignInMyNamespace("rxExpandIfElse", memoise::memoise(rxExpandIfElse.slow));
    utils::assignInMyNamespace("rxSymPyDfDyFull", memoise::memoise(rxSymPyDfDyFull.slow));
    utils::assignInMyNamespace("rxSymPyJacobian", memoise::memoise(rxSymPyJacobian.slow));
    utils::assignInMyNamespace("rxSymPySensitivity", memoise::memoise(rxSymPySensitivity.slow));
} ## nocov end

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
              RxODE.suppress.syntax.info     =c(FALSE, FALSE)
              );

##' Permissive or Strict RxODE sytax options
##'
##' This sets the RxODE syntax to be permissive or strict
##'
##' @param expr Expression to evaluate in the permissive/strict
##'     environment.  If unspecified, set the options for the current
##'     environment.
##' @param silent when true, also silence the syntax errors and
##'     interactive output (useful in testing).
##' @param respect when TRUE, respect any options that are specified.
##'     This is called at startup, but really should not be called
##'     elsewhere, otherwise the options are not changed.
##' @param rxclean when TRUE, call rxClean before and after the expr
##'     is called.
##'
##' @author Matthew L. Fidler
##' @export
rxPermissive <- function(expr, silent=FALSE, respect=FALSE, rxclean=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != -1)){
    args  <- as.list(match.call())[-1];
    args$op.rx <- 2;
    do.call("rxOptions", args, envir=parent.frame(1));
}
##' @rdname rxPermissive
##' @export
rxStrict <- function(expr, silent=FALSE, respect=FALSE, rxclean=(regexpr("/tests/testthat/", getwd(), fixed=TRUE) != -1)){
    args  <- as.list(match.call())[-1];
    args$op.rx <- 1;
    do.call("rxOptions", args, envir=parent.frame(1));
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
rxOptions <- function(expr, op.rx=NULL, silent=FALSE, respect=FALSE,
                      rxclean=(regexpr("/tests/testthat/", getwd(), fixed=TRUE))){
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
            on.exit({options(opOld); if (rxclean){rxClean();}});
        }
        if (respect){
            op <- options();
            w <- !(names(op.rx) %in% names(op))
            if (any(w)) options(op.rx[w]);
        } else {
            options(op.rx);
        }
        if (class(substitute(expr)) == "{"){
            return(eval(substitute(expr), envir=parent.frame(1)));
        }
    }
}
