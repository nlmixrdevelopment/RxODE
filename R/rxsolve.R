##' Solves a ODE equation
##'
##' This uses RxODE family of objects, file, or model specification to
##' solve a ODE system.
##'
##' @param object is a either a RxODE family of objects, or a file-name
##'     with a RxODE model specification, or a string with a RxODE
##'     model specification.
##'
##' @param covs_interpolation specifies the interpolation method for
##'     time-varying covariates. When solving ODEs it often samples
##'     times outside the sampling time specified in \code{events}.
##'     When this happens, the time varying covariates are
##'     interpolated.  Currently this can be:
##'
##' \itemize{
##' \item \code{"linear"} interpolation (the default), which interpolates the covariate
##'     by solving the line between the observed covariates and extrapolating the new
##'     covariate value.
##' \item \code{"constant"} -- Last observation carried forward.
##' \item \code{"NOCB"} -- Next Observation Carried Backward.  This is the same method
##'       that NONMEM uses.
##' \item \code{"midpoint"} Last observation carried forward to midpoint; Next observation
##'   carried backward to midpoint.
##' }
##'
##' @param theta A vector of parameters that will be named THETA[#] and
##'     added to parameters
##'
##' @param eta A vector of parameters that will be named ETA[#] and
##'     added to parameters
##'
##' @param add.cov A boolean indicating if covariates should be added
##'     to the output matrix or data frame. By default this is
##'     disabled.
##'
##' @param params a numeric named vector with values for every
##'     parameter in the ODE system; the names must correspond to the
##'     parameter identifiers used in the ODE specification;
##'
##' @param events an \code{eventTable} object describing the input
##'     (e.g., doses) to the dynamic system and observation sampling
##'     time points (see \code{\link{eventTable}});
##'
##' @param inits a vector of initial values of the state variables
##'     (e.g., amounts in each compartment), and the order in this
##'     vector must be the same as the state variables (e.g., PK/PD
##'     compartments);
##'
##' @param scale a numeric named vector with scaling for ode
##'     parameters of the system.  The names must correstond to the
##'     parameter identifiers in the ODE specification. Each of the
##'     ODE variables will be divided by the scaling factor.  For
##'     example \code{scale=(center=2)} will divide the center ODE
##'     variable by 2.
##'
##' @param covs a matrix or dataframe the same number of rows as the
##'     sampling points defined in the events \code{eventTable}.  This
##'     is for time-varying covariates.
##'
##' @param stiff a logical (\code{TRUE} by default) indicating whether
##'     the ODE system is stiff or not.
##'
##'     For stiff ODE sytems (\code{stiff = TRUE}), \code{RxODE} uses the
##'     LSODA (Livermore Solver for Ordinary Differential Equations)
##'     Fortran package, which implements an automatic method switching
##'     for stiff and non-stiff problems along the integration
##'     interval, authored by Hindmarsh and Petzold (2003).
##'
##'     For non-stiff systems (\code{stiff = FALSE}), \code{RxODE} uses
##'     DOP853, an explicit Runge-Kutta method of order 8(5, 3) of
##'     Dormand and Prince as implemented in C by Hairer and Wanner
##'     (1993).
##'
##' @param transit_abs boolean indicating if this is a transit
##'     compartment absorption
##'
##' @param atol a numeric absolute tolerance (1e-08 by default) used
##'     by the ODE solver to determine if a good solution has been
##'     achieved;
##'
##' @param rtol a numeric relative tolerance (1e-06 by default) used
##'     by the ODE solver to determine if a good solution has been
##'     achieved.
##'
##' @param maxsteps maximum number of (internally defined) steps allowed
##'     during one call to the solver. (5000 by default)
##'
##' @param hmin The minimum absolute step size allowed. The default
##'     value is 0.
##'
##' @param hmax The maximum absolute step size allowed.  The default
##'     checks for the maximum difference in times in your sampling and
##'     events, and uses this value.  The value 0 is equivalent to
##'     infinite maximum absolute step size.
##'
##' @param hini The step size to be attempted on the first step. The
##'     default value is determined by the solver (when hini = 0)
##'
##' @param maxordn The maximum order to be allowed for the nonstiff
##'     (Adams) method.  The default is 12.  It can be between 1 and
##'     12.
##'
##' @param maxords The maximum order to be allowed for the stiff (BDF)
##'     method.  The default value is 5.  This can be between 1 and 5.
##'
##' @param matrix A boolean inticating if a matrix should be returned
##'     instead of the RxODE's solved object
##'
##' @param ... Other arguments including scaling factors for each
##'     compartment.  This includes S# = numeric will scale a compartment
##'     # by a dividing the compartment amount by the scale factor,
##'     like NONMEM.
##'
##' @return An \dQuote{rxSolve} solve object that stores the solved
##'     value in a matrix with as many rows as there are sampled time
##'     points and as many columns as system variables (as defined by
##'     the ODEs and additional assignments in the RxODE model code).
##'     It also stores information about the call to allow dynmaic
##'     updating of the solved object.
##'
##'     The operations for the object are simialar to a data-frame, but
##'     expand the \code{$} and \code{[[""]]} access operators and
##'     assignment operators to resolve based on different parameter
##'     values, initial conditions, solver parameters, or events (by
##'     updaing the \code{time} variable).
##'
##'     You can call the \code{\link{eventTable}} methods on the solved
##'     object to update the event table and resolve the system of
##'     equations.  % Should be able to use roxygen templates...
##'
##' @references
##'
##' Hindmarsh, A. C.
##' \emph{ODEPACK, A Systematized Collection of ODE Solvers}.
##' Scientific Computing, R. S. Stepleman et al. (Eds.),
##' North-Holland, Amsterdam, 1983, pp. 55-64.
##'
##' Petzold, L. R.
##' \emph{Automatic Selection of Methods for Solving Stiff and Nonstiff
##' Systems of Ordinary Differential Equations}.
##' Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
##'
##' Hairer, E., Norsett, S. P., and Wanner, G.
##' \emph{Solving ordinary differential equations I, nonstiff problems}.
##' 2nd edition, Springer Series in Computational Mathematics,
##' Springer-Verlag (1993).
##'
##'
##' @seealso \code{\link{RxODE}}
##' @author Melissa Hallow, Wenping Wang and Matthew Fidler
##' @export
rxSolve <- function(object, params = NULL, events = NULL, inits = NULL, covs = NULL, method = "lsoda", transit_abs = NULL, atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0L, maxordn = 12L, maxords = 5L, cores, covs_interpolation = "linear", add.cov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL, sigmaNcores = 1L, sigmaIsChol = FALSE, amountUnits = NA_character_, timeUnits = "hours", stiff, theta = NULL, eta = NULL, update.object=FALSE){
    if (!missing(stiff) && missing(method)){
        if (rxIs(stiff, "logical")){
            if (stiff){
                method <- "lsoda"
                warning("stiff=TRUE has been replaced with method = \"lsoda\".")
            } else {
                method <- "dop853"
                warning("stiff=FALSE has been replaced with method = \"dop853\".")
            }
        }
    }
    if (missing(cores)){
        if (method == "lsoda"){
            cores <- rxCores();
        } else {
            cores <- 1L;
        }
    }
    nms <- names(as.list(match.call())[-1]);
    .Call(`_RxODE_rxSolveC`, object, nms,
          params, events, inits, covs, method, transit_abs, atol, rtol,
          maxsteps, hmin, hmax, hini, maxordn, maxords, cores,
          covs_interpolation, add.cov, matrix, sigma, sigmaDf,
          sigmaNcores, sigmaIsChol, amountUnits, timeUnits, theta, eta, update.object);
}

##' @export
update.rxSolve <- function(object, ...){
    rxSolve(object, ...);
}


##' @author Matthew L.Fidler
##' @export
print.rxSolve <- function(x, ...){
    if (rxIs(x, "rxSolve")){
        args <- as.list(match.call(expand.dots = TRUE));
        if (any(names(args) == "n")){
            n <- args$n;
        } else {
            n <- 6L;
        }
        if (any(names(args) == "width")){
            width <- args$width;
        } else {
            width <- NULL;
        }
        env <- attr(x, ".RxODE.env");
        rxode <- env$object;
        message("Solved RxODE object");
        is.dplyr <- requireNamespace("dplyr", quietly = TRUE) && RxODE.display.tbl;
        ## cat(sprintf("Dll: %s\n\n", rxDll(x)))
        df <- x$params.single
        if (!is.null(df)){
            message("\nParameters ($params):");
            print(df)
        } else {
            df <- x$pars
            if (!is.null(df)){
                message("\nParameters ($params):");
                if (rxIs(df, "data.frame")){
                    if (!is.dplyr){
                        print(head(as.matrix(df), n = n));
                    } else {
                        print(dplyr::as.tbl(df), n = n, width = width);
                    }
                }
            }
        }
        df <- x$covs;
        if (!is.null(df)){
            message("\nCovariates ($covs):");
            if (!is.dplyr){
                print(head(as.matrix(df), n = n));
            } else {
                print(dplyr::as.tbl(df), n = n, width = width);
            }
        }

        message("\nInitial Conditions ($inits):")
        print(x$inits);
        ## inits <- lst$inits[regexpr(regSens, names(lst$inits)) == -1];
        ## print(inits);
        message("\nFirst part of data (object):")
        if (!is.dplyr){
            print(head(as.matrix(x), n = n));
        } else {
            print(dplyr::as.tbl(x), n = n, width = width);
        }
    } else {
        print.data.frame(x)
    }
}

##' @author Matthew L.Fidler
##' @export
summary.rxSolve <- function(object, ...){
    if (rxIs(x, "rxSolve")){
        message("Model:");
        message("################################################################################");
        message(rxNorm(object));
        message("################################################################################");
        message("Parameters:")
        is.dplyr <- requireNamespace("dplyr", quietly = TRUE) && RxODE.display.tbl;
        df <- object$pars
        if (!is.dplyr){
            print(head(as.matrix(x), n = 6L));
        } else {
            print(dplyr::as.tbl(x$pars), n = 6L, width = width);
        }
        message("\n\nInitial conditions:")
        print(object$inits);
        message("\n\nSummary of solved data:")
        print(summary.data.frame(object))
    } else {
        class(x) <- "data.frame"
        NextMethod("print", x);
    }
}

##' @author Matthew L.Fidler
##' @export
is.rxSolve <- function(x){
    .Call(`_RxODE_rxIs`, x, "rxSolve");
}

##' @author Matthew L.Fidler
##' @export
`$.rxSolve` <-  function(obj, arg, exact = FALSE){
    return(.Call(`_RxODE_rxSolveGet`, obj, arg, exact))
}

##' @author Matthew L.Fidler
##' @export
`[.rxSolve` <- function(x, i, j, drop=if (missing(i))
                                              TRUE
                                          else length(cols) == 1){
    class(x) <- "data.frame";
    NextMethod("[");
}

##' @author Matthew L.Fidler
##' @export
"[[.rxSolve" <- function(obj, arg, exact = TRUE){
    return(.Call(`_RxODE_rxSolveGet`, obj, arg, exact))
}

##' @export
t.rxSolve <- function(x){
    x <- as.matrix(x)
    NextMethod("t", x);
}

##' @export
dimnames.rxSolve <- function(x){
    list(row.names(x), names(x));
}

##' @export
"dimnames<-.rxSolve" <- function(x, value){
    class(x) <- "data.frame";
    "dimnames<-.data.frame"(x, value);
}

##'@export
"[<-.rxSolve" <- function(x, i, j, value){
    if (missing(i) && rxIs(j, "character")){
        message("here")
        ret <- .Call(`_RxODE_rxSolveUpdate`, x, j, value);
        if (is.null(ret)){
            class(x) <- "data.frame";
            return(`[<-.data.frame`(x,, j, value = value))
        } else {
            return(ret);
        }
    } else {
        class(x) <- "data.frame"
        if (nargs() < 4){
            if (missing(j)){
                return(`[<-.data.frame`(x, i, value = value))
            } else {
                return(`[<-.data.frame`(x,, j, value = value))
            }
        } else{
            return(`[<-.data.frame`(x, i, j, value))
        }
    }
}
##'@export
`$<-.rxSolve` <- function(x, name, value){
    ret <- .Call(`_RxODE_rxSolveUpdate`, x, name, value);
    if (is.null(ret)){
        class(x) <- "data.frame"
        return (`$<-.data.frame`(x, name, value));
    } else {
        return(ret);
    }
}
##'@export
"[[<-.rxSolve" <- function(x, i, j, value){
    if (missing(j) && rxIs(i, "character")){
        ret <- .Call(`_RxODE_rxSolveUpdate`, x, i, value);
        if (!is.null(ret)){
            return(ret);
        } else {
            class(x) <- "data.frame"
            if (missing(j)){
                return("[[<-.data.frame"(x, i, value=value))
            } else {
                return("[[<-.data.frame"(x, i, j, value))
            }

        }
    } else {
        class(x) <- "data.frame"
        if (missing(j)){
            return("[[<-.data.frame"(x, i, value=value))
        } else {
            return("[[<-.data.frame"(x, i, j, value))
        }
    }
}

##' Update Solved object with '+'
##'
##' @param solved Solved object
##' @param new New information added tothe table.
##' @return new solved object
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
`+.rxSolve` <- function(solved, new){
    if (rxIs(new,"rx.event")){
        return(update(solved, events=new));
    } else {
        return(as.data.frame(solved) + new);
    }
}

## dim (gets you nrow and ncol), t, dimnames
##
## [1] $<-           [             [[<-          [<-           all.equal
## [6] anyDuplicated as.data.frame as.data.table as.list       as.matrix
## [11] coerce        coerce<-      dcast         dim           dimnames
## [16] dimnames<-    duplicated    format        head          initialize
## [21] is.na         melt          merge         na.omit       names<-
## [26] Ops           print         show          slotsFromS3   split
## [31] subset        tail          transform     unique        within
