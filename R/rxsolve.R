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
##'     added to inits
##'
##' @param eta A vector of parameters that will be named ETA[#] and
##'     added to inits
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
rxSolve <- function(object, params = NULL, events = NULL, inits = NULL, covs = NULL, method = "lsoda", transit_abs = NULL, atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0L, maxordn = 12L, maxords = 5L, cores = 1L, covs_interpolation = "linear", add.cov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL, sigmaNcores = 1L, sigmaIsChol = FALSE, amountUnits = NA_character_, timeUnits = "hours", stiff){
    if (!missing(stiff) && missing(method)){
        if (rxIs(stiff, "logical")){
            if (stiff){
                method <- "lsoda"
                warning("stiff=TRUE has been replaced with method = \"lsoda\".")
            } else {
                method <- "dop"
                warning("stiff=FALSE has been replaced with method = \"dop\".")
            }
        }
    }
    .Call(`_RxODE_rxSolveC`, object, params, events, inits, covs, method, transit_abs, atol, rtol, maxsteps, hmin, hmax, hini, maxordn, maxords, cores, covs_interpolation, add.cov, matrix, sigma, sigmaDf, sigmaNcores, sigmaIsChol, amountUnits, timeUnits);
}


##' @author Matthew L.Fidler
##' @export
print.solveRxODE7 <- function(x, ...){
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
    message("Parameters ($params):");
    df <- x$params.single
    if (!is.null(df)){
        print(df)
    } else {
        df <- x$pars
        if (rxIs(df, "data.frame")){
            if (!is.dplyr){
                print(head(as.matrix(x), n = n));
            } else {
                print(dplyr::as.tbl(x$pars), n = n, width = width);
            }
        }
    }
    ## w <- which((names(lst$params) %in% names(x)))
    ## if (length(w) > 0){
    ##     print(lst$params[-w]);
    ##     message("\nFirst Part of Time Varying Covariates:");
    ##     d <- as.data.frame(env$covs)[, names(lst$params)[w]];
    ##     if (length(w) == 1){
    ##         d <- data.frame(d = d);
    ##         names(d) <- names(lst$params)[w];
    ##     }
    ##     if (!is.dplyr){
    ##         print(head(d), n = n);
    ##     } else {
    ##         print(dplyr::as.tbl(d), n = n, width = width);
    ##     }
    ## }  else {
    ##     p <- lst$params;
    ##     if (length(env$pcov) > 0){
    ##         p2 <- p[-env$pcov];
    ##         print(p2)
    ##         message("\nTime Varying Covariates");
    ##         message(paste(names(p[env$pcov]), collapse=" "));
    ##     } else {
    ##         print(p);
    ##     }
    ## }
    message("\n\nInitial Conditions ($inits):")
    print(x$inits);
    ## inits <- lst$inits[regexpr(regSens, names(lst$inits)) == -1];
    ## print(inits);
    message("\n\nFirst part of data (object):")
    if (!is.dplyr){
        print(head(as.matrix(x), n = n));
    } else {
        print(dplyr::as.tbl(x), n = n, width = width);
    }
}

##' @author Matthew L.Fidler
##' @export
summary.solveRxODE7 <- function(object, ...){
    env <- attr(object, ".env");
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
}

##' @author Matthew L.Fidler
##' @export
`$.solveRxODE7` <-  function(obj, arg, exact = TRUE){
    .Call(`_RxODE_rxSolveGet`, obj, arg);
}

##' @author Matthew L.Fidler
##' @export
`[.solveRxODE` <- function(x, i, j, drop){
    df <- as.data.frame(x);
    assign.names <- NULL
    if (!missing(j) && is(j,"character")){
        nm <- names(df);
        nms <- gsub(regSens, "_sens_\\1_\\2", nm);
        assign.names <- j;
        j <- as.vector(sapply(j, function(x){
            w <- which(nm == x)
            if (length(w) >= 1){
                return(w[1]);
            } else {
                w <- which(nms == x)
                return(w[1]);
            }
        }));
    }
    if (!missing(i) && missing(j) && missing(drop)){
        df <- df[i, ];
    } else if (missing(i) && !missing(j) && missing(drop)){
        df <- df[, j];
        if (!is.null(assign.names) && is(df,"data.frame"))
            names(df) <- assign.names
    } else if (!missing(i) && !missing(j) && missing(drop)){
        df <- df[i, j];
        if (!is.null(assign.names) && is(df,"data.frame"))
            names(df) <- assign.names
    } else if (missing(i) && missing(j) && missing(drop)){
        df <- df[drop = drop];
    } else if (!missing(i) && missing(j) && !missing(drop)){
        df <- df[i,, drop = drop];
    } else if (missing(i) && !missing(j) && !missing(drop)){
        df <- df[, j, drop = drop];
        if (!is.null(assign.names))
            names(df) <- assign.names
    } else if (!missing(i) && !missing(j) && !missing(drop)){
        df <- df[i, j, drop = drop];
        if (!is.null(assign.names) && is(df,"data.frame"))
            names(df) <- assign.names
    } else if (missing(i) && missing(j) && !missing(drop)){
        df <- df[,, drop = drop];
    }
    return(rxTbl(df))
}

##' @author Matthew L.Fidler
##' @export
"[[.solveRxODE" <- function(obj, arg, exact = TRUE, internal = FALSE){
    if (internal){
        env <- attr(obj, ".env");
        tmp <- c(list(params=env$params, inits=env$inits), env$extra.args)
        return(tmp[[arg, exact = exact]]);
    } else {
        `$.solveRxODE`(obj, arg, exact = exact);
    }
}





##' Assign solved objects using the [] syntax
##' @param obj solved object
##' @param arg1 first argument
##' @param arg2 second argument
##' @param value value assumed
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"[<-.solveRxODE" <- function(obj, arg1, arg2, value){
    if (any(arg2 == c(1, "t", "time")) && missing(arg1)){
        obj$time <- value
        return(obj);
    }
    df <- as.data.frame(obj);
    if (missing(arg1) & missing(arg2)){
        df[] <- value;
    } else if (missing(arg1)){
        df[, arg2] <-  value;
    } else if (missing(arg2)){
        df[arg1, ] <-  value;
    } else {
        d3[arg1, arg2] <- value;
    }
    return(rxTbl(df, "assignment"))
}


########################################################
## Updating solved object


##' Update the solved object with any of the new parameters.
##'
##' @param object Object to be updated
##' @param ... Arguments to be updated, and resolved.
##'
##' @author Matthew L.Fidler
##' @export
update.solveRxODE <- function(object, ...){
    env <- attr(object, ".env");
    rxode <- env$env$out;
    args <- c(list(..., params=env$params, inits=env$inits, matrix=FALSE), env$extra.args);
    args <- args[!duplicated(names(args))];
    do.call(rxode$solve, args)
}

##' Update Solved object with '+'
##'
##' @param solved Solved object
##' @param new New information added tothe table.
##' @return new solved object
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
`+.solveRxODE` <- function(solved, new){
    if (is(new,"EventTable")){
        return(update(solved, events=new));
    } else {
        return(as.data.frame(solved) + new);
    }
}

##' $ Assign for RxODE solved objects
##'
##' Assign objects by argumnt as obj$arg <- value
##'
##' This also works as obj[[arg]] <- value
##'
##' @param obj solveRxDll object
##' @param arg Dollar sign name
##' @param value assigned Value.
##' @seealso \code{\link{rxSolve}}
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"$<-.solveRxODE" <- function(obj, arg, value){
    if (arg == "t"){
        arg <- "time";
    }
    env <- attr(obj, ".env")
    rxode <- env$env$out;
    lst <- c(list(params=env$params, inits=env$inits), env$extra.args)
    iarg <- gsub(regIni, "", arg);
    if (arg == "time"){
        if (is(value,"EventTable")){
            cat("Update event table and solved object.\n");
            return(update(obj, events = eventTable))
        } else if (is(value,"data.frame")){
        } else if (is(value,"numeric")){
            rxCat("Updating sampling times in the event table updating object.\n");
            eventTable <- lst$events$copy();
            eventTable$clear.sampling();
            et <- eventTable()
            et$import.EventTable(eventTable$get.EventTable());
            et$add.sampling(value);
            return(update(obj, events = et));
        }
    } else if (arg != iarg && any(rxState(rxode) == iarg) && length(value) == 1){
        rxCat("Updating object with new initial conditions.\n")
        inits <- c(value);
        names(inits) <- gsub(regIni, "", arg);
        return(update(obj, inits = inits));
    } else if (any(rxParams(rxode) == arg)){
        rxCat("Updating object with new parameter values.\n")
        if (length(value) == 1){
            covs <- as.data.frame(lst$covs);
            if (any(names(covs) == arg)){
                cat(sprintf("Changing time-varying covariate \"%s\" to a simple parameter value %s.\n", arg, value));
                ncovs <- names(covs);
                covs <- as.data.frame(covs[, names(covs) != arg]);
                names(covs) <- ncovs[ncovs != arg];
                params <- c(value);
                names(params) <- arg;
                params <- c(params, lst$params);
                return(update(obj, params = params, covs = covs));
            } else {
                params <- c(value);
                names(params) <- arg;
                params <- c(params, lst$params);
                return(update(obj, params = params));
            }
        } else if (length(value) == length(lst$events$get.sampling()$time)){
            cat(sprintf("Changing simple parameter \"%s\" to a time-varying covariate.\n", arg));
            covs <- as.data.frame(lst$covs);
            covs[[arg]] <- value;
            return(update(obj, covs = covs));
        }
    } else if (arg == "params"){
        rxCat("Updating object with new parameter values.\n");
        return(update(obj, params = value));
    } else if (arg == "inits"){
        rxCat("Updating object with new initial conditions.\n")
        return(update(obj, inits = value));
    } else if (any(arg == names(lst))){
        args <- list(value);
        names(args) <- arg;
        rxCat(sprintf("Updating object with new solving argument %s = %s.\n", arg, value))
        return(do.call("update", args, envir = parent.frame(1)));
    } else {
        df <- as.data.frame(obj);
        df <- "$<-.data.frame"(df, arg, value);
        obj <- rxTbl(df, "assignment");
    }
    return(obj);
}


##' Assign solved objects using the [[]] syntax
##' @param obj solved object
##' @param arg element of solved object
##' @param value value assumed
##' @seealso \code{\link{rxSolve}}
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"[[<-.solveRxODE" <- function(obj, arg, value){
    return("$<-.solveRxODE"(obj, arg, value = value))
}


## This causes the whole package to crash...
## ##' @name as.data.table
## ##' @export as.data.table.solveRxDll
## ##'
## ##' @method as.data.table solveRxDll
## ##'
## ##' @title as.data.table for \code{solveRxDll} object
## ##' @description compatability function for tidyr
## ##' @param data Solved ODE, an \code{solveRxDll} object.
## ##' @param ... Additional arguments
## ##'
## as.data.table.solveRxDll <- function(x, ...){
##     call <- as.list(match.call(expand.dots=TRUE))[-1];
##     call$x <- as.data.table(x);
##     return(do.call(getFromNamespace("as.data.table","data.table"), call, envir = parent.frame(1)));
## }

##' Setup Data and Parameters
##'
##' @param dll
##' @inheritParams rxSolve
##' @param sigma Named sigma matrix.
##' @param sigmaDf The degrees of freedom of a t-distribution for
##'     simulation.  By default this is \code{Inf}, or to simulate
##'     from a normal distribution instead of a t.
##' @param sigmaNcores Number of cores for residual simulation.
##'     This, along with the seed, affects both the outcome and speed
##'     of simulation. By default it is one.
##' @param sigmaIsChol Indicates if the \code{sigma} supplied is a
##'     Cholesky decomposed matrix instead of the traditional
##'     symmetric matrix.
##' @return Data setup for running C-based RxODE runs.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export

