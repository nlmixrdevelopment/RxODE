##' @author Matthew L. Fidler
##' @export
##' @keywords internal
as.data.frame.solveRxODE <- function(x, row.names = NULL, optional = FALSE, ...){
    attr(x, ".env") <- NULL;
    class(x) <- "data.frame";
    return(as.data.frame(x, row.names=row.names, optional=optional, ...));
}
##' @name as_data_frame
##' @export as_data_frame.solveRxODE
##' @method as_data_frame solveRxODE
##' @description compatability function for as_data_frame
##' @title as_data_frame for solveRxODE
##' @param x Solved data frame
##' @param ... Additional arguments
##' @author Matthew L. Fidler
##' @keywords internal
as_data_frame.solveRxODE <- function(x, ...){
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    call$x <- as.data.frame(x);
    return(do.call(getFromNamespace("as_data_frame","tibble"), call, envir=parent.frame(1)));
}

##' @name as_data_frame
##' @export as_data_frame.solveRxODE
##' @method as_data_frame solveRxODE
##' @description compatability function for as_data_frame
##' @title as_data_frame for solveRxODE
##' @param x Solved data frame
##' @param ... Additional arguments
##' @author Matthew L. Fidler
##' @keywords internal
as.tbl.solveRxODE <- function(x, ...){
    call <- as.list(match.call(expand.dots=TRUE))[-1];
    call$x <- as.data.frame(x);
    return(do.call(getFromNamespace("as.tbl","dplyr"), call, envir=parent.frame(1)));
}

##' @author Matthew L.Fidler
##' @export
print.solveRxODE <- function(x, ...){
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
    env <- attr(x, ".env");
    rxode <- env$env$out;
    lst <- c(list(params=env$params, inits=env$inits), env$extra.args)
    cat("Solved RxODE object\n");
    cat(sprintf("Dll: %s\n\n", rxDll(rxode)))
    cat("Parameters:\n")
    is.dplyr <- requireNamespace("dplyr", quietly = TRUE) && RxODE.display.tbl;
    w <- which((names(lst$params) %in% names(x)))
    if (length(w) > 0){
        print(lst$params[-w]);
        message("\nFirst Part of Time Varying Covariates:");
        d <- as.data.frame(lst$covs)[, names(lst$params)[w]];
        if (length(w) == 1){
            d <- data.frame(d = d);
            names(d) <- names(lst$params)[w];
        }
        if (!is.dplyr){
            print(head(d), n = n);
        } else {
            print(dplyr::as.tbl(d), n = n, width = width);
        }
    }  else {
        p <- lst$params;
        if (length(env$pcov) > 0){
            p2 <- p[-env$pcov];
            print(p2)
            message("\nTime Varying Covariates");
            message(paste(names(p[env$pcov]), collapse=" "));
        } else {
            print(p);
        }
    }
    cat("\n\nInitial Conditions:\n")
    inits <- lst$inits[regexpr(regSens, names(lst$inits)) == -1];
    print(inits);
    cat("\n\nFirst part of data:\n")
    if (!is.dplyr){
        print(head(as.matrix(x), n = n));
    } else {
        print(dplyr::as.tbl(x), n = n, width = width);
    }
}

##' @author Matthew L.Fidler
##' @export
summary.solveRxODE <- function(object, ...){
    env <- attr(object, ".env");
    rxode <- env$env$out;
    lst <- c(list(params=env$params, inits=env$inits), env$extra.args)
    cat("Solved RxODE object\n");
    cat(sprintf("DLL: %s\n\n", rxDll(rxode)))
    cat("Model:\n");
    cat("################################################################################\n");
    cat(rxode$dll$modVars$model["model"]);
    cat("################################################################################\n");
    cat("Parameters:\n")
    w <- which((names(lst$params) %in% names(object)));
    if (length(w) > 0){
        print(lst$params[-w]);
        cat("\n\nSummary of time-varying covariates:\n");
        d <- as.data.frame(lst$covs)[, names(lst$params)[w]];
        if (length(w) == 1){
            d <- data.frame(d = d);
            names(d) <- names(lst$params)[w];
        }
        print(summary(d));
    }  else {
        print(lst$params);
    }
    cat("\n\nInitial conditions:\n")
    print(lst$inits);
    cat("\n\nSummary of solved data:\n")
    print(summary.data.frame(object))
}

##
## Accessing Solved Object
## accessComp <- function(obj, arg){
##     lst <- obj;
##     class(lst) <- "list";
##     if (any(names(obj) == arg)){
##         return(lst[[arg]]);
##     } else {
##         if (arg == "calcJac"){
##             return(length(rxModelVars(obj)$dfdy) > 0)
##         } else if (arg == "calcSens"){
##             return(length(rxModelVars(obj)$sens) > 0)
##         } else if (any(rxState(obj) == gsub(regIni, "", arg))){
##             arg <- gsub(regIni, "", arg);
##             ret <- rxInits(obj)[arg];
##             if (is.na(ret)){
##                 ret <- NA;
##                 names(ret) <- arg;
##                 return(ret)
##             } else {
##                 return(ret)
##             }
##         } else if (any(rxParams(obj) == arg)){
##             ret <- rxInits(obj)[arg];
##             if (is.na(ret)){
##                 ret <- NA;
##                 names(ret) <- arg;
##                 return(ret)
##             } else {
##                 return(ret)
##             }
##         } else {
##             return(NULL);
##         }
##     }
## }

##' @author Matthew L.Fidler
##' @export
`$.solveRxODE` <-  function(obj, arg, exact = TRUE){
    m <- as.data.frame(obj);
    ret <- m[[arg, exact = exact]];
    if (is.null(ret) & class(arg) == "character"){
        if (nchar(arg) > 6 && substr(arg, 1, 6) == "_sens_"){
            w <- which(gsub(regSens, "_sens_\\1_\\2", names(m)) == arg);
            if (length(w) == 1){
                return(m[, w]);
            } else {
                    return(NULL)
            }
        }
        ##
        ## Slows down
        ##
        ## w <- which(gsub(regSens, "\\1_\\2", names(m)) == arg);
        ## if (length(w) == 1){
        ##     return(m[, w]);
        ## }
        ##
        if (arg == "t"){
            return(m[["time"]]);
        } else {
            env <- attr(obj, ".env");
            tmp <- c(list(params=env$params, inits=env$inits), env$extra.args);
            if (regexpr(regToSens1, arg) != -1){
                ret <- m[[gsub(regToSens1, "d/dt(d(\\1)/d(\\2))", arg)]];
                if (!is.null(ret)){
                    return(ret)
                }
            }
            if (isTRUE(exact)){
                w <- which(names(tmp) == arg);
            } else {
                w <- which(regexpr(rex::rex(start, arg), names(tmp)) != -1)
                if (length(w) == 1 && !any(names(tmp) == arg) && is.na(exact)){
                    warning(sprintf("partial match of '%s' to '%s'", arg, names(tmp)[w]));
                }
            }
            if (length(w) == 1){
                return(tmp[[w]]);
            }
            if (any(names(tmp$param) == arg)){
                return(tmp$param[arg]);
            }
            if (any(names(tmp$init) == gsub(regIni, "", arg))){
                arg <- gsub(regIni, "", arg);
                return(tmp$init[arg]);
            }
            if (any(arg == names(tmp$events))){
                if (substr(arg, 0, 4) == "get."){
                    return(tmp$events[[arg]]);
                } else {
                    call <- as.list(match.call(expand.dots = TRUE));
                    env <- parent.frame();
                    return(function(..., .obj = obj, .objName = toString(call$obj), .objArg = toString(call$arg), .envir = env){
                        return(solveRxODE_updateEventTable(.obj, .objName, .objArg, ..., envir = .envir));
                    });
                }
            }
            return(NULL);
        }
    } else {
        return(rxTbl(ret));
    }
}

solveRxODE_updateEventTable <- function(obj, objName, name, ..., envir = parent.frame()){
    rxCat("Update with new event specification.\n")
    env <- attr(obj, ".env")
    events <- eventTable()
    tmp <- env$extra.args$events$copy()
    events$import.EventTable(tmp$get.EventTable())
    events[[name]](...);
    tmp <- update(obj, events = eval(events));
    assign(objName, tmp, envir = envir);
    ## envir$...RxODE...temp... <- tmp;
    ## eval(parse(text = sprintf("%s <- ...RxODE...temp...;", objName)), envir = envir);
    ## lst <- attr(tmp, "solveRxDll");
    ## envir$...RxODE...temp... <- lst;
    ## eval(parse(text = sprintf("%s <- ...RxODE...temp...;", objName)), envir = envir);
    ## rm("...RxODE...temp...", envir = envir);
    invisible()
}


##' @author Matthew L.Fidler
##' @export
`[.solveRxODE` <- function(x, i, j, drop){
    df <- as.data.frame(x);
    assign.names <- NULL
    if (!missing(j) && class(j) == "character"){
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
        if (!is.null(assign.names) && class(df) == "data.frame")
            names(df) <- assign.names
    } else if (!missing(i) && !missing(j) && missing(drop)){
        df <- df[i, j];
        if (!is.null(assign.names) && class(df) == "data.frame")
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
        if (!is.null(assign.names) && class(df) == "data.frame")
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
    if (class(new) == "EventTable"){
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
        if (class(value) == "EventTable"){
            cat("Update event table and solved object.\n");
            return(update(obj, events = eventTable))
        } else if (class(value) == "data.frame"){
        } else if (class(value) == "numeric"){
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
##'     interpolated.  Currently this can be \code{"linear"}
##'     interpolation (the default), which interpolates the covariate
##'     by solving the line between the observed covariates and
##'     extrapolating the new covariate value. The other possibility is
##'     \code{"constant"}, or Last observation carried forward.  In this
##'     approach, the last observation of the covariate is considered
##'     the current value of the covariate.
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
rxSolve <- function(object,                      # RxODE object
                    params=NULL,                      # Parameter
                    events=NULL,                      # Events
                    inits              = NULL,   # Initial Events
                    scale              = c(), #scale
                    covs               = NULL,   # Covariates
                    stiff              = TRUE,   # Is the system stiff
                    transit_abs        = NULL,  # Transit compartment absorption?
                    atol               = 1.0e-6, # Absoltue Tolerance for LSODA solver
                    rtol               = 1.0e-6, # Relative Tolerance for LSODA solver
                    maxsteps           = 5000,   # Maximum number of steps
                    hmin               = 0,      # Hmin
                    hmax               = NULL,   # Hmax
                    hini               = 0,      # Hini
                    maxordn            = 12,     # maxordn
                    maxords            = 5,      # maxords
                    ...,
                    covs_interpolation = c("linear", "constant"),
                    theta=numeric(), eta=numeric(), add.cov=FALSE) {
    ## rxSolve returns
    UseMethod("rxSolve");
} # end function rxSolve


##' @rdname rxSolve
##' @export
rxSolve.RxODE <- function(object, params=NULL, events=NULL, inits = NULL, scale=c(), covs = NULL, stiff = TRUE, transit_abs = NULL,
                          atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
                          maxords = 5, ..., covs_interpolation = c("linear", "constant"),
                          theta=numeric(), eta=numeric(), matrix=FALSE, add.cov=FALSE){
    return(object$solve(params, events, inits, scale, covs, stiff, transit_abs, atol, rtol, maxsteps, hmin, hmax, hini, maxordn, maxords,...,
                        covs_interpolation = covs_interpolation, theta=theta, eta=eta, matrix=matrix, add.cov=add.cov))
}

##' @rdname rxSolve
##' @export
rxSolve.solveRxODE <- function(object, params=NULL, events=NULL, inits = NULL, scale=c(), covs = NULL, stiff = TRUE, transit_abs = NULL,
                          atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
                          maxords = 5, ..., covs_interpolation = c("linear", "constant"),
                          theta=numeric(), eta=numeric(), matrix=FALSE, add.cov=FALSE){
    env <- attr(object, ".env");
    rxode <- env$env$out;
    return(rxode$solve(params, events, inits, scale, covs, stiff, transit_abs, atol, rtol, maxsteps, hmin, hmax, hini, maxordn, maxords,...,
                       covs_interpolation = covs_interpolation, theta=theta, eta=eta, matrix=matrix, add.cov=add.cov))
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
