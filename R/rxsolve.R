##' Solves a ODE equation
##'
##' This uses RxODE family of objects, file, or model specification to
##' solve a ODE system.
##'
##' @param object is a either a RxODE family of objects, or a file-name
##'     with a RxODE model specification, or a string with a RxODE
##'     model specification.
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
##' @param covs_interpolation specifies the interpolation method for
##'     time-varying covariates. When solving ODEs it often samples
##'     times outside the sampling time specified in \code{events}.
##'     When this happens, the time varying covariates are
##'     interpolated.  Currently this can be \code{"Linear"}
##'     interpolation (the default), which interpolates the covariate
##'     by solving the line between the observed covariates and
##'     extrapolating the new covariate value. The other possibility is
##'     \code{"LOCF"}, or Last observation carried forward.  In this
##'     approach, the last observation of the covariate is considered
##'     the current value of the covariate.
##'
##' @param stiff a logical (\code{TRUE} by default) indicating whether
##'     the ODE system is stifff or not.
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
##' @param atol a numeric absolute tolerance (1e-08 by default);
##'
##' @param rtol a numeric relative tolerance (1e-06 by default).
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
                    params,                      # Parameter
                    events,                      # Events
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
                    covs_interpolation = c("Linear", "LOCF")
                    ) {
    ## rxSolve returns
    UseMethod("rxSolve");
} # end function rxSolve

##' Update the solved object with any of the new parameters.
##'
##' This is a wrapper to the rxSolve method.
##'
##' @param object Object to be updated
##' @param ... Arguments to be updated, and resolved.
##'
##' @author Matthew L.Fidler
##' @export
update.solveRxDll <- function(object, ...){
    rxSolve(object, ...);
}

## FIXME: getCall perhaps

##' @rdname rxSolve
##' @export
rxSolve.solveRxDll <- function(object, params, events, inits, scale , covs, stiff, transit_abs, atol, rtol, maxsteps, hmin,
                               hmax, hini, maxordn, maxords, ..., covs_interpolation= c("Linear", "LOCF")){
    call <- as.list(match.call(expand.dots = TRUE));
    lst <- attr(object, "solveRxDll");
    lst <- lst[names(lst) != "matrix"];
    lst <- lst[names(lst) != "object"];
    for (n in names(call)[c(-1, -2)]){
        if (n == "params"){
            for (n2 in names(eval(call$params))){
                lst$params[n2] <-  eval(call$params)[n2];
            }
        } else if (n == "inits"){
            for (n2 in names(eval(call$inits))){
                lst$inits[n2] <- eval(call$inits)[n2]
            }
        } else {
            lst[[n]] <- call[[n]];
        }
    }
    lst$object <- object$object;
    return(do.call(rxSolve.rxDll, lst, envir = parent.frame(1)));
}

##' @rdname rxSolve
##' @export
rxSolve.RxODE <- function(object, params, events, inits = NULL, scale=c(), covs = NULL, stiff = TRUE, transit_abs = NULL,
                          atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
                          maxords = 5, ..., covs_interpolation = c("Linear", "LOCF")){
    rxSolve.rxDll(object$cmpMgr$rxDll(), params, events, inits, scale, covs, stiff, transit_abs, atol, rtol, maxsteps, hmin,
                  hmax, hini, maxordn, maxords, ..., covs_interpolation = covs_interpolation);
}
##' @rdname rxSolve
##' @export
rxSolve.RxCompilationManager <- function(object, params, events, inits = NULL, scale=c(), covs = NULL, stiff = TRUE,
                                         transit_abs = NULL, atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0,
                                         hmax = NULL, hini = 0, maxordn = 12, maxords = 5, ...,
                                         covs_interpolation = c("Linear", "LOCF")){
    rxSolve.rxDll(object$rxDll(), params, events, inits, scale, covs, stiff, transit_abs, atol, rtol, maxsteps, hmin, hmax,
                  hini, maxordn, maxords, ...,
                  covs_interpolation = covs_interpolation);
}
##' @rdname rxSolve
##' @export
rxSolve.character <- function(object, params, events, inits = NULL, scale=c(), covs = NULL, stiff = TRUE, transit_abs = NULL,
                              atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
                              maxords = 5, ..., covs_interpolation = c("Linear", "LOCF")){
    rxSolve.rxDll(rxCompile(object), params, events, inits, scale, covs, stiff, transit_abs, atol, rtol, maxsteps, hmin, hmax,
                  hini, maxordn, maxords, ..., covs_interpolation = covs_interpolation);
}
##' @rdname rxSolve
##' @export
rxSolve.rxDll <- function(object, params, events, inits = NULL, scale = c(),
                          covs = NULL, stiff = TRUE, transit_abs = NULL,
                          atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
                          maxords = 5, ..., covs_interpolation = c("Linear", "LOCF")){
    ## rxSolve.rxDll returns a solved object
    if (missing(events) && class(params) == "EventTable"){
        events <- params;
        params <- c();
    }
    if (is.null(transit_abs)){
        transit_abs <- rxModelVars(object)$podo;
        if (transit_abs){
            warning("Assumed transit compartment model since 'podo' is in the model.")
        }
    }
    ## Params and inits passed
    extra.args <- list(events = events$copy(),
                       covs = covs, stiff = stiff,
                       transit_abs = transit_abs, atol = atol, rtol = rtol, maxsteps = maxsteps,
                       hmin = hmin, hmax = hmax, hini = hini, maxordn = maxordn, maxords = maxords,
                       covs_interpolation = covs_interpolation, ...);
    extra.args <- extra.args[names(extra.args) != "counts"];
    event.table <- events$get.EventTable()
    if (!is.numeric(maxordn))
        stop("`maxordn' must be numeric")
    if (maxordn < 1 || maxordn > 12)
        stop("`maxordn' must be >1 and < = 12")
    if (!is.numeric(maxords))
        stop("`maxords' must be numeric")
    if (maxords < 1 || maxords > 5)
        stop("`maxords' must be >1 and < = 5")
    if (!is.numeric(rtol))
        stop("`rtol' must be numeric")
    if (!is.numeric(atol))
        stop("`atol' must be numeric")
    if (!is.numeric(hmin))
        stop("`hmin' must be numeric")
    if (hmin < 0)
        stop("`hmin' must be a non-negative value")
    if (is.null(hmax)){
        if (is.null(event.table$time)){
            hmax <- 0;
        } else {
            hmax <- max(abs(diff(event.table$time)))
        }
    }
    if (!is.numeric(hmax))
        stop("`hmax' must be numeric")
    if (hmax < 0)
        stop("`hmax' must be a non-negative value")
    if (hmax == Inf)
        hmax <- 0
    if (!is.null(hini)){
        if (hini < 0)
            stop("`hini' must be a non-negative value")
    } else {
        hini <- 0;
    }
    ## preserve input arguments.
    inits <- rxInits(object, inits, rxState(object), 0);
    params <- rxInits(object, params, rxParams(object), NA, !is.null(covs));

    if (!is.null(covs)){
        cov <- as.matrix(covs);
        pcov <- sapply(dimnames(cov)[[2]], function(x){
            w <- which(x == names(params));
            if (length(w) == 1){
                return(w)
            } else {
                return(0);
            }
        })
        n_cov <- length(pcov);
        ## Now check if there is any unspecified parameters by either covariate or parameter
        w <- which(is.na(params));
        if (!all(names(params)[w] %in% dimnames(cov)[[2]])){
            stop("Some model specified variables were not specified by either a covariate or parameter");
        }
        ## Assign all parameters matching a covariate to zero.
        for (i in pcov){
            if (i > 0){
                params[i] <- 0;
            }
        }
        covnames <- dimnames(cov)[[2]]
    } else {
        ## For now zero out the covariates
        pcov <- c();
        cov <- c();
        n_cov <- 0;
        covnames <- c();
    }
    lhs_vars <- rxLhs(object);
    if (is.null(inits)){
        n <- rxState(object)
        inits <- rep(0.0, length(n));
        names(inits) <- n;
    }
    s <- as.list(match.call(expand.dots = TRUE))
    wh <- grep(pattern = "[Ss]\\d+$", names(s))
    if (length(scale) > 0 && length(wh) > 0){
        stop("Cannot specify both 'scale=c(...)' and S#=, please pick one for to scale the ODE compartments.")
    }
    ## HACK: fishing scaling variables "S1 S2 S3 ..." from params call
    ## to solve(). Maybe define a "scale = c(central = 7.6, ...)" argument
    ## similar to "params = "?
    scaler.ix <- c()
    if (length(wh) > 0) {
        scaler.ix <- as.numeric(substring(names(s)[wh], 2))
        if (any(duplicated(scaler.ix))){
            stop("Duplicate scaling factors found.");
        }
        scale <- unlist(s[wh]);
        if (any(length(inits) < scaler.ix)){
            warning(sprintf("scaler variable(s) above the number of compartments: %s.",
                            paste(paste0("S", scaler.ix[scaler.ix > length(inits)]), collapse=", ")))
            scale <- scale[scaler.ix < length(inits)]
            scaler.ix <- scaler.ix[scaler.ix < length(inits)];
        }
        names(scale) <- rxState(object)[scaler.ix];
    }
    scale <- c(scale);
    if (length(covs_interpolation) > 1){
        isLocf <- 0;
    } else if (covs_interpolation == "LOCF"){
        isLocf <- 1;
    } else if (covs_interpolation != "Linear"){
        stop("Unknown covariate interpolation specified.");
    }
    ## may need to reload (e.g., when we re-start R and
    ## re-instantiate the RxODE object from a save.image.
    ## cmpMgr$dynLoad()
    rxLoad(object);
    if (event.table$time[1] != 0){
        warning(sprintf("The initial conditions are at t = %s instead of t = 0.", event.table$time[1]))
    }
    ## Ensure that inits and params have names.
    names(inits) <- rxState(object);
    names(params) <- rxParams(object);
    ret <- object$.call(rxTrans(object)["ode_solver_sexp"],
                        ## Parameters
                        params,
                        inits,
                        lhs_vars,
                        ## events
                        event.table$time,
                        as.integer(event.table$evid),
                        as.double(event.table$amt[event.table$evid>0]),
                        ## Covariates
                        as.integer(pcov),
                        as.double(cov),
                        as.integer(isLocf),
                        ## Solver options (double)
                        as.double(atol),
                        as.double(rtol),
                        as.double(hmin),
                        as.double(hmax),
                        as.double(hini),
                        ## Solver options ()
                        as.integer(maxordn),
                        as.integer(maxords),
                        as.integer(maxsteps),
                        as.integer(stiff),
                        as.integer(transit_abs),
                        ## Passed to build solver object.
                        object,
                        extra.args);
    rc <- attr(ret, "solveRxDll")$counts["rc"];

    ## Subset to observations only.
    attr(ret, "solveRxDll")$matrix <- attr(ret, "solveRxDll")$matrix[events$get.obs.rec(), ];

    if (rc != 0)
        stop(sprintf("could not solve ODE, IDID = %d (see further messages)", rc))

    for (d in names(scale)){
        attr(ret, "solveRxDll")$matrix[, d] <- attr(ret, "solveRxDll")$matrix[, d] / scale[d];
    }

    return(ret);
} # end function rxSolve.rxDll

##' @author Matthew L.Fidler
##' @export
print.solveRxDll <- function(x, ...){
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
    lst <- attr(x, "solveRxDll")
    cat("Solved RxODE object\n");
    cat(sprintf("Dll: %s\n\n", rxDll(lst$object)))
    cat("Parameters:\n")
    is.dplyr <- requireNamespace("dplyr", quietly = TRUE) && getOption("RxODE.display.tbl", TRUE);
    w <- which((names(lst$params) %in% names(as.data.frame(x))))
    if (length(w) > 0){
        print(lst$params[-w]);
        cat("\n\nFirst Part of Time Varying Covariates:\n");
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
        print(lst$params);
    }
    cat("\n\nInitial Conditions:\n")
    print(lst$inits);
    cat("\n\nFirst part of data:\n")
    if (!is.dplyr){
        print(head(as.matrix(x), n = n));
    } else {
        print(dplyr::as.tbl(x), n = n, width = width);
    }
}

##' @author Matthew L.Fidler
##' @export
summary.solveRxDll <- function(object, ...){
    lst <- attr(object, "solveRxDll")
    cat("Solved RxODE object\n");
    cat(sprintf("Dll: %s\n\n", rxDll(lst$object)))
    cat("Model:\n");
    cat("################################################################################\n");
    cat(rxModelVars(object)$model["model"]);
    cat("################################################################################\n");
    cat("Parameters:\n")
    w <- which((names(lst$params) %in% names(as.data.frame(object))))
    if (length(w) > 0){
        print(lst$params[-w]);
        cat("\n\nSummary of Time Varying Covariates:\n");
        d <- as.data.frame(lst$covs)[, names(lst$params)[w]];
        if (length(w) == 1){
            d <- data.frame(d = d);
            names(d) <- names(lst$params)[w];
        }
        print(summary(d));
    }  else {
        print(lst$params);
    }
    cat("\n\nInitial Conditions:\n")
    print(lst$inits);
    cat("\n\nSummary of solved data:\n")
    print(summary(as.data.frame(object)))
}


##' @author Matthew L.Fidler
##' @export
as.matrix.solveRxDll <- function(x, ...){
    lst <- attr(x, "solveRxDll")
    return(lst$matrix);
}


##' @author Matthew L.Fidler
##' @export
as.data.frame.solveRxDll <- function(x, row.names = NULL, optional = FALSE, ...,
                                     stringsAsFactors = default.stringsAsFactors()){
    return(as.data.frame(as.matrix(x), row.names = row.names, optional = optional, ...,
                         stringsAsFactors = stringsAsFactors));
}

##' Use the as_data_frame method for solved object x
##'
##' @param x Solved RxODE object
##' @return data frame of solved object.
##' @author Matthew L. Fidler
as_data_frame.solveRxDll <- function(x){
    return(tibble::as_data_frame(as.matrix(x)));
}


##' @author Matthew L.Fidler
##' @importFrom utils head
##' @export
head.solveRxDll <- function(x, n = 6L, ...){
    return(utils::head.matrix(as.matrix(x), n = n, ...));
}


##' @author Matthew L.Fidler
##' @importFrom utils tail
##' @export
tail.solveRxDll <- function(x, n = 6L, addrownums = TRUE, ...){
    return(utils::tail.matrix(as.matrix(x), n = n, addrownums = addrownums, ...));
}

##' Convert Solved RxODE object to tbl
##' @param x Solved RxDll object
##' @param ... other arguments (ignored)
##' @author Matthew L.Fidler
##' @export as.tbl.solveRxDll
as.tbl.solveRxDll <- function(x, ...){
    return(dplyr::as.tbl(as.data.frame(x)));
}

solveRxDll_updateEventTable <- function(obj, objName, name, ..., envir = parent.frame()){
    if (getOption("RxODE.verbose", TRUE)){ ## nocov start
        cat("Update with new event specification.\n");
    } ## nocov end
    tmp <- attr(obj, "solveRxDll");
    events <- tmp$events;
    events[[name]](...);
    tmp <- rxSolve.solveRxDll(obj, events = eval(events));
    envir$...RxODE...temp... <- tmp;
    eval(parse(text = sprintf("%s <- ...RxODE...temp...;", objName)), envir = envir);
    lst <- attr(tmp, "solveRxDll");
    envir$...RxODE...temp... <- lst;
    eval(parse(text = sprintf("attr(%s, \"solveRxDll\") <- ...RxODE...temp...;", objName)), envir = envir);
    rm("...RxODE...temp...", envir = envir);
    invisible()
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
"$<-.solveRxDll" <- function(obj, arg, value){
    if (arg == "t"){
        arg <- "time";
    }
    lst <- attr(obj, "solveRxDll")
    iarg <- gsub(regIni, "", arg);
    if (arg == "time"){
        if (class(value) == "EventTable"){
            cat("Update event table and solved object.\n");
            return(rxSolve(obj, events = eventTable))
        } else if (class(value) == "data.frame"){
        } else if (class(value) == "numeric"){
            if (getOption("RxODE.verbose", TRUE)){ ## nocov start
                cat("Updating sampling times in the event table updating object.\n");
            } ## nocov end
            eventTable <- lst$events$copy();
            eventTable$clear.sampling();
            eventTable$add.sampling(value);
            return(rxSolve(obj, events = eventTable));
        }
    } else if (arg != iarg && any(rxState(obj$object) == iarg) && length(value) == 1){
        if (getOption("RxODE.verbose", TRUE)){ ## nocov start
            cat("Updating object with new initial conditions.\n");
        } ## nocov end
        inits <- c(value);
        names(inits) <- gsub(regIni, "", arg);
        return(rxSolve(obj, inits = inits));
    } else if (any(rxParams(obj$object) == arg)){
        if (getOption("RxODE.verbose", TRUE)){ ## nocov start
            cat("Updating object with new paramter values.\n");
        } ## nocov end
        if (length(value) == 1){
            covs <- as.data.frame(lst$covs);
            if (any(names(covs) == arg)){
                cat(sprintf("Changing time varying covariate \"%s\" to a simple parameter value %s\n", arg, value));
                ncovs <- names(covs);
                covs <- as.data.frame(covs[, names(covs) != arg]);
                names(covs) <- ncovs[ncovs != arg];
                params <- c(value);
                names(params) <- arg;
                return(rxSolve.solveRxDll(obj, params = params, covs = covs));
            } else {
                params <- c(value);
                names(params) <- arg;
                return(rxSolve.solveRxDll(obj, params = params));
            }
        } else if (length(value) == length(lst$events$get.sampling()$time)){
            cat(sprintf("Changing simple parameter \"%s\" to a time-varying covariate.\n", arg));
            covs <- as.data.frame(lst$covs);
            covs[[arg]] <- value;
            return(rxSolve.solveRxDll(obj, covs = covs));
        }
    } else if (arg == "params"){
        if (getOption("RxODE.verbose", TRUE)){ ## nocov start
            cat("Updating object with new paramter values.\n");
        } ## nocov end
        return(rxSolve.solveRxDll(obj, params = value));
    } else if (arg == "inits"){
        if (getOption("RxODE.verbose", TRUE)){ ## nocov start
            cat("Updating object with new initial conditions.\n");
        } ## nocov end
        return(rxSolve.solveRxDll(obj, inits = value));
    } else if (any(arg == names(lst))){
        args <- list(obj, value);
        names(args) <- c("object", arg);
        if (getOption("RxODE.verbose", TRUE)){ ## nocov start
            cat(sprintf("Updating object with new solving argument %s = %s.\n", arg, value))
        } ## nocov end
        return(do.call("rxSolve", args, envir = parent.frame(1)))
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
"[[<-.solveRxDll" <- function(obj, arg, value){
    return("$<-.solveRxDll"(obj, arg, value = value))
}

##' Assign solved objects using the [] syntax
##' @param obj solved object
##' @param arg1 first argument
##' @param arg2 second argument
##' @param value value assumed
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"[<-.solveRxDll" <- function(obj, arg1, arg2, value){
    if (any(arg2 == c(1, "t", "time")) && missing(arg1)){
        obj$time <- value
        return(obj);
    }
    lst <- attr(obj, "solveRxDll")
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

##' Assign rownames to rxSolve object
##'
##' row.names(x) <- value;
##'
##' @param x rxode object
##' @param value value assigned.
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"row.names<-.solveRxDll" <- function(x, value){
    lst <- attr(x, "solveRxDll")
    df <- as.data.frame(x);
    m <- as.matrix("row.names<-.data.frame"(df, value));
    lst$matrix <- m;
    attr(x, "solveRxDll") <- lst;
    return(x);
}

##' Assign dimensions of rxSolve object
##'
##' dim(x) <- value
##'
##' @param x rxSolve object
##' @param value dimensions
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"dim<-.solveRxDll" <- function(x, value){
    lst <- attr(x, "solveRxDll")
    m <- as.matrix(x);
    dim(m) <- value;
    lst$matrix <- m;
    attr(x, "solveRxDll") <- lst;
    return(x);
}


##' Assign dimension names for rxSolve object
##'
##' dimnames(x) <- value
##'
##' @param x rxSolve object
##' @param value dimension names assigned
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
"dimnames<-.solveRxDll" <- function(x, value){
    lst <- attr(x, "solveRxDll")
    m <- as.matrix(x);
    dimnames(m) <- value;
    lst$matrix <- m;
    attr(x, "solveRxDll") <- lst;
    return(x);
}

##' @author Matthew L.Fidler
##' @export
"split<-.solveRxDll" <- function(x, f, drop = FALSE, ..., value){
    "split<-"(as.data.frame(x), f, drop, ..., value = value);
}

## FIXME rbind, cbind could be possible...

##' @rdname cbind.solveRxDll
##' @author Matthew L.Fidler
##' @export
rbind.solveRxDll <- function(...){
    stop("rbind is unsupported.  First convert to a data.frame with as.data.frame(x).")
}
##' cbind/rbind solveRxDll
##'
##' Cbind/rbind is disabled for RxOde solved objects.  Use as.data.frame(x)
##' to use cbind/rbind.
##'
##' @param ... ignored parameters
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
cbind.solveRxDll <- function(...){
    stop("cbind is unsupported.  First convert to a data.frame with as.data.frame(x).")
}
