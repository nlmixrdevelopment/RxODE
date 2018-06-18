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
##' @param method The method for solving ODEs.  Currently this supports:
##'
##' \itemize{
##' \item \code{"liblsoda"} thread safe lsoda.  This supports parallel
##'            thread-based solving, and ignores user Jacobian specification.
##' \item \code{"lsoda"} -- LSODA solver.  Does not support parallel thread-based
##'       solving, but allows user Jacobian specification.
##' \item \code{"dop853"} -- DOP853 solver.  Does not support parallel thread-based
##'         solving nor user Jacobain specification
##' }
##'
##' @param transitAbs boolean indicating if this is a transit
##'     compartment absorption
##'
##' @param atol a numeric absolute tolerance (1e-8 by default) used
##'     by the ODE solver to determine if a good solution has been
##'     achieved;  This is also used in the solved linear model to check
##'     if prior doses do not add anything to the solution.
##'
##' @param rtol a numeric relative tolerance (1e-6 by default) used
##'     by the ODE solver to determine if a good solution has been
##'     achieved. This is also used in the solved linear model to check
##'      if prior doses do not add anything to the solution.
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
##' @param cores Number of cores used in parallel ODE solving.  This
##'     defaults to the number or system cores determined by
##'     \code{\link{rxCores}} for methods that support parallel
##'     solving (ie thread-safe methods like "liblsoda").
##'
##' @param covsInterpolation specifies the interpolation method for
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
##' @param addCov A boolean indicating if covariates should be added
##'     to the output matrix or data frame. By default this is
##'     disabled.
##'
##' @param matrix A boolean inticating if a matrix should be returned
##'     instead of the RxODE's solved object.
##'
##' @param sigma Named sigma covariance or Cholesky decomposition of a
##'     covariance matrix.  The names of the columns indicate
##'     parameters that are simulated.  These are simulated for every
##'     observation in the solved system.
##'
##' @param sigmaDf Degrees of freedom of the sigma t-distribution.  By
##'     default it is equivalent to \code{Inf}, or a normal distribution.
##'
##' @param nCoresRV Number of cores used for the simulation of the
##'     sigma variables.  By default this is 1. This uses the package
##'     \code{\link[mvnfast]{rmvn}} and \code{\link[mvnfast]{rmvt}}.
##'     To reproduce the results you need to run on the same platform
##'     with the same number of cores. This is the reason this is set
##'     to be one, regardless of what the number of cores are used in
##'     threaded ODE solving.
##'
##' @param sigmaIsChol Boolean indicating if the sigma is in the
##'     Cholesky decomposition instead of a symmetric covariance
##'
##' @param nDisplayProgress An integer indicating the minimum number
##'     of c-based solves before a progress bar is shown.  By default
##'     this is 10,000.
##'
##' @param amountUnits This supplies the dose units of a data frame
##'     supplied instead of an event table.  This is for importing the
##'     data as an RxODE event table.
##'
##' @param timeUnits This supplies the time units of a data frame
##'     supplied instead of an event table.  This is for importing the
##'     data as an RxODE event table.
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
##' @param theta A vector of parameters that will be named THETA[#] and
##'     added to parameters
##'
##' @param eta A vector of parameters that will be named ETA[#] and
##'     added to parameters
##'
##' @param addDosing Boolean indicating if the solve should add RxODE
##'     evid and amt columns.  This will also include dosing
##'     information and estimates at the doses.  Be default, RxODE
##'     only includes estimates at the observations. (default
##'     \code{FALSE}).
##'
##' @param updateObject This is an internally used flag to update the
##'     RxODE solved object (when supplying an RxODE solved object) as
##'     well as returning a new object.  You probably should not
##'     modify it's \code{FALSE} default unless you are willing to
##'     have unexpected results.
##'
##' @param doSolve Internal flag.  By default this is \code{TRUE},
##'     when \code{FALSE} a list of solving options is returned.
##'
##' @param returnType This tells what type of object is returned.  The currently supported types are:
##' \itemize{
##' \item \code{"rxSolve"} (default) will return a reactive data frame
##'      that can change easily change different pieces of the solve and
##'      update the data frame.  This is the currently standard solving
##'      method in RxODE,  is used for \code{rxSolve(object, ...)}, \code{solev(object,...)},
##' \item \code{"data.frame"} -- returns a plain, non-reactive data
##'      frame; Currently very slightly Faster than \code{returnType=\"matrix\"}
##' \item \code{"matrix"} -- returns a plain matrix with column names attached
##'     to the solved object.  This is what is used \code{object$run} as well as \code{object$solve}
##' }
##'
##' @inheritParams rxSimThetaOmega
##'
##' @inheritParams stats::simulate
##'
##' @param a when using \code{solve}, this is equivalent to the
##'     \code{object} argument.  If you specify \code{object} later in
##'     the argument list it overwrites this parameter.
##'
##' @param b when using \code{solve}, this is equivalent to the
##'     \code{params} argument.  If you specify \code{params} as a
##'     named argument, this overwrites the output
##'
##' @param nsim represents the number of simulations.  For RxODE, if
##'     you supply single subject event tables (created with
##'     eventTable)
##'
##' @param setupOnly Only setup the internal C structure, do not
##'     solve.  After setting it up, and using the structure in C, it
##'     needs to be freed by \code{\link{rxSolveFree}}.
##'
##' @return An \dQuote{rxSolve} solve object that stores the solved
##'     value in a matrix with as many rows as there are sampled time
##'     points and as many columns as system variables (as defined by
##'     the ODEs and additional assignments in the RxODE model code).
##'     It also stores information about the call to allow dynamic
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
##' @seealso \code{\link{RxODE}}
##' @author Matthew Fidler, Melissa Hallow and  Wenping Wang
##' @export
rxSolve <- function(object, ...){
    UseMethod("rxSolve")
}
##' @rdname rxSolve
##' @export
rxSolve.default <- function(object, params=NULL, events=NULL, inits = NULL, scale = NULL,
                    covs = NULL, method = c("liblsoda", "lsoda", "dop853"),
                    transitAbs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
                    maxsteps = 5000L, hmin = 0L, hmax = NULL, hini = 0, maxordn = 12L, maxords = 5L, ...,
                    cores,
                    covsInterpolation = c("linear", "locf", "nocb", "midpoint"),
                    addCov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                    nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress=10000L,
                    amountUnits = NA_character_, timeUnits = "hours", stiff,
                    theta = NULL, eta = NULL, addDosing=FALSE, updateObject=FALSE, doSolve=TRUE,
                    omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                    nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                    nStud = 1L, dfSub=0.0, dfObs=0.0, returnType=c("rxSolve", "matrix", "data.frame"),
                    seed=NULL, nsim=NULL, setupOnly=FALSE){
    .xtra <- list(...);
    if (is.null(transitAbs) && !is.null(.xtra$transit_abs)){
        transitAbs <- .xtra$transit_abs;
    }
    if (missing(updateObject) && !is.null(.xtra$update.object)){
        updateObject <- .xtra$update.object;
    }
    if (missing(doSolve) && !is.null(.xtra$do.solve)){
        doSolve <- .xtra$do.solve;
    }
    if (missing(covsInterpolation) && !is.null(.xtra$covs_interpolation)){
        covsInterpolation <- .xtra$covs_interpolation;
    }
    if (missing(addCov) && !is.null(.xtra$add.cov)){
        addCov <- .xtra$add.cov;
    }
    if (!is.null(seed)){
        set.seed(seed);
    }
    if (!is.null(nsim)){
        if (rxIs(params, "eventTable") || rxIs(events, "eventTable") && nSub == 1L){
            nSub <- nsim;
        } else if (nStud == 1L){
            nStud <- nsim;
        }

    }
    if (!doSolve){
        .modVars <- RxODE::rxModelVars(object);
        .state <- .modVars$state;
        .lhs <- .modVars$lhs;
        .pars <- .modVars$params;
        .stateIgnore <- .modVars$state.ignore
        if (!is.null(params)){
            if (is.null(events) && is(params, "EventTable")){
                events <- params;
                params <- c();
            }
        }
        if (is.null(transitAbs)){
            transitAbs <- .modVars$podo;
            if (transitAbs){
                warning("Assumed transit compartment model since 'podo' is in the model.")
            }
        }
        if (!is(params, "numeric")){
            .n <- names(params);
            params <- as.double(params);
            names(params) <- .n;
        }
        if (missing(stiff)) stiff <- TRUE;
        ## Params and inits passed
        extra.args <- list(events = events$copy(),
                           covs = covs, stiff = stiff,
                           transitAbs = transitAbs, atol = atol, rtol = rtol, maxsteps = maxsteps,
                           hmin = hmin, hmax = hmax, hini = hini, maxordn = maxordn, maxords = maxords,
                           covsInterpolation = covsInterpolation, addCov=addCov, ...);
        params <- c(params, RxODE::rxThetaEta(theta, eta));
        .eventTable <- events$get.EventTable()
        if (!is.numeric(maxordn))
            stop("'maxordn' must be numeric.")
        if (maxordn < 1 || maxordn > 12)
            stop("'maxordn' must be >1 and < = 12.")
        if (!is.numeric(maxords))
            stop("'maxords' must be numeric.")
        if (maxords < 1 || maxords > 5)
            stop("'maxords' must be >1 and < = 5.")
        if (!is.numeric(rtol))
            stop("'rtol' must be numeric.")
        if (!is.numeric(atol))
            stop("'atol' must be numeric.")
        if (!is.numeric(hmin))
            stop("'hmin' must be numeric.")
        if (hmin < 0)
            stop("'hmin' must be a non-negative value.")
        if (is.null(hmax)){
            if (is.null(.eventTable$time) || length(.eventTable$time) == 1){
                hmax <- 0;
            } else {
                hmax <- max(abs(diff(.eventTable$time)))
            }
        }
        if (!is.numeric(hmax))
            stop("'hmax' must be numeric.")
        if (hmax < 0)
            stop("'hmax' must be a non-negative value.")
        if (hmax == Inf)
            hmax <- 0
        if (!is.null(hini)){
            if (hini < 0)
                stop("'hini' must be a non-negative value.")
        } else {
            hini <- 0;
        }
        ## preserve input arguments.
        inits <- RxODE::rxInits(object, inits, .state, 0);
        params <- RxODE::rxInits(object, params, .pars, NA, !is.null(covs));
        if (!is.null(covs)){
            .cov <- as.matrix(covs);
            .covLen <- dim(.cov)[1];
            if (.covLen !=  length(.eventTable$time)){
                .samplingTime <- events$get.sampling()$time;
                if (.covLen != length(.samplingTime))
                    stop("Covariate length need to match the sampling times or all the times in the event table.");
                .lst <- as.matrix(do.call("cbind", lapply(seq(1L, dim(.cov)[2]), function(i){
                                                      f <- stats::approxfun(.samplingTime, .cov[, i])
                                                      return(f(.eventTable$time))
                                                  })))
                dimnames(.lst) <- list(NULL, dimnames(.cov)[[2]]);
                .cov <- .lst;
            }
            .pcov <- sapply(dimnames(.cov)[[2]], function(x){
                .w <- which(x == names(params));
                if (length(.w) == 1){
                    return(.w)
                } else {
                    return(0);
                }
            })
            ## Now check if there is any unspecified parameters by either covariate or parameter
            .w <- which(is.na(params));
            if (!all(names(params)[.w] %in% dimnames(.cov)[[2]])){
                print(params)
                stop("Some model specified variables were not specified by either a covariate or parameter");
            }
            ## Assign all parameters matching a covariate to zero.
            for (i in .pcov){
                if (i > 0){
                    params[i] <- 0;
                }
            }
        } else {
            ## For now zero out the covariates
            .pcov <- c();
            .cov <- c();
        }
        if (is.null(inits)){
            .n <- .state;
            inits <- rep(0.0, length(.n));
            names(inits) <- .n;
        }
        .s <- as.list(match.call(expand.dots = TRUE))
        .wh <- grep(pattern = "[Ss]\\d+$", names(.s))
        if (length(scale) > 0 && length(.wh) > 0){
            stop("Cannot specify both 'scale=c(...)' and S#=, please pick one to scale the ODE compartments.")
        }
        ## HACK: fishing scaling variables "S1 S2 S3 ..." from params call
        ## to solve(). Maybe define a "scale = c(central = 7.6, ...)" argument
        ## similar to "params = "?
        .scalerIx <- c()
        if (length(.wh) > 0) {
            .scalerIx <- as.numeric(substring(names(.s)[.wh], 2))
            if (any(duplicated(.scalerIx))){
                stop("Duplicate scaling factors found.");
            }
            scale <- unlist(.s[.wh]);
            if (any(length(inits) < .scalerIx)){
                warning(sprintf("Scaler variable(s) above the number of compartments: %s.",
                                paste(paste0("S", .scalerIx[.scalerIx > length(inits)]), collapse=", ")))
                scale <- scale[.scalerIx < length(inits)]
                .scalerIx <- .scalerIx[.scalerIx < length(inits)];
            }
            names(scale) <- .state[.scalerIx];
        }
        scale <- c(scale);
        scale <- RxODE::rxInits(object, scale, .state, 1, noini=TRUE);
        .isLocf <- 0L;
        if (length(covsInterpolation) > 1){
            .isLocf <- 0L;
        } else if (covsInterpolation == "constant"){
            .isLocf <- 1L;
        } else if (covsInterpolation == "NOCB"){
            .isLocf <- 2L;
        } else if (covsInterpolation == "midpoint"){
            .isLocf <- 3L;
        } else if (covsInterpolation != "linear"){
            stop("Unknown covariate interpolation specified.");
        }
        ## if (.eventTable$time[1] != 0){
        ##     warning(sprintf("The initial conditions are at t = %s instead of t = 0.", .eventTable$time[1]))
        ## }
        ## Ensure that inits and params have names.
        names(inits) <- .state
        names(params) <- .pars;

        .time <- as.double(.eventTable$time);
        .evid <- as.integer(.eventTable$evid);
        .amt <- as.double(.eventTable$amt[.eventTable$evid>0]);
        ## Covariates
        .pcov <- as.integer(.pcov);
        .cov <- as.double(.cov);
        .isLocf <- as.integer(.isLocf);
        ## Solver options (double)
        atol <- as.double(atol);
        rtol <- as.double(rtol);
        hmin <- as.double(hmin);
        hmax <- as.double(hmax);
        hini <- as.double(hini);
        ## Solver options ()
        maxordn <- as.integer(maxordn);
        maxords <- as.integer(maxords);
        maxsteps <- as.integer(maxsteps);
        stiff <- as.integer(stiff);
        transitAbs <- as.integer(transitAbs);
        .doMatrix <- as.integer(matrix);
        addCov <- as.integer(addCov)
        ret <- list(params=params,
                    inits=inits,
                    lhs_vars=.lhs,
                    ## events
                    time=.time,
                    evid=.evid,
                    amt=.amt,
                    ## Covariates
                    pcov=.pcov,
                    covs=.cov,
                    isLocf=.isLocf,
                    ## Solver options (double)
                    atol=atol,
                    rtol=rtol,
                    hmin=hmin,
                    hmax=hmax,
                    hini=hini,
                    ## Solver options ()
                    maxordn=maxordn,
                    maxords=maxords,
                    maxsteps=maxsteps,
                    stiff=stiff,
                    transitAbs=transitAbs,
                    ## Passed to build solver object.
                    object=object,
                    extra.args=extra.args,
                    scale=scale,
                    events=events,
                    event.table=.eventTable,
                    do.matrix=.doMatrix,
                    addCov=addCov,
                    state.ignore=.stateIgnore);
        return(ret);
    }
    ## stiff = TRUE, transitAbs = NULL,
    ## atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
    ## maxords = 5, ..., covsInterpolation = c("linear", "constant", "NOCB", "midpoint"),
    ## theta=numeric(), eta=numeric(), matrix=TRUE,addCov=FALSE,
    ## inC=FALSE, counts=NULL, doSolve=TRUE
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
    } else {
        method <- match.arg(method);
    }
    if (!missing(returnType)){
        .matrixIdx <- c("rxSolve"=0, "matrix"=1, "data.frame"=2);
        matrix <- .matrixIdx[match.arg(returnType)];
    } else if (!is.null(.xtra$return.type)){
        .matrixIdx <- c("rxSolve"=0, "matrix"=1, "data.frame"=2);
        matrix <- .matrixIdx[.xtra$return.type];
    } else {
        matrix <- as.integer(matrix);
    }
    .methodIdx <- c("lsoda"=1, "dop853"=0, "liblsoda"=2);
    method <- as.integer(.methodIdx[method]);
    if (Sys.info()[["sysname"]] == "SunOS" && method == 2){
        method <- 1;
    }
    if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1];
    covsInterpolation <- tolower(match.arg(covsInterpolation,
                                           c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")))
    if (covsInterpolation == "constant") covsInterpolation <- "locf";
    covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
    if (any(duplicated(names(.xtra)))){
        stop("Duplicate arguments do not make sense.");
    }
    if (missing(cores)){
        cores <- RxODE::rxCores();
    }
    .nms <- names(as.list(match.call())[-1]);
    .Call(`_RxODE_rxSolveCsmall`, object, .nms, .xtra,
          params, events, inits, scale, covs,
          list(method, #0
               transitAbs, #1
               atol, #2
               rtol, #3
               maxsteps, #4
               hmin, #5
               hmax, #6
               hini, #7
               maxordn, #8
               maxords, #9
               cores, #10
               covsInterpolation, #11
               addCov, #12
               matrix, #13
               sigma, #14
               sigmaDf, #15
               nCoresRV, #16
               sigmaIsChol, nDisplayProgress, amountUnits,
               timeUnits, addDosing, theta, eta, updateObject,
               doSolve, omega, omegaDf, omegaIsChol, nSub, thetaMat,
               thetaDf, thetaIsChol, nStud, dfSub, dfObs,
               as.integer(setupOnly)));
}

##' @rdname rxSolve
##' @export
update.rxSolve <- function(object, ...){
    rxSolve(object, ...);
}

##' @rdname rxSolve
##' @export
predict.RxODE <- function(object, ...){
    rxSolve(object, ...);
}

##' @rdname rxSolve
##' @export
predict.rxSolve <- predict.RxODE

##' @importFrom stats simulate

##' @rdname rxSolve
##' @export
simulate.RxODE <- function(object, nsim = 1L, seed = NULL, ...){
    rxSolve(object, ..., seed=seed, nsim=nsim);
}
##' @rdname rxSolve
##' @export
simulate.rxSolve <- simulate.RxODE

##' @rdname rxSolve
##' @export
solve.rxSolve <- function(a, b, ...){
    lst <- as.list(match.call()[-1])
    n <- names(lst)
    if (!missing(a)){
        n[n == "a"] <- "";
    }
    if (!missing(b)){
        n[n == "b"] <- "";
    }
    names(lst) <- n
    do.call("rxSolve", lst, envir=parent.frame(1))
}

##' @rdname rxSolve
##' @export
solve.RxODE <- solve.rxSolve

.sharedPrint <- function(x, n, width, bound=""){
    .isDplyr <- requireNamespace("dplyr", quietly = TRUE) && RxODE.display.tbl;
    ## cat(sprintf("Dll: %s\n\n", rxDll(x)))
    df <- x$params.single
    pars.msg <- cli::rule(left=paste0(crayon::bold("Parameters"), " (",
                                      crayon::yellow(bound), crayon::bold$blue("$params"), "):"));
    if (!is.null(df)){
        cat(pars.msg, "\n");
        print(df)
    } else {
        df <- x$pars
        if (!is.null(df)){
            cat(pars.msg, "\n");
            if (rxIs(df, "data.frame")){
                if (!.isDplyr){
                    print(head(as.matrix(df), n = n));
                } else {
                    print(dplyr::as.tbl(df), n = n, width = width);
                }
            }
        }
    }
    df <- x$covs;
    if (!is.null(df)){
        cat(cli::rule(left=paste0(crayon::bold("Covariates"), " (",
                                  crayon::yellow(bound), crayon::bold$blue("$covs"), "):")), "\n");
        if (!.isDplyr){
            print(head(as.matrix(df), n = n));
        } else {
            print(dplyr::as.tbl(df), n = n, width = width);
        }
    }

    cat(cli::rule(left=paste0(crayon::bold("Initial Conditions"),
                              " (", crayon::yellow(bound), crayon::bold$blue("$inits"), "):")), "\n")
    print(x$inits);
    return(invisible(.isDplyr));
}

##' @author Matthew L.Fidler
##' @export
print.rxSolve <- function(x, ...){
    if (rxIs(x, "rxSolve")){
        bound <- .getBound(x, parent.frame(2));
        cat(cli::rule(center=crayon::bold("Solved RxODE object"), line="bar2"), "\n");
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
        .isDplyr <- .sharedPrint(x, n, width, bound)
        ## inits <- lst$inits[regexpr(regSens, names(lst$inits)) == -1];
        ## print(inits);
        cat(cli::rule(left=crayon::bold("First part of data (object):")), "\n");
        if (!.isDplyr){
            print(head(as.matrix(x), n = n));
        } else {
            print(dplyr::as.tbl(x), n = n, width = width);
        }
        message(cli::rule(line="bar2"))
    } else {
        print.data.frame(x)
    }
}

##' @author Matthew L.Fidler
##' @export
summary.rxSolve <- function(object, ...){
    if (rxIs(object, "rxSolve")){
        bound <- .getBound(object, parent.frame(2));
        cat(cli::rule(center=crayon::bold("Summary of Solved RxODE object"), line="bar2"), "\n");
        cat(cli::rule(left=paste0(crayon::bold("Model"),
                                  " (", crayon::yellow(bound), crayon::bold$blue("$model"), "):")), "\n");
        cat(rxNorm(object), "\n");
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
        .sharedPrint(object, n, width, bound)
        cat(cli::rule(left=crayon::bold("Summary of solved data:")), "\n");
        print(summary.data.frame(object))
        cat(cli::rule(line="bar2"), "\n")
    } else {
        class(object) <- "data.frame"
        NextMethod("summary", object);
    }
}

##' Check to see if this is an rxSolve object.
##'
##' @param x object to check to see if it is rxSolve
##'
##' If this is an rxSolve object that has expired strip all rxSolve
##' information.
##'
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
`[.rxSolve` <- function(x, i, j, drop){
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
    if (rxIs(new, "rx.event")){
        return(update(solved, events=new));
    } else {
        return(as.data.frame(solved) + new);
    }
}

##'@export
print.rxModelText <- function(x, ...){
    cat(cli::rule(center=crayon::bold("RxODE Model Syntax"), line="bar2"), "\n");
    cat(as.vector(x), "\n");
    cat(cli::rule(line="bar2"), "\n");
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
