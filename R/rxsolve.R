##'@rdname rxSolve
##'@export
rxControl <- function(scale = NULL,
                      method = c("liblsoda", "lsoda", "dop853"),
                      transitAbs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
                      maxsteps = 70000L, hmin = 0L, hmax = NA, hmaxSd= 0, hini = 0, maxordn = 12L, maxords = 5L, ...,
                      cores,
                      covsInterpolation = c("locf", "linear", "nocb", "midpoint"),
                      addCov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                      sigmaLower=-Inf, sigmaUpper=Inf,
                      nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress=10000L,
                      amountUnits = NA_character_, timeUnits = "hours", stiff,
                      theta = NULL,
                      thetaLower=-Inf, thetaUpper=Inf,
                      eta = NULL, addDosing=FALSE,
                      stateTrim=Inf, updateObject=FALSE,
                      omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                      omegaLower=-Inf, omegaUpper=Inf,
                      nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                      nStud = 1L, dfSub=0.0, dfObs=0.0, returnType=c("rxSolve", "matrix", "data.frame", "data.frame.TBS", "data.table", "tbl", "tibble"),
                      seed=NULL, nsim=NULL,
                      minSS=10L, maxSS=1000L,
                      infSSstep=12,
                      strictSS=TRUE,
                      params=NULL,events=NULL,
                      istateReset=TRUE,
                      subsetNonmem=TRUE,
                      linLog=FALSE,
                      maxAtolRtolFactor=0.1,
                      from=NULL,
                      to=NULL,
                      by=NULL,
                      length.out=NULL,
                      iCov=NULL,
                      keep=NULL,
                      drop=NULL,
                      idFactor=TRUE,
                      mxhnil=0,
                      hmxi=0.0,
                      warnIdSort=TRUE,
                      warnDrop=TRUE,
                      ssAtol = 1.0e-8, ssRtol = 1.0e-6,
                      safeZero=TRUE){
    .xtra <- list(...);
    if (is.null(transitAbs) && !is.null(.xtra$transit_abs)){
        transitAbs <- .xtra$transit_abs;
    }
    if (missing(updateObject) && !is.null(.xtra$update.object)){
        updateObject <- .xtra$update.object;
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
    ## stiff = TRUE, transitAbs = NULL,
    ## atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
    ## maxords = 5, ..., covsInterpolation = c("linear", "constant", "NOCB", "midpoint"),
    ## theta=numeric(), eta=numeric(), matrix=TRUE,addCov=FALSE,
    ## inC=FALSE, counts=NULL, doSolve=TRUE
    if (!missing(stiff) && missing(method)){
        if (rxIs(stiff, "logical")){
            if (stiff){
                method <- "lsoda"
                .Deprecated("method = \"lsoda\"",old="stiff=TRUE")
            } else {
                method <- "dop853"
                .Deprecated("method = \"dop853\"",old="stiff=FALSE")
            }
        }
    } else {
        if (!rxIs(method, "integer")){
            method <- match.arg(method);
        }
    }
    .matrixIdx <- c("rxSolve"=0, "matrix"=1, "data.frame"=2, "data.frame.TBS"=3, "data.table"=4,
                    "tbl"=5, "tibble"=5);
    if (!missing(returnType)){
        matrix <- .matrixIdx[match.arg(returnType)];
    } else if (!is.null(.xtra$return.type)){
        matrix <- .matrixIdx[.xtra$return.type];
    } else {
        matrix <- as.integer(matrix);
    }
    if (!rxIs(method, "integer")){
        .methodIdx <- c("lsoda"=1, "dop853"=0, "liblsoda"=2);
        method <- as.integer(.methodIdx[method]);
    }
    if (Sys.info()[["sysname"]] == "SunOS" && method == 2){
        method <- 1;
    }
    if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1];
    if (!rxIs(covsInterpolation, "integer")){
        covsInterpolation <- tolower(match.arg(covsInterpolation,
                                               c("linear", "locf", "LOCF", "constant",
                                                 "nocb", "NOCB", "midpoint")))
        if (covsInterpolation == "constant") covsInterpolation <- "locf";
        covsInterpolation  <- as.integer(which(covsInterpolation ==
                                               c("linear", "locf", "nocb", "midpoint")) - 1);
    }
    if (any(duplicated(names(.xtra)))){
        stop("Duplicate arguments do not make sense.");
    }
    if (any(names(.xtra)=="covs")){
        stop("Covariates can no longer be specified by 'covs' include them in the event dataset.");
    }
    if (missing(cores)){
        cores <- RxODE::rxCores();
    }
    .ret <- list(scale=scale,
                 method=method,
                 transitAbs=transitAbs,
                 atol=atol,
                 rtol=rtol,
                 maxsteps=maxsteps,
                 hmin=hmin,
                 hmax=hmax,
                 hini=hini,
                 maxordn=maxordn,
                 maxords=maxords,
                 covsInterpolation=covsInterpolation,
                 cores=cores,
                 addCov=addCov,
                 matrix=matrix,
                 sigma=lotri(sigma),
                 sigmaDf=sigmaDf,
                 nCoresRV=nCoresRV,
                 sigmaIsChol=sigmaIsChol,
                 nDisplayProgress=nDisplayProgress,
                 amountUnits=amountUnits,
                 timeUnits=timeUnits,
                 theta=theta,
                 eta=eta,
                 addDosing=addDosing,
                 stateTrim=stateTrim,
                 updateObject=updateObject,
                 omega=lotri(omega),
                 omegaDf=omegaDf,
                 omegaIsChol=omegaIsChol,
                 nSub=nSub,
                 thetaMat=thetaMat,
                 thetaDf=thetaDf,
                 thetaIsChol=thetaIsChol,
                 nStud=nStud,
                 dfSub=dfSub,
                 dfObs=dfObs,
                 seed=seed,
                 nsim=nsim,
                 minSS=minSS, maxSS=maxSS,
                 strictSS=as.integer(strictSS),
                 infSSstep=as.double(infSSstep),
                 istateReset=istateReset,
                 subsetNonmem=subsetNonmem,
                 linLog=linLog, hmaxSd=hmaxSd,
                 maxAtolRtolFactor=maxAtolRtolFactor,
                 from=from,
                 to=to,
                 by=by,
                 length.out=length.out,
                 iCov=iCov,
                 keep=keep, keepF=character(0), keepI=character(0),
                 drop=drop,
                 warnDrop=warnDrop,
                 omegaLower=omegaLower, omegaUpper=omegaUpper,
                 sigmaLower=sigmaLower, sigmaUpper=sigmaUpper,
                 thetaLower=thetaLower, thetaUpper=thetaUpper,
                 idFactor=idFactor,
                 mxhnil=mxhnil, hmxi=hmxi, warnIdSort=warnIdSort,
                 ssAtol=ssAtol, ssRtol = ssRtol, safeZero=as.integer(safeZero));
    return(.ret)
}

##' Solving \& Simulation of a ODE/solved system (and solving options) equation
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
##' @param iCov A data frame of individual non-time varying covariates
##'     to combine with the \code{params} to form a parameter
##'     data.frame.
##'
##' @param scale a numeric named vector with scaling for ode
##'     parameters of the system.  The names must correspond to the
##'     parameter identifiers in the ODE specification. Each of the
##'     ODE variables will be divided by the scaling factor.  For
##'     example \code{scale=c(center=2)} will divide the center ODE
##'     variable by 2.
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
##' @param hmax The maximum absolute step size allowed.  When
##'     \code{hmax=NA} (default), uses the average difference (+hmaxSd*sd) in times
##'     and sampling events. When \code{hmax=NULL} RxODE uses the
##'     maximum difference in times in your sampling and events.  The
##'     value 0 is equivalent to infinite maximum absolute step size.
##'
##' @param hmaxSd The number of standard deviations of the time
##'     difference to add to hmax. The default is 0
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
##' @param mxhnil maximum number of messages printed (per problem)
##'     warning that T + H = T on a step (H = step size).  This must
##'     be positive to result in a non-default value.  The default
##'     value is 0 (or infinite).
##'
##' @param hmxi inverse of the maximum absolute value of H to be used.
##'     hmxi = 0.0 is allowed and corresponds to an infinite hmax
##'     (default).  hmin and hmxi may be changed at any time, but will
##'     not take effect until the next change of H is considered.
##'     This option is only considered with method=liblsoda.
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
##' @param matrix A boolean indicating if a matrix should be returned
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
##'     For stiff ODE systems (\code{stiff = TRUE}), \code{RxODE} uses the
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
##' @param stateTrim When amounts/concentrations in one of the states
##'     are above this value, trim them to be this value. By default
##'     Inf.  Also trims to -stateTrim for large negative
##'     amounts/concentrations.  If you want to trim between a range
##'     say `c(0, 2000000)` you may specify 2 values with a lower and
##'     upper range to make sure all state values are in the
##'     reasonable range.
##'
##' @param updateObject This is an internally used flag to update the
##'     RxODE solved object (when supplying an RxODE solved object) as
##'     well as returning a new object.  You probably should not
##'     modify it's \code{FALSE} default unless you are willing to
##'     have unexpected results.
##'
##' @param returnType This tells what type of object is returned.  The currently supported types are:
##' \itemize{
##' \item \code{"rxSolve"} (default) will return a reactive data frame
##'      that can change easily change different pieces of the solve and
##'      update the data frame.  This is the currently standard solving
##'      method in RxODE,  is used for \code{rxSolve(object, ...)}, \code{solve(object,...)},
##' \item \code{"data.frame"} -- returns a plain, non-reactive data
##'      frame; Currently very slightly faster than \code{returnType="matrix"}
##' \item \code{"matrix"} -- returns a plain matrix with column names attached
##'     to the solved object.  This is what is used \code{object$run} as well as \code{object$solve}
##' \item \code{"data.table"} -- returns a \code{data.table}; The \code{data.table} is
##'     created by reference (ie \code{setDt()}), which should be fast.
##' \item \code{"tbl"} or \code{"tibble"} returns a tibble format.
##' }
##'
##' @param seed an object specifying if and how the random number
##'    generator should be initialized
##'
##' @param omega Estimate of Covariance matrix. When omega is a list,
##'     assume it is a block matrix and convert it to a full matrix
##'     for simulations.
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
##' @param minSS Minimum number of iterations for a steady-state dose
##'
##' @param maxSS Maximum number of iterations for a steady-state dose
##'
##' @param strictSS Boolean indicating if a strict steady-state is
##'     required. If a strict steady-state is (\code{TRUE}) required
##'     then at least \code{minSS} doses are administered and the
##'     total number of steady states doses will continue until
##'     \code{maxSS} is reached, or \code{atol} and \code{rtol} for
##'     every compartment have been reached.  However, if ODE solving
##'     problems occur after the \code{minSS} has been reached the
##'     whole subject is considered an invalid solve. If
##'     \code{strictSS} is \code{FALSE} then as long as \code{minSS}
##'     has been reached the last good solve before ODE solving
##'     problems occur is considered the steady state, even though
##'     either \code{atol}, \code{rtol} or \code{maxSS} have not
##'     been achieved.
##'
##' @param infSSstep Step size for determining if a constant infusion
##'     has reached steady state.  By default this is large value,
##'     420.
##'
##' @param ssAtol Steady state atol convergence factor.  Can be
##'     a vector based on each state.
##'
##' @param ssRtol Steady state rtol convergence factor.  Can be a
##'     vector based on each state.
##'
##' @param istateReset When \code{TRUE}, reset the \code{ISTATE} variable to 1 for
##'     lsoda and liblsoda with doses, like \code{deSolve}; When \code{FALSE}, do
##'     not reset the \code{ISTATE} variable with doses.
##'
##' @param addDosing Boolean indicating if the solve should add RxODE
##'     EVID and related columns.  This will also include dosing
##'     information and estimates at the doses.  Be default, RxODE
##'     only includes estimates at the observations. (default
##'     \code{FALSE}). When \code{addDosing} is \code{NULL}, only
##'     include \code{EVID=0} on solve and exclude any model-times or
##'     \code{EVID=2}. If \code{addDosing} is \code{NA} the classic
##'     \code{RxODE} EVID events. When \code{addDosing} is \code{TRUE}
##'     add the event information in NONMEM-style format; If
##'     \code{subsetNonmem=FALSE} RxODE will also extra event types
##'     (\code{EVID}) for ending infusion and modeled times:
##'
##' \itemize{
##'
##' \item \code{EVID=-1} when the modeled rate infusions are turned
##' off (matches \code{rate=-1})
##'
##' \item \code{EVID=-2} When the modeled duration infusions are
##' turned off (matches \code{rate=-2})
##'
##' \item \code{EVID=-10} When the specified \code{rate} infusions are
##' turned off (matches \code{rate>0})
##'
##' \item \code{EVID=-20} When the specified \code{dur} infusions are
##' turned off (matches \code{dur>0})
##'
##' \item \code{EVID=101,102,103,...} Modeled time where 101 is the
##' first model time, 102 is the second etc.
##'
##' }
##'
##' @param subsetNonmem subset to NONMEM compatible EVIDs only.  By default TRUE.
##'
##' @param linLog Boolean indicating if linear compartment models be
##'     calculated more accurately in the log-space (slower) By
##'     default this is off (\code{FALSE})
##'
##' @param maxAtolRtolFactor The maximum \code{atol}/\code{rtol} that
##'     FOCEi and other routines may adjust to.  By default 0.1
##'
##' @param from When there is no observations in the event table,
##'     start observations at this value. By default this is zero.
##'
##' @param to When there is no observations in the event table, end
##'     observations at this value. By default this is 24 + maximum
##'     dose time.
##'
##' @param length.out The number of observations to create if there
##'     isn't any observations in the event table. By default this is 200.
##'
##' @param by When there are no observations in the event table, this
##'     is the amount to increment for the observations between \code{from}
##'     and \code{to}.
##'
##' @param keep Columns to keep from either the input dataset or the
##'     \code{iCov} dataset.  With the \code{iCov} dataset, the column
##'     is kept once per line.  For the input dataset, if any records
##'     are added to the data LOCF (Last Observation Carried forward)
##'     imputation is performed.
##'
##' @param drop Columns to drop from the output
##'
##' @param idFactor This boolean indicates if original ID values
##'     should be maintained. This changes the default sequentially
##'     ordered ID to a factor with the original ID values in the
##'     original dataset.  By default this is enabled.
##'
##' @param warnIdSort Warn if the ID is not present and RxODE assumes
##'     the order of the parameters/iCov are the same as the order of
##'     the parameters in the input dataset.
##'
##' @param warnDrop Warn if column(s) were supposed to be dropped, but
##'     were not present.
##'
##' @param safeZero Use safe zero divide and log routines.  By default
##'     this is turned on but you may turn it off if you wish.
##'
##' @return An \dQuote{rxSolve} solve object that stores the solved
##'     value in a matrix with as many rows as there are sampled time
##'     points and as many columns as system variables (as defined by
##'     the ODEs and additional assignments in the RxODE model code).
##'     It also stores information about the call to allow dynamic
##'     updating of the solved object.
##'
##'     The operations for the object are similar to a data-frame, but
##'     expand the \code{$} and \code{[[""]]} access operators and
##'     assignment operators to resolve based on different parameter
##'     values, initial conditions, solver parameters, or events (by
##'     updating the \code{time} variable).
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
rxSolve.default <- function(object, params=NULL, events=NULL, inits = NULL, ...){
    on.exit({
        .clearPipe();
    });
    .applyParams <- FALSE
    .rxParams <- NULL
    if (rxIs(object, "rxEt")){
        if (!is.null(events)){
            stop("Events in pipeline AND in solving arguments, please provide just one.")
        }
        if (is.null(.pipelineRx)){
            stop("Need an RxODE compiled model as the start of the pipeline");
        } else {
            events <- object
            object <- .pipelineRx
        }
    } else if (rxIs(object, "rxParams")){
        .applyParams <- TRUE
        if (is.null(params) && !is.null(object$params)){
            params <- object$params;
        }
        if (is.null(.pipelineRx)){
            stop("Need an RxODE compiled model as the start of the pipeline");
        } else {
            .rxParams <- object
            object <- .pipelineRx
        }
        if (is.null(.pipelineEvents)){
            stop("Need an RxODE events as a part of the pipeline")
        } else {
            events <- .pipelineEvents;
            assignInMyNamespace(".pipelineEvents", NULL);
        }

    }
    if (!is.null(.pipelineEvents) && is.null(events) && is.null(params)){
        events <- .pipelineEvents;
    } else if (!is.null(.pipelineEvents) && !is.null(events)){
        stop("'events' in pipeline AND in solving arguments, please provide just one.")
    } else if (!is.null(.pipelineEvents) && !is.null(params) &&
               rxIs(params, "event.data.frame")){
        stop("'events' in pipeline AND in solving arguments, please provide just one.")
    }

    if (!is.null(.pipelineParams) && is.null(params)){
        params <- .pipelineParams;
    } else if (!is.null(.pipelineParams) && !is.null(params)){
        stop("'params' in pipeline AND in solving arguments, please provide just one.")
    }

    if (!is.null(.pipelineInits) && is.null(inits)){
        inits <- .pipelineInits;
    } else if (!is.null(.pipelineInits) && !is.null(inits)){
        stop("'inits' in pipeline AND in solving arguments, please provide just one.")
    }

    if (.applyParams){
        if (!is.null(.rxParams$inits)){
            inits <- .rxParams$inits
        }
    }
    .xtra <- list(...);
    if (any(duplicated(names(.xtra)))){
        stop("Duplicate arguments do not make sense.");
    }
    if (any(names(.xtra)=="covs")){
        stop("Covariates can no longer be specified by 'covs';\n  include them in the event dataset.\n\nIndividual covariates: Can be specified by a 'iCov' dataset\n each each individual covariate has a value\n\nTime varying covariates: modify input event data-frame or\n  eventTable to include covariates(https://tinyurl.com/y52wfc2y)\n\nEach approach needs the covariates named to match the variable in the model\n");
    }
    .nms <- names(as.list(match.call())[-1]);
    .lst <- list(...);
    .setupOnly <- 0L
    if (any(names(.lst)==".setupOnly")){
        .setupOnly <- .lst$.setupOnly;
    }
    .ctl <- rxControl(...,events=events,params=params);
    .n1 <- setdiff(intersect(tolower(names(params)),tolower(names(.ctl$iCov))),"id")
    .n2 <- c(.n1, setdiff(intersect(tolower(names(events)),tolower(names(.ctl$iCov))),"id"))
    .n1 <- unique(c(.n1, .n2))
    if (length(.n1) > 0){
        stop(sprintf("iCov has information contained in parameters/event data;\nDuplicate columns: %s", paste(.n1, collapse=", ")));
    }
    if (!is.null(.pipelineThetaMat) && is.null(.ctl$thetaMat)){
        .ctl$thetaMat <- .pipelineThetaMat;
    }
    if (!is.null(.pipelineOmega) && is.null(.ctl$omega)){
        .ctl$omega <- .pipelineOmega;
    }
    if (!is.null(.pipelineSigma) && is.null(.ctl$sigma)){
        .ctl$sigma <- .pipelineSigma;
    }
    if (!is.null(.pipelineDfObs) && .ctl$dfObs==0){
        .ctl$dfObs <- .pipelineDfObs;
    }
    if (!is.null(.pipelineDfSub) && .ctl$dfSub==0){
        .ctl$dfSub <- .pipelineDfSub;
    }
    if (!is.null(.pipelineNSub) && .ctl$nSub==1){
        .ctl$nSub <- .pipelineNSub;
    }
    if (!is.null(.pipelineNStud) && .ctl$nStud==1){
        .ctl$nStud <- .pipelineNStud;
    }
    if (!is.null(.pipelineICov) && is.null(.ctl$iCov)){
        .ctl$iCov <- .pipelineICov;
    }
    if (!is.null(.pipelineKeep) && is.null(.ctl$keep)){
        .ctl$keep <- .pipelineKeep;
    }
    if (.applyParams){
        if (!is.null(.rxParams$thetaMat) && is.null(.ctl$thetaMat)){
            .ctl$thetaMat <- .rxParams$thetaMat;
        }
        if (!is.null(.rxParams$omega) && is.null(.ctl$omega)){
            .ctl$omega <- .rxParams$omega;
        }
        if (!is.null(.rxParams$sigma) && is.null(.ctl$sigma)){
            .ctl$sigma <- .rxParams$sigma;
        }
        if (!is.null(.rxParams$dfSub)){
            if (.ctl$dfSub== 0){
                .ctl$dfSub <- .rxParams$dfSub;
            }
        }
        if (!is.null(.rxParams$nSub)){
            if (.ctl$nSub== 1){
                .ctl$nSub <- .rxParams$nSub;
            }
        }
        if (!is.null(.rxParams$nStud)){
            if (.ctl$nStud== 1){
                .ctl$nStud <- .rxParams$nStud;
            }
        }
        if (!is.null(.rxParams$dfObs)){
            if (.ctl$dfObs == 0){
                .ctl$dfObs <- .rxParams$dfObs;
            }
        }
        if (!is.null(.rxParams$iCov)){
            if (is.null(.ctl$iCov)){
                .ctl$iCov <- .rxParams$iCov;
            }
        }
        if (!is.null(.rxParams$keep)){
            if (is.null(.ctl$keep)){
                .ctl$keep <- .rxParams$keep;
            }
        }
    }
    if (.ctl$nSub==1 && inherits(.ctl$iCov, "data.frame")){
        .ctl$nSub <- length(.ctl$iCov[,1])
    } else if (.ctl$nSub !=1 && .ctl$nStud !=1 && inherits(.ctl$iCov, "data.frame")){
        if (.ctl$nSub !=length(.ctl$iCov[,1])){
            stop("'nSub' does not match the number of subjects in iCov");
        }
    } else if (.ctl$nSub !=1 && !.ctl$nStud !=1 && inherits(.ctl$iCov, "data.frame")){
        if (.ctl$nSub*.ctl$nStud !=length(.ctl$iCov[,1])){
            stop("'nSub'*'nStud' does not match the number of subjects in iCov");
        }
    }
    ## Prefers individual keep over keeping from the input data
    .keepI <- character(0)
    .keepF <- character(0)
    if (!is.null(.ctl$keep)){
        .mv <- rxModelVars(object);
        .vars <- c(.mv$lhs, .mv$state);
        .keepF <- setdiff(.ctl$keep, .vars)
        if (!is.null(.ctl$iCov)){
            .keepI <- intersect(.keepF, names(.ctl$iCov));
            .keepF <- setdiff(.keepF, .keepI);
        }
    }
    .ctl$keepI <- .keepI
    .ctl$keepF <- .keepF
    .ret <- rxSolve_(object, .ctl, .nms, .xtra,
                     params, events, inits,setupOnly=.setupOnly);
    if (.ctl$matrix == 4L){
        data.table::setDT(.ret);
    } else if (.ctl$matrix == 5L){
        .ret <- tibble::as_tibble(.ret);
    }
    return(.ret)
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

##' @rdname rxSolve
##' @export
predict.rxEt <- predict.RxODE

##' @rdname rxSolve
##' @export
predict.rxParams <- predict.RxODE

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
simulate.rxParams <- simulate.RxODE

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

##' @rdname rxSolve
##' @export
solve.rxParams <- solve.rxSolve

##' @rdname rxSolve
##' @export
solve.rxEt <- solve.rxSolve

.sharedPrint <- function(x, n, width, bound=""){
    ## nocov start
    .isDplyr <- requireNamespace("tibble", quietly = TRUE) && RxODE.display.tbl;
    ## cat(sprintf("Dll: %s\n\n", rxDll(x)))
    df <- x$params.single
    if (length(df) > 0){
        pars.msg <- .cliRule(left=paste0(crayon::bold("Parameters"), " (",
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
                        print(tibble::as_tibble(df), n = n, width = width);
                    }
                }
            }
        }
    }
    df <- x$covs;
    if (!is.null(df)){
        .cliRule(left=paste0(crayon::bold("Covariates"), " (",
                                  crayon::yellow(bound), crayon::bold$blue("$covs"), "):"))
        if (!.isDplyr){
            print(head(as.matrix(df), n = n));
        } else {
            print(tibble::as_tibble(df), n = n, width = width);
        }
    }

    .cliRule(left=paste0(crayon::bold("Initial Conditions"),
                              " (", crayon::yellow(bound), crayon::bold$blue("$inits"), "):"))
    print(x$inits);
    if (any(names(x) == "sim.id")){
        .uncert <- character(0)
        if (!is.null(x$thetaMat)){
            .uncert <- c(.uncert, paste0("parameters (", crayon::yellow(bound), crayon::bold$blue("$thetaMat"), " for changes)"))
        }
        if (!is.null(x$omegaList)){
            .uncert <- c(.uncert, paste0("omega matrix (", crayon::yellow(bound), crayon::bold$blue("$omegaList"), ")"))
        }
        if (!is.null(x$sigmaList)){
            .uncert <- c(.uncert, paste0("sigma matrix (", crayon::yellow(bound), crayon::bold$blue("$sigmaList"), ")"))
        }
        if (length(.uncert) == 0L){
            cat(paste0("\nSiulation ", crayon::bold("without uncertainty"), " in parameters, omega or sigma matricies\n\n"));
        } else if (length(.uncert) == 1L){
            cat(paste0("\nSimulation ", crayon::bold("with uncertainty"), " in ", paste(.uncert, collapse=", "), "\n\n"));
        } else {
            cat(paste0("\nSimulation ", crayon::bold("with uncertainty"), " in:\n  - ", paste(.uncert, collapse="\n  - "), "\n\n"));
        }
    }
    return(invisible(.isDplyr));
    ## nocov end
}

##' @author Matthew L.Fidler
##' @export
print.rxSolve <- function(x, ...){
    ##nocov start
    if (rxIs(x, "rxSolve")){
        bound <- .getBound(x, parent.frame(2));
        .cliRule(center=crayon::bold("Solved RxODE object"), line="bar2")
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
        .cliRule(left=crayon::bold("First part of data (object):"))
        if (!.isDplyr){
            print(head(as.matrix(x), n = n));
        } else {
            print(tibble::as_tibble(x), n = n, width = width);
        }
        .cliRule(line="bar2")
    } else {
        print.data.frame(x)
    }
    ##nocov end
}

##' @author Matthew L.Fidler
##' @export
summary.rxSolve <- function(object, ...){
    ## nocov start
    if (rxIs(object, "rxSolve")){
        bound <- .getBound(object, parent.frame(2));
        .cliRule(center=crayon::bold("Summary of Solved RxODE object"), line="bar2")
        .cliRule(left=paste0(crayon::bold("Model"),
                                  " (", crayon::yellow(bound), crayon::bold$blue("$model"), "):"));
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
        .cliRule(left=crayon::bold("Summary of solved data:"))
        print(summary.data.frame(object))
        .cliRule(line="bar2")
    } else {
        class(object) <- "data.frame"
        NextMethod("summary", object);
    }
    ## nocov end
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
##' @param new New information added to the table.
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
    ## nocov start
    .cliRule(center=crayon::bold("RxODE Model Syntax"), line="bar2");
    .code <- deparse(body(eval(parse(text=paste("function(){",as.vector(x),"}")))))
    .code[1]  <- "RxODE({"
    .code[length(.code)]  <- "})";
    cat(paste(.code,collapse="\n"), "\n");
    .cliRule(line="bar2")
    ## nocov end
}

##'@export
drop_units.rxSolve <- function(x){
    dropUnitsRxSolve(x);
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


.SD <- NULL
`:=` <- function (...) {
    stop("This is only used in data.table")
}


##'@export
confint.rxSolve <- function(object, parm=NULL, level = 0.95, ...){
    RxODE::rxReq("data.table")
    ## p1 <-eff <-Percentile <-sim.id <-id <-p2 <-p50 <-p05 <- p95 <- . <- time <- trt <- NULL
    ## RxODE::rxReq("dplyr")
    ## RxODE::rxReq("tidyr")
    if (level <=0 || level >=1){
        stop("simulation summaries must be between 0 and 1");
    }
    .stk <- rxStack(object, parm)
    data.table::setDT(.stk)
    .a <- (1-level)/2;
    .p <- c(.a, 0.5, 1-.a);
    .lst <- list(lvl=paste0("p",.p*100),
                 parm=levels(.stk$trt));
    class(.lst) <- "rxHidden";
    if (object$env$args$nStud <= 1){
        if (object$env$args$nSub < 2500){
            warning("In order to put confidence bands around the intervals, you need at least 2500 simulations.")
            message("summarizing data...", appendLF=FALSE)
            .stk <- .stk[, list(p1=.p, eff=stats::quantile(.SD$value, probs=.p, na.rm=TRUE), Percentile=sprintf("%s%%",.p*100)),
                 by=c("time", "trt")]
            if (requireNamespace("tibble", quietly = TRUE)){
                .stk <- tibble::as_tibble(.stk)
            }
            .cls <- c("rxSolveConfint1", class(.stk));
            attr(.cls, ".rx") <- .lst
            class(.stk) <- .cls
            message("done.")
            return(.stk)
        } else {
            .n <- round(sqrt(object$env$args$nSub));
        }
    } else {
        .n <- object$env$args$nStud;
    }
    message("Summarizing data")
    .ret <- .stk[, id := sim.id %% .n
         ][, list(p1=.p, eff=stats::quantile(.SD$value, probs=.p, na.rm=TRUE)), by = c("id", "time", "trt")
           ][,setNames(as.list(stats::quantile(.SD$eff, probs=.p, na.rm=TRUE)),
                       sprintf("p%s",.p*100)),
             by = c("p1", "time", "trt")
             ]
    .ret$Percentile <- factor(sprintf("%s%%",.ret$p1*100))
    if (requireNamespace("tibble", quietly = TRUE)){
        .ret <- tibble::as_tibble(.ret)
    }
    message("done.")
    .cls <- c("rxSolveConfint2", class(.ret));
    attr(.cls, ".rx") <- .lst
    class(.ret) <- .cls
    return(.ret);
}
