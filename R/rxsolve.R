## nlmixr-style parsing of omega matrix expressions.
##' Easily Parse block-diagonal matrices with lower triangular info
##'
##' @param x list, matrix or expression, see details
##'
##' @param ... Other arguments treated as a list that will be
##'     concatenated then reapplied to this function.
##'
##' @return named symmetric matrix useful in RxODE simulations (and
##'     perhaps elsewhere)
##'
##' @details
##'
##'  This can take an R matrix, a list including matrices or expressions, or expressions
##'
##'  Expressions can take the form
##'
##'  name ~ estimate
##'
##'  Or the lower triangular matrix when "adding" the names
##'
##'  name1 + name2 ~ c(est1,
##'                    est2, est3)
##'
##'  The matricies are concatenated into a block diagonal matrix, like
##'  \code{\link[Matrix]{bdiag}}, but allows expressions to specify
##'  matrices easier.
##'
##'
##' @examples
##'
##' ## A few ways to specify the same matrix
##' rxMatrix({et2 + et3 + et4 ~ c(40,
##'                                0.1, 20,
##'                                0.1, 0.1, 30)})
##'
##' ## You  do not need to enclose in {}
##' rxMatrix(et2 + et3 + et4 ~ c(40,
##'                                0.1, 20,
##'                                0.1, 0.1, 30),
##'           et5 ~ 6)
##' ## But if you do enclose in {}, you can use multi-line matrix specifications:
##'
##' rxMatrix({et2 + et3 + et4 ~ c(40,
##'                                0.1, 20,
##'                                0.1, 0.1, 30);
##'           et5 ~ 6;
##'           })
##'
##' ## You can also add lists or actual R matrices as in this example:
##' rxMatrix(list(et2 + et3 + et4 ~ c(40,
##'                                   0.1, 20,
##'                                   0.1, 0.1, 30),
##'               matrix(1,dimnames=list("et5","et5"))))
##'
##' ## Overall this is a flexible way to specify symmetric block diagonal matrices.
##'
##' @author Matthew L Fidler
##' @seealso \code{\link{rxSolve}}
##' @export
rxMatrix  <- function(x, ...){
    .lst  <- list(...);
    if (is.null(x)){
        .ret  <- NULL;
    } else if (is.list(x)){
        omega  <- lapply(x, rxMatrix);
        if (is(omega, "list")){
            .omega <- as.matrix(Matrix::bdiag(omega));
            .d <- unlist(lapply(seq_along(omega), function(x){dimnames(omega[[x]])[2]}))
            dimnames(.omega) <- list(.d, .d);
            omega <- .omega;
        }
        .ret  <- omega
    } else if (is.matrix(x)) {
        .ret  <- x
    } else {
        .env  <- new.env(parent=emptyenv());
        .env$df  <- NULL;
        .env$eta1 <- 0L;
        .f  <- function(x, env){
            if (is.name(x)){
                return(character())
            } else if (is.call(x)){
                if (identical(x[[1]], quote(`~`))){
                    if (length(x[[3]]) == 1){
                        ## et1 ~ 0.2
                        env$netas <- 1;
                        env$eta1 <- env$eta1 + 1;
                        env$names  <- c(env$names, as.character(x[[2]]));
                        env$df  <- rbind(env$df,
                                         data.frame(i=env$eta1, j=env$eta1, x=as.numeric(eval(x[[3]]))));
                    } else {
                        ## et1+et2+et3~c() lower triangular matrix
                        ## Should fixed be allowed????
                        if (any(tolower(as.character(x[[3]][[1]])) == c("c", "fix", "fixed"))){
                            if (any(tolower(as.character(x[[3]][[1]])) == c("fix", "fixed"))){
                                stop("fix/fixed are not allowed with RxODE omega specifications.");
                            }
                            env$netas <- length(x[[3]]) - 1;
                            .num <- sqrt(1+env$netas*8)/2-1/2
                            if (round(.num) == .num){
                                .n <- unlist(strsplit(as.character(x[[2]]), " +[+] +"));
                                .n <- .n[.n != "+"];
                                if(length(.n) == .num){
                                    env$names  <- c(env$names, .n);
                                    .r <- x[[3]][-1];
                                    .r <- sapply(.r, function(x){
                                        return(as.numeric(eval(x)));
                                    });
                                    .i <- 0
                                    .j <- 1;
                                    for (.k in seq_along(.r)){
                                        .v <- .r[.k];
                                        .i <- .i + 1;
                                        if (.i==.j){
                                            env$df  <- rbind(env$df,
                                                             data.frame(i=env$eta1+.i, j=env$eta1+.i, x=.v));
                                            .j <- .j + 1;
                                            .i <- 0;
                                        } else {
                                            env$df  <- rbind(env$df,
                                                             data.frame(i=c(env$eta1+.i, env$eta1+.j),
                                                                        j=c(env$eta1+.j, env$eta1+.i), x=.v));
                                        }
                                    }
                                    env$eta1 <- env$eta1 + .num;
                                }  else {
                                    stop("The left handed side of the expression must match the number of items in the lower triangular matrix.");
                                }
                            }
                        } else {
                            stop("Omega expression should be 'name ~ c(lower-tri)'");
                        }
                    }
                } else if (identical(x[[1]], quote(`{`))){
                    lapply(x, .f, env=env)
                } else if (identical(x[[1]], quote(`quote`))){
                    lapply(x[[2]], .f, env=env)
                } else {
                    stop("Omega expression should be 'name ~ c(lower-tri)'")
                }
            } else if (is.pairlist(x)) {
                lapply(x, .f, env=env);
            } else if (is.atomic(x)){
                stop("Something wrong with matrix parsing, likely a number is on a line without an identifier.");
            } else {
                stop("Unknown element in matrix parsing")
            }
        }
        .sX  <- substitute(x)
        if (is.call(.sX)){
            if (identical(.sX[[1]], quote(`[[`))){
                .sX  <- x
            }
        }
        .doParse  <- TRUE
        if (is.call(.sX)){
            if (identical(.sX[[1]], quote(`matrix`))){
                return(x);
            } else if (identical(.sX[[1]], quote(`list`))){
                omega  <- lapply(x, rxMatrix);
                if (is(omega, "list")){
                    .omega <- as.matrix(Matrix::bdiag(omega));
                    .d <- unlist(lapply(seq_along(omega), function(x){dimnames(omega[[x]])[2]}))
                    dimnames(.omega) <- list(.d, .d);
                    omega <- .omega;
                }
                .ret  <- omega
                .doParse  <- FALSE
            }
        }
        if (.doParse){
            .f(.sX,.env);
            .ret <- diag(.env$eta1);
            for (.i in seq_along(.env$df$i)){
                .ret[.env$df$i[.i], .env$df$j[.i]]  <- .env$df$x[.i];
            }
            dimnames(.ret)  <- list(.env$names,.env$names)
        }
    }
    if (length(.lst)==0) return(.ret)
    else return(rxMatrix(c(list(.ret),.lst)))
}

##'@rdname rxSolve
##'@export
rxControl <- function(scale = NULL,
                      method = c("liblsoda", "lsoda", "dop853"),
                      transitAbs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
                      maxsteps = 70000L, hmin = 0L, hmax = NA, hmaxSd= 0, hini = 0, maxordn = 12L, maxords = 5L, ...,
                      cores,
                      covsInterpolation = c("locf", "linear", "nocb", "midpoint"),
                      addCov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                      nCoresRV = 1L, sigmaIsChol = FALSE, nDisplayProgress=10000L,
                      amountUnits = NA_character_, timeUnits = "hours", stiff,
                      theta = NULL, eta = NULL, addDosing=FALSE,
                      stateTrim=Inf, updateObject=FALSE,
                      omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                      nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                      nStud = 1L, dfSub=0.0, dfObs=0.0, returnType=c("rxSolve", "matrix", "data.frame", "data.frame.TBS"),
                      seed=NULL, nsim=NULL,
                      minSS=7, maxSS=70000,
                      atolSS=1e-9, rtolSS=1e-9,
                      params=NULL,events=NULL,
                      istateReset=TRUE,
                      subsetNonmem=TRUE,
                      linLog=FALSE){
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
                warning("stiff=TRUE has been replaced with method = \"lsoda\".")
            } else {
                method <- "dop853"
                warning("stiff=FALSE has been replaced with method = \"dop853\".")
            }
        }
    } else {
        if (!rxIs(method, "integer")){
            method <- match.arg(method);
        }
    }
    .matrixIdx <- c("rxSolve"=0, "matrix"=1, "data.frame"=2, "data.frame.TBS"=3);
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
                 sigma=rxMatrix(sigma),
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
                 omega=rxMatrix(omega),
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
                 atolSS=atolSS[1], rtolSS=rtolSS[1],
                 istateReset=istateReset,
                 subsetNonmem=subsetNonmem,
                 linLog=linLog, hmaxSd=hmaxSd);
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
##' @param scale a numeric named vector with scaling for ode
##'     parameters of the system.  The names must correstond to the
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
##'
##' @param stateTrim When amounts/concentrations in one of the states
##'     are above this value, trim them to be this value. By default
##'     Inf.  Also trims to -stateTrim for lage negative
##'     amounts/concentrations
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
##'      frame; Currently very slightly Faster than \code{returnType="matrix"}
##' \item \code{"matrix"} -- returns a plain matrix with column names attached
##'     to the solved object.  This is what is used \code{object$run} as well as \code{object$solve}
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
##' @param atolSS Absolute tolerance to check if a solution arrived at
##'     steady state.
##'
##' @param rtolSS Relative tolerance to check if a solution arrived at
##'     steady state.
##'
##' @param istateReset When TRUE, reset the ISTATE variable to 1 for
##'     lsoda and liblsoda with doses, like deSolve; When FALSE, do
##'     not reset the ISTATE variable with doses.
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
rxSolve.default <- function(object, params=NULL, events=NULL, inits = NULL, ...){
    on.exit({rxSolveFree()});
    .xtra <- list(...);
    ## if (is.null(transitAbs) && !is.null(.xtra$transit_abs)){
    ##     transitAbs <- .xtra$transit_abs;
    ## }
    ## if (missing(updateObject) && !is.null(.xtra$update.object)){
    ##     updateObject <- .xtra$update.object;
    ## }
    ## if (missing(covsInterpolation) && !is.null(.xtra$covs_interpolation)){
    ##     covsInterpolation <- .xtra$covs_interpolation;
    ## }
    ## if (missing(addCov) && !is.null(.xtra$add.cov)){
    ##     addCov <- .xtra$add.cov;
    ## }
    ## if (!is.null(seed)){
    ##     set.seed(seed);
    ## }
    ## if (!is.null(nsim)){
    ##     if (rxIs(params, "eventTable") || rxIs(events, "eventTable") && nSub == 1L){
    ##         nSub <- nsim;
    ##     } else if (nStud == 1L){
    ##         nStud <- nsim;
    ##     }

    ## }
    ## ## stiff = TRUE, transitAbs = NULL,
    ## ## atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 50000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
    ## ## maxords = 5, ..., covsInterpolation = c("linear", "constant", "NOCB", "midpoint"),
    ## ## theta=numeric(), eta=numeric(), matrix=TRUE,addCov=FALSE,
    ## ## inC=FALSE, counts=NULL, doSolve=TRUE
    ## if (!missing(stiff) && missing(method)){
    ##     if (rxIs(stiff, "logical")){
    ##         if (stiff){
    ##             method <- "lsoda"
    ##             warning("stiff=TRUE has been replaced with method = \"lsoda\".")
    ##         } else {
    ##             method <- "dop853"
    ##             warning("stiff=FALSE has been replaced with method = \"dop853\".")
    ##         }
    ##     }
    ## } else {
    ##     if (!rxIs(method, "integer")){
    ##         method <- match.arg(method);
    ##     }
    ## }
    ## .matrixIdx <- c("rxSolve"=0, "matrix"=1, "data.frame"=2, "data.frame.TBS"=3);
    ## if (!missing(returnType)){
    ##     matrix <- .matrixIdx[match.arg(returnType)];
    ## } else if (!is.null(.xtra$return.type)){
    ##     matrix <- .matrixIdx[.xtra$return.type];
    ## } else {
    ##     matrix <- as.integer(matrix);
    ## }
    ## if (!rxIs(method, "integer")){
    ##     .methodIdx <- c("lsoda"=1, "dop853"=0, "liblsoda"=2);
    ##     method <- as.integer(.methodIdx[method]);
    ## }
    ## if (Sys.info()[["sysname"]] == "SunOS" && method == 2){
    ##     method <- 1;
    ## }
    ## if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1];
    ## covsInterpolation <- tolower(match.arg(covsInterpolation,
    ##                                        c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")))
    ## if (covsInterpolation == "constant") covsInterpolation <- "locf";
    ## covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
    if (any(duplicated(names(.xtra)))){
        stop("Duplicate arguments do not make sense.");
    }
    if (any(names(.xtra)=="covs")){
        stop("Covariates can no longer be specified by 'covs' include them in the event dataset.");
    }
    .nms <- names(as.list(match.call())[-1]);
    ## if (rxIs(object, "rxSolve")){
    ##     if (rxIs(params, "rx.event") && is.null(events)){
    ##         .nms <- c("events",.nms);
    ##         events <- params
    ##         params <- NULL;
    ##     }
    ## }
    .lst <- list(...);
    .setupOnly <- 0L
    if (any(names(.lst)==".setupOnly")){
        .setupOnly <- .lst$.setupOnly;
    }
    rxSolve_(object, rxControl(...,events=events,params=params), .nms, .xtra,
             params, events, inits,setupOnly=.setupOnly);
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
    if (any(names(x) == "sim.id")){
        .uncert <- character(0)
        if (!is.null(x$thetaMat)){
            .uncert <- c(.uncert, paste0("parameters (", crayon::yellow(bound), crayon::bold$blue("$thetaMat"), " for changes)"))
        }
        if (!is.null(x$omegaList)){
            .uncert <- c(.uncert, paste0("omega matrix (", crayon::yellow(bound), crayon::bold$blue("$omegaList"), ")"))
        }
        if (!is.null(x$omegaList)){
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
    .code <- deparse(body(eval(parse(text=paste("function(){",as.vector(x),"}")))))
    .code[1]  <- "RxODE({"
    .code[length(.code)]  <- "})";
    cat(paste(.code,collapse="\n"), "\n");
    cat(cli::rule(line="bar2"), "\n");
}

##'@export
plot.rxSolve <- function(x,y,...){
    .cmts <- c(as.character(substitute(y)),
               names(sapply(as.character(as.list(match.call()[-(1:3)])),`c`)))
    if (length(.cmts)==1 &&.cmts[1]==""){
        .cmts <- NULL
    }
    .dat <- rxStack(x,.cmts);
    time <- value <- id <- sim.id  <- NULL
    if (any(names(.dat)=="id")){
        .dat$id <- factor(.dat$id);
        if (length(.cmts)==1){
            ggplot(.dat,ggplot2::aes(time,value,color=id))+
                geom_line(size=1.2) + ylab(.cmts)
        } else {
            ggplot(.dat,ggplot2::aes(time,value,color=id))+
                geom_line(size=1.2) +facet_wrap( ~ trt, scales="free_y")
        }
    } else if (any(names(.dat)=="sim.id")){
        .dat$sim.id <- factor(.dat$sim.id);
        if (length(.cmts)==1){
            ggplot(.dat,ggplot2::aes(time,value,color=sim.id))+
                geom_line(size=1.2) + ylab(.cmts)
        } else {
            ggplot(.dat,ggplot2::aes(time,value,color=sim.id))+
                geom_line(size=1.2) +facet_wrap( ~ trt, scales="free_y")
        }
    } else {
        if (length(.cmts)==1){
            ggplot(.dat,ggplot2::aes(time,value))+
                geom_line(size=1.2) + ylab(.cmts)
        } else {
            ggplot(.dat,ggplot2::aes(time,value))+
                geom_line(size=1.2) +facet_wrap(~ trt, scales="free_y")
        }
    }
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
