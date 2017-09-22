rex::register_shortcuts("RxODE");
##' Create an ODE-based model specification
##'
##' Create a dynamic ODE-based model object suitably for translation
##' into fast C code
##'
##' @param model This is the ODE model specification.  It can be:
##' \itemize{
##'  \item a string containing the set of ordinary differential
##'     equations (ODE) and other expressions defining the changes in
##'     the dynamic system.
##'  \item a file name where the ODE system equation is contained
##'  \item An ODE expression enclosed in \code{\{\}}
##'   }
##'
##' (see also the \code{filename} argument). For
##'     details, see the sections \dQuote{Details} and
##'     \dQuote{\code{RxODE Syntax}} below.
##'
##' @param modName a string to be used as the model name. This string
##'     is used for naming various aspects of the computations,
##'     including generating C symbol names, dynamic libraries,
##'     etc. Therefore, it is necessary that \code{modName} consists of
##'     simple ASCII alphanumeric characters starting with a letter.
##'
##' @param wd character string with a working directory where to
##'     create a subdirectory according to \code{modName}. When
##'     specified, a subdirectoy named after the
##'     \dQuote{\code{modName.d}} will be created and populated with a
##'     C file, a dynamic loading library, plus various other working
##'     files. If missing, the files are created (and removed) in the
##'     temporary directory, and the RxODE DLL for the model is
##'     created in the current directory named \code{rx_????_platform}, for
##'     example \code{rx_129f8f97fb94a87ca49ca8dafe691e1e_i386.dll}
##'
##' @param filename A file name or connection object where the
##'     ODE-based model specification resides. Only one of \code{model}
##'     or \code{filename} may be specified.
##'
##' @param do.compile logical specifying whether to carry out the
##'     parsing of the ODE into C code and its compilation. Default is
##'     \code{TRUE}
##'
##' @param extraC  Extra c code to include in the model.  This can be
##'     useful to specify functions in the model.  These C functions
##'     should usually take \code{double} precision arguments, and
##'     return \code{double} precision values.
##'
##' @param debug is a boolean indicating if the executable should be
##'     compiled with verbose debugging information turned on.
##'
##' @param calcSens boolean indicating if RxODE will calculate the
##'     sennsitivities according to the specified ODEs.
##'
##' @param calcJac boolean indicating if RxODE will calculate the
##'     Jacobain according to the specified ODEs.
##'
##' @param collapseModel boolean indicating if RxODE will remove all
##'     LHS variables when calculating sensitivites.
##'
##' @param ... any other arguments are passed to the function
##'     \code{\link{readLines}}, (e.g., encoding).
##'
##' The \dQuote{Rx} in the name \code{RxODE} is meant to suggest the
##' abbreviation \emph{Rx} for a medical prescription, and thus to
##' suggest the package emphasis on pharmacometrics modeling, including
##' pharmacokinetics (PK), pharmacodynamics (PD), disease progression,
##' drug-disease modeling, etc.
##'
##' The ODE-based model specification may be coded inside a character
##' string or in a text file, see Section \emph{RxODE Syntax} below for
##' coding details.  An internal \code{RxODE} compilation manager
##' object translates the ODE system into C, compiles it, and
##' dynamically loads the object code into the current R session.  The
##' call to \code{RxODE} produces an object of class \code{RxODE} which
##' consists of a list-like structure (closure) with various member
##' functions (see Section \emph{Value} below).
##'
##' For evaluating \code{RxODE} models, two types of inputs may be
##' provided: a required set of time points for querying the state of
##' the ODE system and an optional set of doses (input amounts).  These
##' inputs are combined into a single \emph{event table} object created
##' with the function \code{\link{eventTable}}.
##'
##' @section RxODE Syntax:
##'
##' An \code{RxODE} model specification consists of one or more
##' statements terminated by semi-colons, \sQuote{\code{;}}, and
##' optional comments (comments are delimited by \code{#} and an
##' end-of-line marker).  \strong{NB:} Comments are not allowed inside
##' statements.
##'
##' A block of statements is a set of statements delimeted by curly
##' braces, \sQuote{\code{\{ ... \}}}. Statements can be either
##' assignments or conditional \code{if} statements. Assignment
##' statements can be: (1) \dQuote{simple} assignmets, where the left
##' hand is an identifier (i.e., variable), (2) special
##' \dQuote{time-derivative} assignments, where the left hand specifies
##' the change of that variable with respect to time e.g.,
##' \code{d/dt(depot)}, or (3) special \dQuote{jacobian} assignments,
##' where the left hand specifies the change of of the ODE with respect
##' to one of the parameters, e.g. \code{df(depot)/dy(kel)}.  The
##' \dQuote{jacobian} assignments are not required, and are only useful
##' for very stiff differential systems.
##'
##' Expressions in assignment and \sQuote{\code{if}} statements can be
##' numeric or logical (no character expressions are currently
##' supported). Numeric expressions can include the following numeric
##' operators (\sQuote{\code{+}}, \sQuote{\code{-}}, \sQuote{\code{*}},
##' \sQuote{\code{/}}, \sQuote{\code{^}}), and those mathematical
##' functions defined in the C or the R math libraries (e.g.,
##' \code{fabs}, \code{exp}, \code{log}, \code{sin}).  (Notice that the
##' modulo operator \sQuote{\code{\%}} is currently not supported.)
##'
##' Identifiers in an \code{RxODE} model specification can refer to:
##' \itemize{
##'    \item state variables in the dynamic system (e.g., compartments in a
##'          pharmacokinetics/pharamcodynamics model);
##'    \item implied input variable, \code{t} (time),
##'    \code{podo} (oral dose, for absorption models), and
##'    \code{tlast} (last time point);
##'    \item model parameters, (\code{ka} rate of absorption, \code{CL}
##'        clearance, etc.);
##'    \item \code{pi}, for the constant pi.
##'    \item others, as created by assignments as part of the model
##'          specification.
##' }
##'
##' Identifiers consists of case-sensitive alphanumeric characters,
##' plus the underscore \sQuote{_} character.  \strong{NB:} the dot
##' \sQuote{.} character is \strong{not} a valid character identifier.
##'
##' The values of these variables at pre-specified time points are
##' saved as part of the fitted/integrated/solved model (see
##' \code{\link{eventTable}}, in particular its member function
##' \code{add.sampling} that defines a set of time points at which to
##' capture a snapshot of the syste via the values of these variables).
##'
##' The ODE specification mini-language is parsed with the help of the
##' open source tool \emph{DParser}, Plevyak (2015).
##'
##'
##' @return An object (closure) of class \dQuote{\code{RxODE}} (see Chambers and Temple Lang (2001))
##'      consisting of the following list of strings and functions:
##'
##'     \item{modName}{the name of the model (a copy of the input argument).}
##'     \item{model}{a character string holding the source model specification.}
##'     \item{get.modelVars}{a function that returns a list with 3 character
##'         vectors, \code{params}, \code{state}, and \code{lhs} of variable names used in the model
##'         specification. These will be output when the model is computed (i.e., the ODE solved by integration).}
##'
##'       \item{solve}{this function solves (integrates) the ODE. This
##'           is done by passing the code to \code{\link{rxSolve}}.
##'           This is as if you called \code{rxSolve(RxODEobject, ...)},
##'           but returns a matrix instead of a rxSolve object.
##'
##'           \code{params}: a numeric named vector with values for every parameter
##'           in the ODE system; the names must correspond to the parameter
##'           identifiers used in the ODE specification;
##'
##'           \code{events}: an \code{eventTable} object describing the
##'           input (e.g., doses) to the dynamic system and observation
##'           sampling time points (see  \code{\link{eventTable}});
##'
##'           \code{inits}: a vector of initial values of the state variables
##'           (e.g., amounts in each compartment), and the order in this vector
##'           must be the same as the state variables (e.g., PK/PD compartments);
##'
##'
##'           \code{stiff}: a logical (\code{TRUE} by default) indicating whether
##'           the ODE system is stifff or not.
##'
##'           For stiff ODE sytems (\code{stiff = TRUE}), \code{RxODE} uses
##'           the LSODA (Livermore Solver for Ordinary Differential Equations)
##'           Fortran package, which implements an automatic method switching
##'           for stiff and non-stiff problems along the integration interval,
##'           authored by Hindmarsh and Petzold (2003).
##'
##'           For non-stiff systems (\code{stiff = FALSE}), \code{RxODE} uses DOP853,
##'           an explicit Runge-Kutta method of order 8(5, 3) of Dormand and Prince
##'           as implemented in C by Hairer and Wanner (1993).
##'
##'           \code{trans_abs}: a logical (\code{FALSE} by default) indicating
##'           whether to fit a transit absorption term
##'           (TODO: need further documentation and example);
##'
##'           \code{atol}: a numeric absolute tolerance (1e-08 by default);
##'
##'           \code{rtol}: a numeric relative tolerance (1e-06 by default).e
##'
##'           The output of \dQuote{solve} is a matrix with as many rows as there
##'           are sampled time points and as many columns as system variables
##'           (as defined by the ODEs and additional assigments in the RxODE model
##'               code).}
##'
##'       \item{isValid}{a function that (naively) checks for model validity,
##'           namely that the C object code reflects the latest model
##'           specification.}
##'       \item{version}{a string with the version of the \code{RxODE}
##'           object (not the package).}
##'       \item{dynLoad}{a function with one \code{force = FALSE} argument
##'           that dynamically loads the object code if needed.}
##'       \item{dynUnload}{a function with no argument that unloads
##'           the model object code.}
##'       \item{cmpMgr}{a \dQuote{compilation manager} object, see
##'           \code{\link{rx.initCmpMgr}}.}
##'       \item{delete}{removes all created model files, including C and DDL files.
##'           The model object is no longer valid and should be removed, e.g.,
##'           \code{rm(m1)}.}
##'       \item{run}{deprecated, use \code{solve}.}
##'       \item{parse}{deprecated.}
##'       \item{compile}{deprecated.}
##'       \item{get.index}{deprecated.}
##'       \item{getObj}{internal (not user callable) function.}
##'
##' @references
##'
##' Chamber, J. M. and Temple Lang, D. (2001)
##' \emph{Object Oriented Programming in R}.
##' R News, Vol. 1, No. 3, September 2001.
##' \url{https://cran.r-project.org/doc/Rnews/Rnews_2001-3.pdf}.
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
##' Plevyek, J.
##' \emph{Dparser}, \url{http://dparser.sourceforge.net}. Web. 12 Oct. 2015.
##'
##' @author Melissa Hallow, Wenping Wang and Matthew Fidler
##'
##' @seealso \code{\link{eventTable}}
##'
##' @examples
##' # Step 1 - Create a model specification
##' ode <- "
##'    # A 4-compartment model, 3 PK and a PD (effect) compartment
##'    # (notice state variable names 'depot', 'centr', 'peri', 'eff')
##'
##'    C2 = centr/V2;
##'    C3 = peri/V3;
##'    d/dt(depot) =-KA*depot;
##'    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
##'    d/dt(peri)  =                    Q*C2 - Q*C3;
##'    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
##' "
##'
##' m1 <- RxODE(model = ode, modName = "m1")
##' print(m1)
##'
##' # Step 2 - Create the model input as an EventTable,
                                        ##' # including dosing and observation (sampling) events
##'
##' # QD (once daily) dosing for 5 days.
##'
##' qd <- eventTable(amount.units = "ug", time.units = "hours")
##' qd$add.dosing(dose = 10000, nbr.doses = 5, dosing.interval = 24)
##'
##' # Sample the system hourly during the first day, every 8 hours
##' # then after
##'
##' qd$add.sampling(0:24)
##' qd$add.sampling(seq(from = 24+8, to = 5*24, by = 8))
##'
##' # Step 3 - set starting parameter estimates and initial
##' # values of the state
##'
##' theta <-
##'     c(KA = .291, CL = 18.6,
##'       V2 = 40.2, Q = 10.5, V3 = 297.0,
##'       Kin = 1.0, Kout = 1.0, EC50 = 200.0)
##'
##' # init state variable
##' inits <- c(0, 0, 0, 1);
##'
##' # Step 4 - Fit the model to the data
##'
##' qd.cp <- m1$solve(theta, events = qd, inits)
##'
##' head(qd.cp)
##'
##' # This returns a matrix.  Note that you can also
##' # solve using name initial values. For example:
##'
##' inits <- c(eff = 1);
##'
##' qd.cp <- solve(m1, theta, events = qd, inits);
##'
##' @keywords models nonlinear
##' @concept Nonlinear regression
##' @concept ODE models
##' @concept Ordinary differential equations
##' @concept Pharmacokinetics (PK)
##' @concept Pharmacodynamics (PD)
##' @useDynLib RxODE, .registration=TRUE
##' @importFrom Rcpp evalCpp
##' @importFrom utils getFromNamespace assignInMyNamespace download.file head sessionInfo
##' @importFrom stats setNames update
##' @importFrom methods signature is
##' @export
RxODE <- function(model, modName = basename(wd), wd = ifelse(RxODE.cache.directory == ".", getwd(), RxODE.cache.directory),
                  filename = NULL, do.compile = NULL, extraC = NULL,
                  debug = FALSE,
                  calcJac=NULL, calcSens=NULL,
                  collapseModel=FALSE, ...) {
    if (missing(modName) && RxODE.delete.unnamed){
        do.delete <- TRUE;
    } else {
        do.delete <- FALSE;
    }
    if (!missing(model) && !missing(filename))
        stop("Must specify exactly one of 'model' or 'filename'.")
    if (missing(model) && !missing(filename)){
        model <- filename;
    }
    if (!missing(model) && missing(filename)){
        if (class(substitute(model)) == "{"){
            model <- deparse(substitute(model));
            if (model[1] == "{"){
                model <- model[-1];
                model <- model[-length(model)];
            }
            model <- paste(model, collapse="\n");
        } else if (class(model) == "RxODE"){
            model <- rxModelVars(model)$model["model"];
            if (!is.null(calcJac) && is.null(calcSens)){
                calcSens <- FALSE;
            }
        }
        ## else if ((class(model) == "function" || class(model) == "call")){
        ##     model <- deparse(body(model))[-1];
        ##     model <- paste(model[-length(model)], collapse="\n");
        ## }
    }
    model <- rxLinCmtTrans(model);
    ## RxODE compilation manager (location of parsed code, generated C,  shared libs, etc.)

    cmpMgr <- rx.initCmpMgr(model, modName, wd,  extraC, debug, missing(modName),
                            calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel);
    ## NB: the set of model variables (modelVars) is only available
    ## after parsing, thus it needs to be dynamically computed in cmpMgr
    get.modelVars <- cmpMgr$get.modelVars

    .version <- rxVersion()["version"]; # object version
    .last.solve.args <- NULL   # to be populated by solve()

    .c <- function(...){.C(...)};

    sexp.address <- NULL
    .sexp <- function(...){
        if (is.null(sexp.address)){
            sexp.address <<- getNativeSymbolInfo(as.vector(out$dll$modVars$trans["ode_solver_sexp"]),gsub("_$","",out$dll$modVars$trans["prefix"]))$address
        }
        .Primitive(".Call")(sexp.address,...);
    }

    solve <- function(params=NULL, events=NULL, inits = NULL, scale = c(),
                      covs = NULL, stiff = TRUE, transit_abs = NULL,
                      atol = 1.0e-8, rtol = 1.0e-6, maxsteps = 5000, hmin = 0, hmax = NULL, hini = 0, maxordn = 12,
                      maxords = 5, ..., covs_interpolation = c("linear", "constant"),
                      theta=numeric(), eta=numeric(), matrix=TRUE,add.cov=FALSE,
                      inC=FALSE, counts=NULL, do.solve=TRUE){
        env <- environment(.c);
        modVars <- dll$modVars;
        trans <- modVars$trans
        state <- modVars$state;
        lhs <- modVars$lhs;
        pars <- modVars$params;
        state.ignore <- modVars$state.ignore
        if (!is.null(params)){
            if (is.null(events) && class(params) == "EventTable"){
                events <- params;
                params <- c();
            }
        }
        if (is.null(transit_abs)){
            transit_abs <- modVars$podo;
            if (transit_abs){
                warning("Assumed transit compartment model since 'podo' is in the model.")
            }
        }
        if (class(params) != "numeric"){
            n <- names(params);
            params <- as.double(params);
            names(params) <- n;
        }
        ## Params and inits passed
        extra.args <- list(events = events$copy(),
                           covs = covs, stiff = stiff,
                           transit_abs = transit_abs, atol = atol, rtol = rtol, maxsteps = maxsteps,
                           hmin = hmin, hmax = hmax, hini = hini, maxordn = maxordn, maxords = maxords,
                           covs_interpolation = covs_interpolation, add.cov=add.cov, ...);
        params <- c(params, rxThetaEta(theta, eta));
        event.table <- events$get.EventTable()
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
            if (is.null(event.table$time) || length(event.table$time) == 1){
                hmax <- 0;
            } else {
                hmax <- max(abs(diff(event.table$time)))
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
        inits <- rxInits(dll, inits, state, 0);
        params <- rxInits(dll, params, pars, NA, !is.null(covs));
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
                print(params)
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
        lhs_vars <- lhs
        if (is.null(inits)){
            n <- state;
            inits <- rep(0.0, length(n));
            names(inits) <- n;
        }
        s <- as.list(match.call(expand.dots = TRUE))
        wh <- grep(pattern = "[Ss]\\d+$", names(s))
        if (length(scale) > 0 && length(wh) > 0){
            stop("Cannot specify both 'scale=c(...)' and S#=, please pick one to scale the ODE compartments.")
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
                warning(sprintf("Scaler variable(s) above the number of compartments: %s.",
                                paste(paste0("S", scaler.ix[scaler.ix > length(inits)]), collapse=", ")))
                scale <- scale[scaler.ix < length(inits)]
                scaler.ix <- scaler.ix[scaler.ix < length(inits)];
            }
            names(scale) <- state[scaler.ix];
        }
        scale <- c(scale);
        scale <- rxInits(dll, scale, state, 1, noini=TRUE);
        isLocf <- 0;
        if (length(covs_interpolation) > 1){
            isLocf <- 0;
        } else if (covs_interpolation == "constant"){
            isLocf <- 1;
        } else if (covs_interpolation != "linear"){
            stop("Unknown covariate interpolation specified.");
        }
        ## if (event.table$time[1] != 0){
        ##     warning(sprintf("The initial conditions are at t = %s instead of t = 0.", event.table$time[1]))
        ## }
        ## Ensure that inits and params have names.
        names(inits) <- state
        names(params) <- pars;

        time <- event.table$time;
        evid <- as.integer(event.table$evid);
        amt <- as.double(event.table$amt[event.table$evid>0]);
        ## Covariates
        pcov=as.integer(pcov);
        cov=as.double(cov);
        isLocf=as.integer(isLocf);
        ## Solver options (double)
        atol=as.double(atol);
        rtol=as.double(rtol);
        hmin=as.double(hmin);
        hmax=as.double(hmax);
        hini=as.double(hini);
        ## Solver options ()
        maxordn=as.integer(maxordn);
        maxords=as.integer(maxords);
        maxsteps=as.integer(maxsteps);
        stiff=as.integer(stiff);
        transit_abs=as.integer(transit_abs);
        do.matrix=as.integer(matrix);
        add.cov = as.integer(add.cov)
        if (do.solve){
            ret <- try({ret <- .sexp(## Parameters
                            params,
                            inits,
                            as.double(scale),
                            lhs_vars,
                            ## events
                            time,
                            evid,
                            amt,
                            ## Covariates
                            pcov,
                            cov,
                            isLocf,
                            ## Solver options (double)
                            atol,
                            rtol,
                            hmin,
                            hmax,
                            hini,
                            ## Solver options ()
                            maxordn,
                            maxords,
                            maxsteps,
                            stiff,
                            transit_abs,
                            ## Passed to build solver object.
                            env,
                            extra.args,
                            do.matrix,
                            add.cov,
                            state.ignore)
                rc <- ret[[2]];
                ret <- ret[[1]];
                ## attr(ret, "solveRxDll")$matrix <- attr(ret, "solveRxDll")$matrix[events$get.obs.rec(), ];
                ## Change sensitivities to be d/dt(d(A)/d(B)) form.
                ## dim <- dimnames(attr(ret, "solveRxDll")$matrix);
                ## dim[[2]] <- gsub(regSens,"d/dt(d(\\1)/d(\\2))",dim[[2]]);
                ## dimnames(attr(ret, "solveRxDll")$matrix) <- dim;
                if (rc != 0)
                    stop(sprintf("could not solve ODE, IDID = %d (see further messages)", rc))
                ret
            })
            if (inherits(ret, "try-error")){
                ## Error solving, try the other solver.
                ## errs <- paste(suppressWarnings({readLines(sink.file)}), collapse="\n");
                stiff <- 1L - stiff;
                ## sink(sink.file);
                try({ret <- .sexp(## Parameters
                         params,
                         inits,
                         as.double(scale),
                         lhs_vars,
                         ## events
                         time,
                         evid,
                         amt,
                         ## Covariates
                         pcov,
                         cov,
                         isLocf,
                         ## Solver options (double)
                         atol,
                         rtol,
                         hmin,
                         hmax,
                         hini,
                         ## Solver options ()
                         maxordn,
                         maxords,
                         maxsteps,
                         stiff,
                         transit_abs,
                         ## Passed to build solver object.
                         env,
                         extra.args,
                         do.matrix,
                         add.cov,
                         state.ignore)
                    rc <- ret[[2]];
                    ret <- ret[[1]];
                    ## attr(ret, "solveRxDll")$matrix <- attr(ret, "solveRxDll")$matrix[events$get.obs.rec(), ];
                    ## Change sensitivities to be d/dt(d(A)/d(B)) form.
                    ## dim <- dimnames(attr(ret, "solveRxDll")$matrix);
                    ## dim[[2]] <- gsub(regSens,"d/dt(d(\\1)/d(\\2))",dim[[2]]);
                    ## dimnames(attr(ret, "solveRxDll")$matrix) <- dim;

                    if (rc != 0)
                        stop(sprintf("Could not solve ODE, IDID = %d (see further messages).", rc))
                    ret
                })
                ## sink();
                if (inherits(ret, "try-error")){
                    stop("Tried both LSODA and DOP853, but could not solve the system.")
                } else {
                    if (stiff == 1L){
                        warning("Originally tried DOP853, but it failed to solve, so used LSODA instead.")
                    } else {
                        warning("Originally tried LSODA, but it failed to solve, so used DOP853 instead.")
                    }
                }
            }
            ## Now do scaling
        } else {
            ret <- list(params=params,
                        inits=inits,
                        lhs_vars=lhs_vars,
                        ## events
                        time=time,
                        evid=evid,
                        amt=amt,
                        ## Covariates
                        pcov=pcov,
                        cov=cov,
                        isLocf=isLocf,
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
                        transit_abs=transit_abs,
                        ## Passed to build solver object.
                        object=env,
                        extra.args=extra.args,
                        scale=scale,
                        events=events,
                        event.table=event.table,
                        do.matrix=do.matrix,
                        add.cov=add.cov,
                        state.ignore=state.ignore);

            if (inC){
                .Call("RxODE_ode_setup", inits, lhs_vars, time, evid, amt, pcov, cov, isLocf, atol, rtol, hmin, hmax,
                      hini, maxordn, maxords, maxsteps, stiff, transit_abs,
                      PACKAGE="RxODE");
            }
        }
        return(ret);
    }
    force <- FALSE
    if (class(do.compile) == "logical"){
        if (do.compile)
            force <- TRUE;
    }
    if (is.null(do.compile)){
        do.compile <- TRUE
    }
    ptr.address <- NULL
    assignPtr <- function(){
        cmpMgr$dynLoad();
        if (is.null(ptr.address)){
            ptr.address <<- getNativeSymbolInfo(as.vector(out$dll$modVars$trans["ode_solver_ptr"]),gsub("_$","",out$dll$modVars$trans["prefix"]))$address;
        }
        .Primitive(".Call")(ptr.address);
    }
    if (do.compile){
        cmpMgr$compile(force);
        ## Add for backward compatibility.
        tmp <- rxDll(cmpMgr$rxDll());
        names(tmp) <- NULL;
        cmpMgr$dllfile <- tmp;
        tmp <- rxTrans(cmpMgr$rxDll(), calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel)["ode_solver"]
        names(tmp) <- NULL
        cmpMgr$ode_solver <- tmp;
        model <- rxModelVars(cmpMgr$rxDll())$model["model"];
        names(model) <- NULL;
        cmpMgr$model <- model
    }
    out <-
       list(modName = modName,
            model = model,           # actual model code
            get.modelVars = get.modelVars,  # extract model variables (pars, lhs, etc)
            solve = solve,
            cmpMgr = cmpMgr,
            dynLoad = cmpMgr$dynLoad,
            dynUnload = cmpMgr$dynUnload,
            isValid = cmpMgr$isValid,
            version = .version,
            delete = cmpMgr$delete,
            ## the next is for backward compatibility and will be deprecated
            parse = cmpMgr$parse,
            compile = cmpMgr$compile,
            get.index = cmpMgr$get.index,
            run = solve,
            getObj = function(obj) get(obj, envir = environment(solve)),
            assignPtr=assignPtr)
    if (do.compile){
        vars <- rxModelVars(cmpMgr$rxDll());
        out$calcJac <- (length(vars$dfdy) > 0);
        out$calcSens <- (length(vars$sens) > 0)
        out$dll <- cmpMgr$rxDll();
        out$dll$old <- getNativeSymbolInfo(as.vector(out$dll$modVars$trans["ode_solver"]),gsub("_$","",out$dll$modVars$trans["prefix"]))$address
        dll <- out$dll;
    } else {
        out$dll <- NULL;
    }
    if (do.delete){
        reg.finalizer(environment(solve), rxodeGc, onexit=TRUE);
    }
    class(out) <- "RxODE"
   out
}
##' Get model properties without compiling it.
##'
##' @param model RxODE specification
##' @inheritParams RxODE
##' @return RxODE trans list
##' @author Matthew L. Fidler
##' @keywords internal
rxGetModel <- function(model, calcSens=FALSE, calcJac=FALSE, collapseModel=FALSE){
    if (class(substitute(model)) == "call"){
        model <- model;
    }
    if (class(substitute(model)) == "{"){
        model <- deparse(substitute(model))
        if (model[1] == "{"){
            model <- model[-1];
            model <- model[-length(model)];
        }
        model <- paste(model, collapse="\n");
    } else if (class(model) == "function" || class(model) == "call"){
        model <- deparse(body(model));
        if (model[1] == "{"){
            model <- model[-1];
            model <- model[-length(model)];
        }
        model <- paste(model, collapse="\n");
    } else if (class(model) == "character"){
        if (file.exists(model)){
            ret$use_model_name <- TRUE;
        }
    } else if (class(model) == "name"){
        model <- eval(model);
    } else {
        stop(sprintf("Can't figure out how to handle the model argument (%s).", class(model)));
    }
    parseModel <- tempfile();
    cFile <- tempfile();
    on.exit({unlink(parseModel); unlink(cFile)});
    sink(parseModel);
    cat(model);
    cat("\n");
    sink();
    return(rxTrans(parseModel, cFile, calcSens=calcSens, calcJac=calcJac, collapseModel=collapseModel, modVars=TRUE));
}

rxGetModel.slow <- NULL

rxAdd <- function(rx, pre, post, ...){
    base <- rxNorm(rx);
    if (!missing(pre)){
        pre <- rxNorm(rxGetModel(pre))
    } else {
        pre <- NULL;
    }
    if (!missing(post)){
        post <- rxNorm(rxGetModel(post))
    } else {
        post <- NULL;
    }
    return(RxODE(paste(c(pre, base, post), collapse="\n"),...));
}

##' Predict an RxODE object
##'
##' \code{predict} solves the odinary differential equations specified by
##'  RxODE object
##'
##' @param object An RxODE object
##'
##' @param ... Solve arguments sent to \code{rxSolve}.  See
##'     \code{\link{rxSolve}}.
##'
##' @author Matthew L.Fidler
##' @export
predict.RxODE <- function(object, ...){
    rxSolve(object, ...);
}

##' Solve RxODE objects
##'
##' @param ... Additional arguments sent to \code{rxSolve}
##'
##' @seealso \code{\link{rxSolve}}
##'
##' @author Matthew L.Fidler
##' @export
solve.RxODE <- function(...){
    rxSolve(...)
}
##' @rdname solve.RxODE
##' @export
solve.RxCompilationManager <- function(...){
    rxSolve(...)
}
##' @rdname solve.RxODE
##' @export
solve.solveRxDll <- function(...){
    rxSolve(...)
}
##' @rdname solve.RxODE
##' @export
solve.character <- function(...){
    rxSolve(...)
}

##' @rdname solve.RxODE
##' @export
solve.rxDll <- function(...){
    rxSolve(...);
}

##' Add item to solved system of equations
##'
##' @title rxChain  Chain or add item to solved system of equations
##'
##' @param obj1 Solved object.
##'
##' @param obj2 New object to be added/piped/chained to solved object.
##'
##' @return When \code{newObject} is an event table, return a new
##'     solved object with the new event table.
##'
##' @author Matthew L. Fidler
##'
##' @export
rxChain <- function(obj1, obj2) {
    args <- rev(as.list(match.call())[-1]);
    names(args) <- c("obj", "solvedObject");
    return(do.call("rxChain2", args, envir = parent.frame(1)));
}

##' @rdname rxChain
##' @export
'+.solveRxDll' <- function(obj1, obj2){
    return(rxChain(obj1, obj2))
}

##' Second command in chaining commands
##'
##' This is s3 method is called internally with \code{+} and \code{\%>\%} operators.
##'
##' @param obj the object being added/chained/piped to the solved object
##' @param solvedObject the solved object
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
rxChain2 <- function(obj, solvedObject){
    UseMethod("rxChain2")
}

##' @rdname rxChain2
##' @export
rxChain2.default <- function(obj, solvedObject){
    args <- as.list(match.call());
    stop(sprintf("Do not know how to add %s to RxODE solved object %s.", toString(args[[2]]), toString(args[[3]])))
}

##' @rdname rxChain2
##' @export
rxChain2.EventTable <- function(obj, solvedObject){
    args <- rev(as.list(match.call())[-1]);
    names(args) <- c("object", "events");
    return(do.call("rxSolve", args, envir = parent.frame(1)));
}

##' Print information about the RxODE object.
##'
##' This prints the model name and its status for being able to be solved
##'
##' @param x An rxode object
##' @param ... Ignored parameters
##' @author Matthew L.Fidler
##' @export
print.RxODE <-
    function(x, ...)
{
    valid <- x$cmpMgr$isValid()
    if (!valid){
        .msg <- "invalid object, needs to be re-created"
    } else {
        .ready <- x$cmpMgr$getObj(".compiled")
        .msg <- if (.ready) "ready to run" else "needs compilation"
    }
    dll <- basename(rxDll(x));
    dll <- substr(dll, 1, nchar(dll) - nchar(.Platform$dynlib.ext) - nchar(.Platform$r_arch) - 1)
    cat(sprintf('RxODE model named "%s" (%s).\n', dll, .msg))
    if (!any(names(list(...)) == "rxSuppress") && valid){
        cat(sprintf('States: %s\n', paste(rxState(x), collapse=", ")))
        cat(sprintf('Params: %s\n', paste(rxParams(x), collapse=", ")))
    }
    invisible(x)
}

##' Print expanded information about the RxODE object.
##'
##' This prints the expanded information about the RxODE object.
##'
##' @param object RxODE object
##' @param ... Ignored parameters
##' @author Matthew L.Fidler
##' @export
summary.RxODE <- function(object, ...)
{
    print.RxODE(object, rxSuppress=TRUE);
    summary.rxDll(object$cmpMgr$rxDll(), noprint = TRUE)
    invisible(object)
}

##' @rdname summary.RxODE
##' @export
summary.RxCompilationManager <- function(object, ...)
{
    print.RxCompilationManager(object);
    summary.rxDll(object$rxDll(), noprint = TRUE)
    invisible(object);
}



##' Return the RxODE coefficients
##'
##' This returns the parameters , state variables
##'
##' @param object is an RxODE object
##' @param ... ignored arguments
##'
##' @return a rxCoef object with the following
##'
##' \item{params}{ is a list of strings for parameters for the RxODE object}
##' \item{state}{ is a list of strings for the names of each state in
##'     the RxODE object.}
##' \item{ini}{ is the model specified default values for the
##'     parameters.}
##' \item{RxODE}{ is the referring RxODE object}
##' @author Matthew L.Fidler
##' @importFrom stats coef
##' @export
coef.RxODE <- function(object,
                       ...){
    ret <- rxModelVars(object)[c("params", "state", "ini", "sens", "fn.ini")];
    ret$RxODE <- object;
    class(ret) <- "rxCoef";
    return(ret);
}

##' @rdname coef.RxODE
##' @export
coef.RxCompilationManager <- function(...){
    coef.RxODE(...);
}

##' @rdname coef.RxODE
##' @export
coef.solveRxODE <- function(object, ...){
    env <- attr(object, ".env");
    rxode <- env$env$out;
    lst <- c(list(params=env$params, inits=env$inits), env$extra.args)
    ret <- list();
    ret$params <- lst$params;
    ret$state <- lst$inits;
    ret$sens <- rxModelVars(rxode)["sens"];
    ret$state <- ret$state[regexpr(regSens, names(ret$state)) == -1]
    ret$RxODE <- rxode;
    class(ret) <- "rxCoefSolve";
    return(ret);
}

##' @rdname coef.RxODE
##' @export
coef.rxDll <- function(...){
    coef.RxODE(...);
};

##' Print the rxCoef object
##'
##' This prints out the user supplied arguments for rxCoef object
##'
##' @param x rxCoef object
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
print.rxCoef <- function(x, ...){
    rxDllObj <- x$RxODE;
    if (length(rxParams(rxDllObj)) > 0){
        cat("\nUser supplied parameters:\n");
        print(rxInits(rxDllObj, c(), rxParams(rxDllObj), NA, TRUE))
        cat("\nUser initial conditions:\n");
        tmp <- rxInits(rxDllObj, c(), rxState(rxDllObj), 0, TRUE);
        if (length(x$sens) > 0){
            tmp <- tmp[regexpr(regSens, names(tmp)) == -1];
        }
        tmp <- tmp[!(names(tmp) %in% x$fn.ini)];
        print(tmp);
    }
    if (length(x$fn.ini) > 0){
        cat("\nParameter-based initial conditions:\n");
        print(x$fn.ini);
    }
    cat("\nCompartments:\n");
    tmp <- rxState(rxDllObj);
    if (length(tmp) > 0){
        names(tmp) <- paste0("cmt=", 1:length(tmp));
        if (length(x$sens) > 0){
            tmp1 <- tmp[regexpr(regSens, tmp) == -1];
            print(tmp1);
            cat("\nSensitivities:\n");
            tmp2 <- gsub(regSens, "d/dt(d(\\1)/d(\\2))", tmp[regexpr(regSens, tmp) != -1]);
            print(tmp2);
        } else {
            print(tmp);
        }
    } else {
        cat("No ODEs in this DLL.\n");
    }
    return(invisible());
}

##' Print the rxCoefSolve object
##'
##' This prints out the user supplied arguments for the rxCoef object
##'
##' @param x rxCoefSolve object
##' @param ... Other (ignored) parameters.
##' @keywords Internal
##' @author Matthew L.Fidler
##' @export
print.rxCoefSolve <- function(x, ...){
    cat("\nUser supplied parameters ($params):\n");
    print(x$params);
    cat("\nUser initial conditions ($state):\n");
    print(x$state);
    rxDllObj <- x$RxODE;
    if (length(rxInits(rxDllObj)) > 0){
        cat("\nDefault parameter values:\n")
        print(rxInits(rxDllObj));
    }
    cat("\nCompartments:\n");
    tmp <- rxState(rxDllObj);
    names(tmp) <- paste0("cmt=", 1:length(tmp));
    sens <- tmp[regexpr(regSens, tmp) != -1];
    if (length(sens) > 0){
        tmp <- tmp[regexpr(regSens, tmp) == -1]
    }
    print(tmp);
    if (length(sens) > 0){
        sens <- gsub(regSens, "d/dt(d(\\1)/d(\\2))", sens);
        cat("\nSensitivities:\n");
        print(sens)
    }
    return(invisible());
}

##' A compilation manager for RxODE models
##'
##' This function parses, compiles, links, and loads the shared object
##' that implements an RxODE model.
##'
##' The function parses and compiles (if needed) the \code{RxODE}
##' model specified in the string \code{model} into a dynamic link
##' library (DLL on Windows) or a shared object (\code{*.so} on
##' Unix-like systems).
##'
##' It then dynamically loads this code into the current R
##' session. (Models previously parsed and compiled in previous R
##' sessions only need to be dynamically loaded into the current R session.)
##'
##' @inheritParams RxODE
##'
##' @param mmod A boolean telling if the modName from \code{RxODE} was
##'     missing.  This affects how the model is created and used.
##' @return  An object (closure) with the following member functions:
##' \item{parse}{
##'     this function parses (translates) the ODE-based model
##'     specification and generates a C file to implement the model.}
##' \item{compile}{
##'     compiles the generated C file for the ODE system
##'     and dynamically loads the machine code in the shared object.}
##' \item{dynLoad}{
##'     if needed, dynamically loads the dynamic library
##'     produced in the \code{compile()} step.  Note that the shared
##'     object persists across R sessions, thus the \code{dynLoad} needs
##'     to be issued as needed.}
##' \item{dynUnload}{
##'     this function unloads the previously dynamically loaded
##'     model object code.  Mainly for internal use.}
##' \item{ode_solver}{
##'     a string with the name of the C symbol for this model solver.}
##' \item{dllfile}{
##'     a string with the name of the dynamic link (or shared object) file.}
##' \item{get.modelVars}{
##'     function that returns a list with 3 character
##'     vectors, \code{params}, \code{state}, and \code{lhs} of variable
##'     names (identifiers) used in the model specification.
##'     These will be output when the model is computed (i.e., the ODE solved).}
##' \item{isValid}{
##'     a function that (naively) checks for model validity,
##'     namely that the C object code reflects the latest model
##'     specification.}
##' \item{get.index}{
##'     helper function to extract the index of one or
##'     more system variables (state, parameter, or other).}
##' \item{getObj}{
##'     internal (not user callable) function.}
##' @examples
##' \dontrun{
##'   cmpMgt <- rx.initCmpMgr(model, "tst1", wd = ".")
##' }
##' @keywords internal models ODE
##' @concept ordinary differential equations
##' @seealso \code{\link{RxODE}}
##' @author Matthew L.Fidler
##' @export rx.initCmpMgr
rx.initCmpMgr <-
    function(model, modName, wd, extraC = NULL, debug = TRUE, mmod = FALSE, calcJac=NULL, calcSens=NULL, collapseModel=FALSE)
{
    ## Initialize the RxODE compilation manager.  This is a stub
    ## function for backward compatability.
    .model <- model;
    .mmod <- mmod;
    .modName <- modName;
    .wd <- wd;
    .parsed <- FALSE;
    .compiled <- FALSE;
    .calcJac <- calcJac;
    .calcSens <- calcSens;
    .collapseModel <- collapseModel
    ## model-specific directory under .md (default current dir)
    if (mmod){
        .mdir <- .wd;
    } else {
        .mdir <- file.path(.wd, sprintf("%s.d", .modName));
    }
    .rxDll <- NULL;
    .debug <- debug;
    .extraC <- extraC;

    ## filenames and shell command required for parsing (these are unique
    ## within a given model directory

    ## files needed for compiling the C output of the parsed ODE
    parse <- function(force = FALSE){
        do.it <- force || !.parsed;
        if (!do.it)
            return(invisible(.parsed));
        .parsed <<- TRUE;
        .compiled <<- FALSE;
        invisible(.parsed);
    }

    compile <- function(force = FALSE){
        do.it <- force|| !.compiled;
        if (!do.it)
            return(invisible(.compiled));
        lwd <- getwd();
        if (!file.exists(.wd))
            dir.create(.wd, recursive = TRUE)
        if (!file.exists(.wd))
            setwd(.wd);
        on.exit(setwd(lwd));
        if (.mmod){
            .rxDll <<- rxCompile(.model, extraC = .extraC, debug = .debug, calcJac=.calcJac, calcSens=.calcSens, collapseModel=.collapseModel);
        } else {
            .rxDll <<- rxCompile(.model, .mdir, extraC = .extraC, debug = .debug, modName = .modName,  calcJac=.calcJac, calcSens=.calcSens, collapseModel=.collapseModel);
        }
        if (class(.rxDll) == "rxDll"){
            .compiled <<- TRUE;
        }
        invisible(.compiled);
    }

    dynLoad <- function(force = FALSE){
        ## NB: we may need to reload (e.g., when we re-start R and
        ## re-instantiate the RxODE object from a save.image.
        if (!rxDllLoaded(.rxDll) && !file.exists(rxDll(.rxDll)) && RxODE.compile.on.load){
            compile();
        }
        rxLoad(.rxDll);
        return();
    }

    dynUnload <- function(){
        return(rxUnload(.rxDll));
    }
    delete <- function(){
        .parsed <<- FALSE;
        .compiled <<- FALSE;
        rxDelete(.rxDll);
        ## TODO: should we remove all objects in the closure?
        ## as the object is no longer valid. Need a valid.object()
        ##unlink(.mdir, recursive = TRUE) # leave dir
    }

    isValid <- function(){
        return(file.exists(rxDll(.rxDll)));
    }

    get.index <- function(s) {
        ## return the (one) state varible index
        if (.compiled){
            return(rxState(.rxDll, s));
        } else {
            stop("Needs to be compiled first.");
        }
    }

   out <-
       list(parse         = parse,
            compile       = compile,
            model         = model,
            dynLoad       = dynLoad,
            dynUnload     = dynUnload,
            ode_solver    = "Need to compile",   # name of C function
            modelDir      = .mdir, # model directory
            dllfile       = "Need to compile",
            get.modelVars = function(){
                mv <- rxModelVars(.rxDll);
                ret <- rxModelVars(.rxDll)[c("params", "state", "lhs")];
                init <- rxInit(.rxDll);
                ret$params <- ret$params[!(ret$params %in% names(init))]
                return(ret);
            },
            isValid       = isValid,
            delete        = delete,
            get.index     = get.index,
            getObj        = function(obj) get(obj, envir = environment(parse)),
            extraC        = .extraC,
            rxDll         = function() .rxDll
            )
    class(out) <- "RxCompilationManager"
   out
}

##' @rdname print.RxODE
##' @export
"print.RxCompilationManager" <-
function(x, ...)
{
    modName <- x$getObj(".modName")
    cat(sprintf("RxCompilationManager for RxODE model '%s'\n", modName))
   invisible(x)
}

rxPrefix <- function(model,          # Model or file name of model
                     modName = NULL, # Model name, overrides calculated model name.
                     calcJac=NULL,
                     calcSens=NULL,
                     collapseModel=FALSE,
                     ...){
    ## rxPrefix returns a prefix for a model.
    if (!is.null(modName)){
        modelPrefix <- sprintf("%s_", gsub("\\W", "_", modName));
    } else if (file.exists(model)){
        modelPrefix <- sprintf("%s_", gsub("\\W", "_", gsub("[.].*$", "", base::basename(model))));
    } else {
        parseModel <- tempfile();
        cFile <- tempfile();
        on.exit({unlink(parseModel); unlink(cFile)});
        sink(parseModel);
        cat(model);
        cat("\n");
        sink();
        trans <- rxTrans(parseModel, cFile, calcSens=calcSens, calcJac=calcJac, collapseModel=collapseModel)
        modelPrefix <- sprintf("rx_%s_", trans["parsed_md5"]);
    }
    modelPrefix <- sprintf("%s%s_", modelPrefix, .Platform$r_arch);
    return(modelPrefix);
} # end function rxPrefix



##' Return the md5 of an RxODE object or file
##'
##' This md5 is based on the model and possibly the extra c code
##' supplied for the model.  In addition the md5 is based on syntax
##' options, compiled RxODE library md5, and the RxODE
##' version/repository.
##'
##' @inheritParams RxODE
##'
##'
##' @param ... ignored arguments
##'
##' @return If this is a RxODE object, return a named list:
##'
##' \item{\code{file_md5}}{ is the model's file's md5}
##'
##' \item{\code{parsed_md5}}{ is the parsed model's file's md5.}
##'
##' Otherwise return the md5 based on the arguments provided
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export rxMd5
rxMd5 <- function(model,         # Model File
                  extraC  = NULL, # Extra C
                  calcJac =NULL,
                  calcSens=NULL,
                  collapseModel=FALSE,
                  ...){
    ## rxMd5 returns MD5 of model file.
    ## digest(file = TRUE) includes file times, so it doesn't work for this needs.
    if (class(model) == "character"){
        if (length(model) == 1){
            if (file.exists(model)){
                ret <- suppressWarnings({readLines(model)});
                mod <- paste(ret, collapse = "\n");
            } else {
                ret <- model
                mod <- model;
            }
        } else {
            ret <- model;
            mod <- paste(ret, collapse="\n");
        }
        if (class(extraC) == "character"){
            if (file.exists(extraC)){
                ret <- c(ret, gsub(rex::rex(or(any_spaces, any_newlines)), "", readLines(extraC), perl = TRUE));
            }
        }
        tmp <- c(RxODE.syntax.assign, RxODE.syntax.star.pow, RxODE.syntax.require.semicolon, RxODE.syntax.allow.dots, RxODE.syntax.allow.ini0, RxODE.syntax.allow.ini, RxODE.calculate.jacobian, RxODE.calculate.sensitivity, calcJac, calcSens, collapseModel);
        ret <- c(ret, tmp);
        tmp <- getLoadedDLLs()$RxODE;
        class(tmp) <- "list";
        ## new RxODE DLLs gives different digests.
        ret <- c(ret, digest::digest(tmp$path,file=TRUE, algo="md5"));
        ## Add version and github repository information
        ret <- c(ret, rxVersion());
        return(list(text = mod,
                    digest = digest::digest(ret, algo="md5")));
    } else {
        rxModelVars(model)$md5;
    }
} # end function rxMd5
##' Translate the model to C code if needed
##'
##' This function translates the model to C code, if needed
##'
##'
##' @inheritParams RxODE
##'
##'
##' @param cFile The C file where the code should be output
##'
##' @param modelPrefix Prefix of the model functions that will be
##'     compiled to make sure that multiple RxODE objects can coexist
##'     in the same R session.
##'
##' @param md5 Is the md5 of the model before parsing, and is used to
##'     embed the md5 into DLL, and then provide for functions like
##'     \code{\link{rxModelVars}}.
##'
##' @param ... Ignored parameters.
##'
##'
##' @param modVars returns the model variables instead of the named
##'     vector of translated properties.
##'
##'
##'
##' @return a named vector of translated model properties
##'       including what type of jacobian is specified, the \code{C} function prefixes,
##'       as well as the \code{C} functions names to be called through the compiled model.
##' @seealso \code{\link{RxODE}}, \code{\link{rxCompile}}.
##' @author Matthew L.Fidler
##' @export
rxTrans <- function(model,
                    cFile       = sprintf("%s.c", gsub("[.][^.]*$", "", model)), # C file output
                    extraC      = NULL,                                       # Extra C file(s)
                    modelPrefix = "",                                         # Model Prefix
                    md5         = "",                                         # Md5 of model
                    modName     = NULL,                                       # Model name for DLL
                    modVars     = FALSE,                                      # Return modVars
                    calcSens=NULL,
                    calcJac=NULL,
                    collapseModel=FALSE,
                    ...){
    UseMethod("rxTrans");
} # end function rxTrans


##' @rdname rxTrans
##' @export
rxTrans.default <- function(model,
                            cFile       = sprintf("%s.c", gsub("[.][^.]*$", "", model)), # C file output
                            extraC      = NULL,                                       # Extra C file(s)
                            modelPrefix = "",                                         # Model Prefix
                            md5         = "",                                         # Md5 of model
                            modName     = NULL,                                       # Model name for DLL
                            modVars     = FALSE,                                      # Return modVars
                            calcSens=NULL,
                            calcJac=NULL,
                            collapseModel=FALSE,
                            ...){
    mv <- rxModelVars(model)
    if (modVars){
        return(mv);
    } else {
        return(c(mv$trans, mv$md5));
    }
}

##' @rdname rxTrans
##' @export
rxTrans.character <- function(model,
                              cFile       = sprintf("%s.c", gsub("[.][^.]*$", "", model)), # C file output
                              extraC      = NULL,                                       # Extra C file(s)
                              modelPrefix = "",                                         # Model Prefix
                              md5         = "",                                         # Md5 of model
                              modName     = NULL,                                       # Model name for DLL
                              modVars     = FALSE,                                      # Return modVars
                              calcSens=NULL,
                              calcJac=NULL,
                              collapseModel=FALSE,
                              ...){
    ## rxTrans returns a list of compiled properties
    if (is.null(calcSens)){
        calcSens <- RxODE.calculate.sensitivity
    }
    if (is.null(calcJac)){
        calcJac <- RxODE.calculate.jacobian;
    }
    if (missing(modelPrefix)){
        modelPrefix <- rxPrefix(model, modName, calcSens, calcJac, collapseModel);
    }
    if (file.exists(model)){
        if (missing(md5)){
            md5 <- rxMd5(model, extraC, calcJac, calcSens, collapseModel)$digest;
        }
    } else {
        stop("This only translates a file (currently; Try rxCompile).");
    }
    parseModel <- tempfile();
    out3 <- tempfile();
    on.exit(unlink(parseModel));
    rxReq("dparser");
    ret <- .Call(trans, model, model, cFile, extraC, modelPrefix, md5, parseModel, out3, PACKAGE="RxODE");
    md5 <- c(file_md5 = md5, parsed_md5 = rxMd5(c(ret$model["normModel"],
                                                  ret$ini,
                                                  ret$state,
                                                  ret$params,
                                                  ret$lhs), extraC, calcJac, calcSens, collapseModel)$digest);
    ret$md5 <- md5
    ## dparser::dpReload();
    ## rxReload()
    if (file.exists(cFile)){
        ret$md5 <- md5
        if (class(calcSens) == "logical"){
            if (!calcSens){
                calcSens <- NULL;
            }
        }
        if (!is.null(calcSens)){
            if (length(rxState(ret)) <= 0){
                stop("Sensitivities do not make sense for models without ODEs.")
            }
            new <- rxSymPySensitivity(rxModelVars(rxNorm(ret)), calcSens=calcSens, calcJac=calcJac,
                                      collapseModel=collapseModel);
            expandModel <- tempfile();
            sink(expandModel);
            cat(new);
            cat("\n");
            sink()
            ret <- try({.Call(trans, model, expandModel, cFile, extraC, modelPrefix, md5, parseModel, out3)}, silent=TRUE);
            md5 <- c(file_md5 = md5, parsed_md5 = rxMd5(ret$model["normModel"], extraC, calcJac, calcSens, collapseModel)$digest)
            ret$md5 <- md5
            if (inherits(ret, "try-error")){
                rxCat("\n", paste(readLines(expandModel), collapse="\n"), "\n");
                stop("Syntax Errors")
            }
            ## dparser::dpReload();
            ## rxReload();
            unlink(expandModel);
            ret$md5 <- md5;
        } else if (calcJac){
            if (length(rxState(ret)) <= 0){
                stop("Jacobians do not make sense for models without ODEs.")
            }
            new <- rxSymPyJacobian(rxModelVars(rxNorm(ret)));
            expandModel <- tempfile();
            sink(expandModel);
            cat(new);
            cat("\n");
            sink()
            ret <- .Call("trans", model, expandModel, cFile, extraC, modelPrefix, md5, parseModel, out3, PACKAGE="RxODE");
            md5 <- c(file_md5 = md5, parsed_md5 = rxMd5(ret$model["normModel"], extraC, calcJac, calcSens, collapseModel)$digest)
            ret$md5 <- md5
            ## dparser::dpReload();
            ## rxReload();
            unlink(expandModel);
            ret$md5 <- md5;
        }
        if (modVars){
            return(ret)
        } else {
            return(c(ret$trans, ret$md5));
        }
    } else {
        rxCat(model);
        stop("Syntax Error (see above)");
    }
}
rxTransMakevars <- function(rxProps,                                                                              # rxTrans translation properties
                            rxDll, # Dll of file
                            cFile,
                            compileFlags =c("parsed_md5", "ode_solver", "ode_solver_sexp", "ode_solver_0_6",
                                            "ode_solver_ptr", "inis",
                                            "model_vars", "calc_lhs", "calc_jac", "dydt"), # List of compile flags
                            debug        = FALSE,                                                                 # Debug compile?
                            ...){
    ## rxTransCompileFlags returns a string for the compiler options
    ctxt <- suppressWarnings({readLines(cFile)});
    neededProps <- c("jac", compileFlags);
    if (all(neededProps %in% names(rxProps))){
        ret <- "";
        if (rxProps["jac"] == "fulluser"){
            ctxt <- gsub("__JT__", "1", ctxt)
            ctxt <- gsub("__MF__", "21", ctxt)
        } else if (rxProps["jac"] == "fullint"){
            ctxt <- gsub("__JT__", "2", ctxt)
            ctxt <- gsub("__MF__", "22", ctxt)
        }
        tmp <- rxProps[compileFlags];
        for (x in c("parsed_md5", "ode_solver", "ode_solver_sexp", "ode_solver_0_6", "ode_solver_ptr", "ode_solver_focei_outer",
                    "model_vars")){
            tmp[sprintf("%s_str", x)] <- sprintf("\"%s\"", tmp[x]);
        }
        tmp["lib_str"] <- sprintf("\"%s\"", gsub(.Platform$dynlib.ext, "", basename(rxDll)));
        for (flg in names(tmp)){
            ctxt <- gsub(sprintf("__%s__", toupper(flg)), tmp[flg], ctxt);
        }
        ctxt <- gsub("__R_INIT__", sprintf("R_init_%s", gsub(.Platform$dynlib.ext, "", basename(rxDll))), ctxt)
        writeLines(ctxt, cFile);
        ret <- ""
        if (debug){
            ret <- sprintf("%s -D__DEBUG__", ret);
        }
        ret <- sprintf("PKG_CPPFLAGS=%s\nPKG_LIBS=$(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)", ret);
        cat(ret);
        return(ret);
    } else {
        cat("Needed Variables: %s\n", paste(neededProps, collapse=","));
        stop(sprintf("Cannot compile, only found %s.", paste(names(rxProps), collapse=",")));
    }
} # end function rxTransCompileFlags

##' Determine if the rxDll is loaded or not.
##'
##' @param x is a RxODE family of objects
##'
##' @param retry is a flag to retry to load if the function can't
##'     determine if the object is loaded or not...
##'
##' @return a boolean stating if the DLL is loaded
##'
##' @author Matthew L.Fidler
##' @export
rxDllLoaded <- function(x, retry = TRUE){
    if (is.null(x)){
        return(FALSE);
    }
    m <- rxModelVars(x)$trans;
    if (any(names(m) == "ode_solver")){
        return(is.loaded(m["ode_solver"]));
    } else if (retry) {
        m <- rxCompile(x, force = FALSE);
        return(rxDllLoaded(m, retry = FALSE))
    } else {
        print(m);
        options(RxODE.echo.compile = TRUE);
        rxSyncOptions()
        on.exit({options(RxODE.echo.compile = TRUE);
            rxSyncOptions();})
        m <- rxCompile(x, force = FALSE);
        stop(sprintf("Can't figure out if the object is loaded (%s)...", .Platform$dynlib.ext));
    }
}
##' Compile a model if needed
##'
##' This is the compilation workhorse creating the RxODE model DLL
##' files.
##'
##' @inheritParams RxODE
##'
##' @param dir This is the model directory where the C file will be
##'     stored for compiling.
##'
##'     If unspecified, the C code is stored in a temporary directory,
##'     then the model is compiled and moved to the current directory.
##'     Afterwards the C code is removed.
##'
##'     If specified, the C code is stored in the specified directory
##'     and then compiled in that directory.  The C code is not removed
##'     after the DLL is created in the same directory.  This can be
##'     useful to debug the c-code outputs.
##'
##' @param prefix is a string indicating the prefix to use in the C
##'     based functions.  If missing, it is calculated based on file
##'     name, or md5 of parsed model.
##'
##'
##' @param force is a boolean stating if the (re)compile should be
##'     forced if RxODE detects that the models are the same as already
##'     generated.
##'
##' @param ... Other arguments sent to the \code{\link{rxTrans}}
##'     function.
##'
##' @return An rxDll object that has the following components
##'
##' \item{dll}{DLL path}
##' \item{model}{model specification}
##' \item{.c}{A function to call C code in the correct context from the DLL
##'          using the \code{\link{.C}} function.}
##' \item{.call}{A function to call C code in the correct context from the DLL
##'          using the \code{\link{.Call}} function.}
##' \item{args}{A list of the arguments used to create the rxDll object.}
##' @seealso \code{\link{RxODE}}
##' @author Matthew L.Fidler
##' @export
rxCompile <- function(model, dir, prefix, extraC = NULL, force = FALSE, modName = NULL,
                      calcJac=NULL, calcSens=NULL, collapseModel=FALSE,
                      ...){
    UseMethod("rxCompile")
}
##' @rdname rxCompile
##' @export
rxCompile.character <-  function(model,           # Model
                                 dir,             # Directory
                                 prefix=NULL,     # Prefix
                                 extraC  = NULL,  # Extra C File.
                                 force   = FALSE, # Force compile
                                 modName = NULL,  # Model Name
                                 calcJac=NULL, # Calculate Jacobian
                                 calcSens=NULL, # Calculate Sensitivity
                                 collapseModel=FALSE,
                                 ...){
    ## rxCompile returns the DLL name that was created.
    dllCopy <- FALSE;
    if (missing(dir)){
        dir <- tempfile();
        dllCopy <-  TRUE;
        on.exit(unlink(dir, recursive = TRUE))
    }
    if (missing(prefix)){
        prefix <- rxPrefix(model, modName, calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel);
    }
    if (!file.exists(dir))
        dir.create(dir, recursive = TRUE)
    cFile <- file.path(dir, sprintf("%s.c", substr(prefix, 0, nchar(prefix)-1)));
    cDllFile <- file.path(dir, sprintf("%s%s", substr(prefix, 0, nchar(prefix)-1), .Platform$dynlib.ext));
    if (dllCopy){
        finalDll <- file.path(getwd(), basename(cDllFile));
        finalC <- file.path(getwd(), basename(cFile));
    } else {
        finalDll <-  cDllFile;
        finalC <-  cFile;
    }
    if (!file.exists(model)){
        mFile <- sprintf("%s.rx", substr(cFile, 0, nchar(cFile)-2));
        sink(mFile);
        cat(model);
        cat("\n");
        sink();
    } else {
        mFile <- model;
    }
    md5 <- rxMd5(mFile, extraC, calcJac, calcSens, collapseModel);
    allModVars <- NULL;
    needCompile <- TRUE
    if (file.exists(finalDll)){
        try(dyn.load(finalDll, local = FALSE), silent = TRUE);
        modVars <- sprintf("%smodel_vars", prefix);
        if (is.loaded(modVars)){
            allModVars <- eval(parse(text = sprintf(".Call(\"%s\")", modVars)), envir = .GlobalEnv)
            modVars <- allModVars$md5;
        }
        if (!any(names(modVars) == "file_md5")){
            needCompile <- FALSE;
        } else {
            if (modVars["file_md5"] == md5$digest){
                needCompile <- FALSE;
            }
        }
    }
    if (force || needCompile){
        Makevars <- file.path(dir, "Makevars");
        trans <- rxTrans(mFile, cFile = cFile, md5 = md5$digest, extraC = extraC, ..., modelPrefix = prefix, calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel);
        if (file.exists(finalDll)){
            if (modVars["parsed_md5"] == trans["parsed_md5"]){
                rxCat("Don't need to recompile, minimal change to model detected.\n");
                needCompile <- FALSE;
            }
        }
        if (force || needCompile){
            ## Setup Makevars
            owd <- getwd();
            on.exit({if (file.exists(Makevars)){
                         unlink(Makevars);
                     };
                         setwd(owd);
            });
            if (file.exists(Makevars)){
                unlink(Makevars);
            }
            ## Now create C file
            rxTrans(mFile, cFile = cFile, md5 = md5$digest, extraC = extraC, ..., modelPrefix = prefix, calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel)
            sink(Makevars);
            cat(rxTransMakevars(trans, finalDll, cFile, ...));
            sink();
            if (dllCopy){
                file.copy(cFile, finalC);
            }
            sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths
            ## Change working directory
            setwd(dir);
            try(dyn.unload(finalDll), silent = TRUE);
            try(unlink(finalDll));
            cmd <- sprintf("%s/bin/R CMD SHLIB %s",
                           Sys.getenv("R_HOME"), base::basename(cFile));
            if (RxODE.echo.compile){
                cat(sprintf("%s\n", cmd));
            }
            compileFile <- tempfile();
            stdErrFile <- tempfile();
            rc <- tryCatch(rx.do.call(sh, list(cmd, ignore.stdout = !RxODE.echo.compile, ignore.stderr = !RxODE.echo.compile)),
                           error = function(e) "error",
                           warning = function(w) "warning");
            if (any(rc == c("error", "warning"))){
                try(rx.do.call(sh, list(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)),
                    silent = FALSE)
                rxCat("\n\nModel:\n", paste(readLines(mFile), collapse="\n"), "\n")
                stop(sprintf("error compiling %s", cFile));
            }
            if (dllCopy){
                file.copy(cDllFile, finalDll);
            }
            try(dyn.load(finalDll, local = FALSE), silent = TRUE);
            modVars <- sprintf("%smodel_vars", prefix);
            if (is.loaded(modVars)){
                allModVars <- eval(parse(text = sprintf(".Call(\"%s\")", modVars)), envir = .GlobalEnv)
            }
        }
    }
    .call <- function(...){return(.Call(...))};
    .c <- function(...){return(.C(...))};
    args <- list(model = model, dir = dir, prefix = prefix,
                 extraC = extraC, force = force, modName = modName,
                 ...);
    ret <- list(dll     = finalDll,
                model   = allModVars$model["model"],
                extra   = extraC,
                modVars = allModVars,
                .call   = .call,
                .c      = .c,
                args    = args);
    class(ret) <- "rxDll";
    return(ret);
}

##' @rdname rxCompile
##' @export
rxCompile.rxDll <- function(model, ...){
    args <- as.list(match.call(expand.dots = TRUE));
    rxDllArgs <- model$args;
    if (any(names(rxDllArgs) == "dir")){
        args$dir <- rxDllArgs$dir;
    }
    if (any(names(rxDllArgs) == "prefix")){
        args$prefix <- rxDllArgs$prefix;
    }
    if (any(names(rxDllArgs) == "extraC")){
        args$extraC <- rxDllArgs$extraC;
    }
    if (any(names(rxDllArgs) == "force")){
        args$force <- rxDllArgs$force;
    }
    if (any(names(rxDllArgs) == "modName")){
        args$modName <- rxDllArgs$modName;
    }
    args$model = rxDllArgs$model;
    return(do.call(getFromNamespace("rxCompile", "RxODE"), args, envir = parent.frame(1)));
}

##' @rdname rxCompile
##' @export
rxCompile.RxODE <- function(model, ...){
    model$cmpMgr$compile()
}


##' Return the DLL associated with the RxODE object
##'
##' This will return the dynamic load library or shared object used to
##' run the C code for RxODE.
##'
##' @param obj A RxODE family of objects or a character string of the
##'     model specification or location of a file with a model
##'     specification.
##'
##' @return a path of the library
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
rxDll <- function(obj, ...){
    UseMethod("rxDll");
}

##' @rdname rxDll
##' @export
rxDll.character <- function(obj, ...){
    return(rxDll(rxCompile(obj, ...)))
}

##' @rdname rxDll
##' @export
rxDll.rxDll <- function(obj, ...){
    return(obj$dll)
}

##' @rdname rxDll
##' @export
rxDll.RxODE <- function(obj, ...){
    return(rxDll(obj$cmpMgr$rxDll()))
}

##' Load the DLL for the object
##'
##' This loads the DLL into the current R session to allow C functions
##' to be called in R.
##'
##' @param obj a RxODE family of objects
##'
##' @author Matthew L.Fidler
##' @export
rxLoad <- function(obj){
    if (!(rxDllLoaded(obj))){
        dll <- obj$cmpMgr$dllfile;
        rc <- try(dyn.load(dll), silent = TRUE);
        if (inherits(rc, "try-error")){
            if (RxODE.compile.on.load){
                rxCompile(obj);
                rc <- try(dyn.load(dll), silent = TRUE);
                if (inherits(rc, "try-error")){ #nocov start
                    ## Should not get here.
                    ## rxDelete(obj);
                    obj <- RxODE(obj);
                    if (rxDllLoaded(obj)){
                        return(invisible());
                    }
                    stop(sprintf("Error loading DLL file '%s', even after trying to recompile.", dll));
                } # nocov end
            } else {
                stop(sprintf("Error loading DLL file '%s'.", dll));
            }
        }
    }
    return(invisible());
}

##' Unload the DLL for the object
##'
##' This unloads the DLL in the R session so that the DLL can be
##' deleted.  All the c functions will no longer be accessible.
##'
##' @param obj a RxODE family of objects
##'
##' @author Matthew L.Fidler
##' @export
rxUnload <- function(obj){
    if ((rxDllLoaded(obj))){
        dll <- rxDll(obj);
        rc <- try(dyn.unload(dll), silent = TRUE)
        if (inherits(rc, "try-error"))
            stop(sprintf("Error unloading DLL file '%s'.", dll));
    }
    return(invisible());
}

##' Delete the DLL for the model
##'
##' This function deletes the DLL, but doesn't delete the model
##' information in the object.
##'
##' @param obj RxODE family of objects
##'
##' @return A boolean stating if the operation was successful.
##'
##' @author Matthew L.Fidler
##' @export
rxDelete <- function(obj){
    if (class(obj) == "RxODE"){
        obj$delete();
    } else {
        dll <- rxDll(obj);
        cFile <- sub(rex::rex(.Platform$dynlib.ext,end), ".c", dll)
        if (file.exists(cFile)){
            unlink(cFile);
        }
        rxUnload(obj)
        unlink(dll);
        return(!file.exists(dll));
    }
}

##' Parameters specified by the model
##'
##' This return the model's parameters that are required to solve the
##' ODE system.
##'
##' @inheritParams rxModelVars
##'
##' @return a character vector listing the parameters in the model.
##'
##' @author Matthew L.Fidler
##' @export
rxParams <- function(obj){
    return(rxModelVars(obj)$params);
}

##' @rdname rxParams
##' @export
rxParam <- rxParams


##' Jacobain and parameter derivates
##'
##' Return Jacobain and parameter derivates
##'
##' @inheritParams rxModelVars
##'
##' @return A list of the jacobian parameters defined in this RxODE
##'     object.
##' @author Matthew L. Fidler
##' @export
rxDfdy <- function(obj){
    return(rxModelVars(obj)$dfdy);
}

##' State variables
##'
##' This returns the model's compartments or states.
##'
##' @inheritParams rxModelVars
##'
##' @param state is a string indicating the state or compartment that
##'     you would like to lookup.
##'
##' @return If state is missing, return a character vector of all the states.
##'
##' If state is a string, return the compartment number of the named state.
##'
##' @seealso \code{\link{RxODE}}
##'
##' @author Matthew L.Fidler
##' @export
rxState <- function(obj, state){
    if (missing(state)){
        return(rxModelVars(obj)$state);
    } else {
        objState <- rxState(obj);
        if (length(objState) == 1)
            warning("Only one state variable should be input.", immediate = TRUE);
        w <- which(objState == state)
        if (length(w) != 1){
            stop(sprintf("Cannot locate compartment \"%s\".", state));
        }
        return(w);
    }
}
##' Left handed Variables
##'
##' This returns the model calculated variables
##'
##' @inheritParams rxModelVars
##'
##' @return a character vector listing the calculated parameters
##' @seealso \code{\link{RxODE}}
##'
##' @author Matthew L.Fidler
##' @export
rxLhs <- function(obj){
    return(rxModelVars(obj)$lhs);
}

rxConditionLst <- list();
##' Current Condition for RxODE object
##'
##' @param obj RxODE object
##' @param condition If specified and is one of the conditions in the
##'     RxODE object (as determined by \code{\link{rxExpandIfElse}}),
##'     assign the RxODE current condition to this parameter.  If the
##'     condition is not one of the known condition, the condition is
##'     set to \code{NULL}, implying no conditioning currently used.
##' @return Current condition for RxODE object
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxCondition <- function(obj, condition=NULL){
    key <- digest::digest(rxNorm(obj,FALSE),algo="md5");
    if (!missing(condition) && is.null(condition)){
        condition <- FALSE;
    }
    if (is.null(condition)){
        return(rxConditionLst[[key]]);
    } else if (any(condition == rxNorm(obj,TRUE))){
        lst <- rxConditionLst;
        lst[[key]] <- condition;
        assignInMyNamespace("rxConditionLst", lst);
        return(rxConditionLst[[key]]);
    } else {
        lst <- rxConditionLst;
        lst[[key]] <- NULL;
        assignInMyNamespace("rxConditionLst", lst);
        return(rxConditionLst[[key]]);
    }
}

##' Get the normalized model
##'
##'
##' This get the syntax prefered model for processing
##'
##' @inheritParams rxModelVars
##' @param condition Character string of a logical condition to use
##'     for subsetting the normalized model.  When missing, and a
##'     condition is not set via \code{rxCondition}, return the whole
##'     code with all the conditional settings intact.  When a
##'     condition is set with \code{rxCondition}, use that condition.
##' @param removeInis A boolean indicating if paramter initilizations will be removed from the model
##' @param removeJac A boolean indicating if the Jacobians will be removed.
##' @param removeSens A boolean indicating if the sensitivities will be removed.
##' @return Normalized Normal syntax (no comments)
##' @author Matthew L. Fidler
##' @export
rxNorm <- function(obj, condition=NULL, removeInis, removeJac, removeSens){
    if (!missing(removeInis) || !missing(removeJac) || !missing(removeSens)){
        ret <- strsplit(rxNorm(obj, condition), "\n")[[1]];
        if (missing(removeInis)){
            removeInis <- FALSE
        }
        if (missing(removeJac)){
            removeJac <- FALSE
        }
        if (missing(removeSens)){
            removeSens <- FALSE
        }
        if (removeInis){
            ret <- rxRmIni(ret)
        }
        if (removeJac){
            ret <- rxRmJac(ret)
        }
        if (removeSens){
            ret <- rxRmSens(ret)
        }
        return(paste(ret, collapse="\n"))
    } else {
        if (class(condition) == "logical"){
            if (!condition){
                condition <- NULL;
            } else {
                tmp <- rxExpandIfElse(obj)
                return(names(tmp))
            }
        } else if (is.null(condition)){
            condition <- rxCondition(obj);
        }
        if (is.null(condition)){
            tmp <- rxModelVars(obj)$model["normModel"]
            names(tmp) <- NULL;
            return(tmp)
        } else {
            if (class(condition) == "character"){
                tmp <- rxExpandIfElse(obj)[condition];
                names(tmp) <- NULL;
                return(tmp)
            } else {
                return(rxNorm(obj, FALSE));
            }
        }
    }
}

##' All model variables for a RxODE object
##'
##' Return all the known model variables for a specified RxODE object
##'
##' These items are only calculated after compilation; they are
##' built-into the RxODE compiled DLL.
##'
##' @param obj RxODE family of objects
##'
##' @return A list of RxODE model properties including:
##'
##' \item{params}{ a character vector of names of the model parameters}
##' \item{lhs}{ a character vector of the names of the model calculated parameters}
##' \item{state}{ a character vector of the compartments in RxODE object}
##' \item{trans}{ a named vector of translated model properties
##'       including what type of jacobian is specified, the \code{C} function prefixes,
##'       as well as the \code{C} functions names to be called through the compiled model.}
##' \item{md5}{a named vector that gives the digest of the model (\code{file_md5}) and the parsed model
##'      (\code{parsed_md5})}
##' \item{model}{ a named vector giving the input model (\code{model}),
##'    normalized model (no comments and standard syntax for parsing, \code{normModel}),
##'    and interim code that is used to generate the final C file \code{parseModel}}
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
rxModelVars <- function(obj){
    UseMethod("rxModelVars");
}

##' @rdname rxModelVars
##' @export
rxModelVars.list <- function(obj){
    if (all(c("params", "lhs", "state", "trans", "ini", "model", "md5", "podo", "dfdy") %in% names(obj))){
        return(obj);
    } else {
        stop("Cannot figure out the model variables.")
    }
}

##' @rdname rxModelVars
##' @export
rxModelVars.rxDll <- function(obj){
    return(obj$modVars)
}

##' @rdname rxModelVars
##' @export
rxModelVars.RxCompilationManager <- function(obj){
    return(rxModelVars.rxDll(obj$rxDll()))
}

##' @rdname rxModelVars
##' @export
rxModelVars.RxODE <- function(obj){
    return(rxModelVars.rxDll(obj$cmpMgr$rxDll()))
}

##' @rdname rxModelVars
##' @export
rxModelVars.solveRxODE <- function(obj){
    env <- attr(obj, ".env");
    rxode <- env$env$out;
    return(rxModelVars.RxODE(rxode));
}

##' @rdname rxModelVars
##' @export
rxModelVars.character <- function(obj){
    if (length(obj) == 1){
        cFile <- tempfile();
        if (file.exists(obj)){
            parseModel <- obj;
            on.exit({unlink(cFile)});
        } else {
            parseModel <- tempfile();
            sink(parseModel);
            cat(paste(obj, collapse="\n"));
            sink()
            on.exit({unlink(parseModel); unlink(cFile)});
        }
        ret <- rxTrans(parseModel, cFile, modVars=TRUE);
        return(ret);
    } else {
        rxModelVars.character(paste(obj, collapse="\n"));
    }
}

rxModelVars.character.slow <- NULL;

##' Print rxDll object
##'
##' This tells if the rxDll is loaded, ready and/or deleted.
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
print.rxDll <- function(x, ...){
    if (file.exists(x$dll)){
        cat(sprintf("RxODE DLL named \"%s\"", basename(x$dll)));
        if (rxDllLoaded(x)){
            cat(" is loaded and ready to use.\n");
        } else {
            cat(" is not loaded now.\n");
        }
    } else {
        cat(sprintf("RxODE DLL named \"%s\" has been deleted.\n", basename(x$dll)));
    }
    invisible(x);
}

##' Summary of rxDll object
##'
##' This gives expanded information about the rxDll object
##'
##' @param object RxDll object
##'
##' @param ... Other arguments.  Includes \code{noprint}, which is a
##'     logical telling if the object should print the rxDll object
##'     first. By default this is FALSE
##'
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
summary.rxDll <- function(object, ...){
    args <- as.list(match.call(expand.dots = TRUE));
    if (any(names(args) == "noprint")){
        noprint <- args$noprint;
    } else {
        noprint <- FALSE;
    }
    if (!noprint)
        print(object);
    cat(sprintf("DLL: %s\n", rxDll(object)));
    cat(sprintf("Jacobian: %s\n", ifelse(rxModelVars(object)$jac == "fulluser", "Full User Specified", "Full Internally Caluclated")));
    print(coef(object));
    if (length(rxLhs(object)) > 0){
        cat("\nCalculated Variables:\n");
        print(rxLhs(object));
    }
    cat("\nModel:\n")
    cat(rxModelVars(object)$model["model"]);
    cat("\n");

    return(invisible(object))
}
##' Format theta and eta for parameter estimate values in RxODE
##'
##' @param theta A vector of theta estimates
##' @param eta A vector of eta estimates
##' @return A named vector for initial values.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxThetaEta <- function(theta=NULL, eta=NULL){
    ret <- c(theta, eta)
    ltheta <- length(theta);
    leta <- length(eta)
    if (ltheta > 0){
        ntheta <- paste0("THETA[",1:length(theta),"]")
    } else  {
        ntheta <- character()
        theta <- numeric()
    }
    if (leta > 0){
        neta <- paste0("ETA[", 1:length(eta), "]");
    } else {
        neta <- character();
        eta <- numeric();
    }
    ret <- c(theta, eta);
    names(ret) <- c(ntheta, neta)
    return(ret)
}

##' Initial Values and State values for a RxODE object
##'
##' Returns the initial values of the rxDll object
##'
##' @param rxDllObj rxDll, RxODE, or named vector representing default
##'     initial arguments
##'
##' @param vec If supplied, named vector for the model.
##'
##' @param req Required names, and the required order for the ODE solver
##'
##' @param default a number or NA representing the default value for
##'     parameters missing in \code{vec}, but required in \code{req}.
##'
##' @param noerror is a boolean specifying if an error should be thrown
##'     for missing parameter values when \code{default} = \code{NA}
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
rxInits <- function(rxDllObj,        # rxDll object
                    vec,             # Supplied parameters
                    req,             # Required names, and order
                    default = 0,     # Default value; If NA, then throw error if doesn't exist
                    noerror = FALSE, # no error thrown for NA
                    noini=FALSE,
                    ...){
    ## rxInits returns the inits of rxDllObj
    ##
    ## - When vec is not specified, returns the model specified
    ##   initialization values.
    ##
    ## - When vec is specified, replace rxDll inis with the values in
    ##   this named vector, and augment with any new variables.
    ##
    ## - When req is specified, make sure that the required variables
    ##   are included in the output in the order specified.
    ##
    ## - The default value specified, anything missing req names is
    ##   replaced with the default value.  If the default value is NA,
    ##   then throw an error if the values are not specified in either
    ##   the vec or the rxDllObj.
    if (noini){
        ini <- c();
    } else {
        ini <- rxModelVars(rxDllObj)$ini;
    }
    miss <- c();
    if (!missing(req)){
        if ((is.na(default) && noerror) || !is.na(default)){
            miss <- structure(rep(default, length(req)), .Names=req);
        }
    }
    if (!missing(vec)){
        if (length(vec) > 0){
            nv <- names(vec)
            if (is.null(nv)){
                if (!missing(req) && length(req) == length(vec)){
                    warning(sprintf("Assumed order of inputs: %s", paste(req, collapse = ", ")))
                    return(vec)
                } else {
                    stop(sprintf("Length mismatch\nreq: c(%s)\nvec: c(%s)\n%s", paste(req, collapse = ", "), paste(vec, collapse = ", "), rxModelVars(rxDllObj)))
                }
            } else {
                vec <- c(vec, ini, miss)
            }
        } else {
            vec <- c(ini, miss);
        }
    } else {
        vec <- c(ini, miss);
    }
    if (!(missing(req))){
        diff <- setdiff(req, names(vec))
        if (length(diff) && !noerror){
            stop(paste0("Missing the following parameter(s):\n",
                       paste(diff, collapse=" ")));
        } else {
            vec <- vec[req];
        }
    } else {
        vec <- vec[!duplicated(names(vec))]
    }
    return(vec);

} # end function rxInits

##' @rdname rxInits
##' @export
rxInit <- rxInits;

##' Reload RxODE DLL
##'
##' Can be useful for debugging
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxReload <- function(){
    tmp <- getLoadedDLLs()$RxODE
    class(tmp) <- "list";
    dyn.unload(tmp$path);
    ret <- is.null(getLoadedDLLs()$RxODE)
    dyn.load(tmp$path);
    ret <- ret && !is.null(getLoadedDLLs()$RxODE)
    return(ret)
}

##' Garbage Collection for RxODE objects
##'
##' When objects are created without a model name/DLL name and RxODE
##' option RxODE.delete.unnamed is TRUE, remove the associated DLL
##' when R is done with them (i.e. when the RxODE object is deleted or
##' R exits normally.)  If R crashes, thes files will still remain.
##'
##' @param env RxODE environment for cleanup.
##' @return Nothing
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxodeGc <- function(env){
    rxDelete(env$out);
}
