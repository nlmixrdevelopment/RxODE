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
##'     LHS variables when calculating sensitivities.
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
##' @importFrom mvnfast rmvn
##' @importFrom PreciseSums fsum
##' @importFrom Rcpp evalCpp
##' @importFrom utils getFromNamespace assignInMyNamespace download.file head sessionInfo
##' @importFrom stats setNames update
##' @importFrom methods signature is
##' @importFrom memoise memoise
##' @export
RxODE <- function(model, modName = basename(wd),
                  wd = ifelse(RxODE.cache.directory == ".", getwd(), RxODE.cache.directory),
                  filename = NULL, extraC = NULL, debug = FALSE, calcJac=NULL, calcSens=NULL,
                  collapseModel=FALSE, ...) {
    if (!missing(model) && !missing(filename))
        stop("Must specify exactly one of 'model' or 'filename'.")
    if (missing(model) && !missing(filename)){
        model <- filename;
    }
    if (!missing(model) && missing(filename)){
        if (is(substitute(model), "{")){
            model <- deparse(substitute(model));
            if (model[1] == "{"){
                model <- model[-1];
                model <- model[-length(model)];
            }
            model <- paste(model, collapse="\n");
        } else if (rxIs(model, "RxODE")){
            model <- RxODE::rxModelVars(model)$model["model"];
            if (!is.null(calcJac) && is.null(calcSens)){
                calcSens <- FALSE;
            }
        }
        ## else if ((is(model,"function") || is(model,"call"))){
        ##     model <- deparse(body(model))[-1];
        ##     model <- paste(model[-length(model)], collapse="\n");
        ## }
    }
    model <- RxODE::rxLinCmtTrans(model);
    class(model) <- "rxModelText"
    ## RxODE compilation manager (location of parsed code, generated C,  shared libs, etc.)
    .env <- new.env(parent=baseenv())
    .env$missing.modName <- missing(modName);
    wd <- .normalizePath(wd, "/", mustWork=F)
    if (.env$missing.modName){
        if (RxODE.tempfiles){
            .env$mdir <- .rxTempDir();
        } else {
            .env$mdir <- wd
        }
    } else {
        .env$mdir <- file.path(wd, sprintf("%s.d", modName));
    }

    if (!file.exists(wd))
        dir.create(wd, recursive = TRUE);

    .env$modName <- modName;
    .env$model <- model;
    .env$extraC <- extraC;
    .env$debug <- debug
    .env$calcJac <- calcJac;
    .env$calcSens <- calcSens;
    .env$collapseModel <- collapseModel;

    .env$wd <- wd;
    .env$compile <- eval(bquote(function(){
        with(.(.env), {
            model <- as.vector(model);
            .lwd <- getwd();
            if (!file.exists(wd))
                dir.create(wd, recursive = TRUE)

            if (!file.exists(wd))
                setwd(wd);
            on.exit(setwd(.lwd));
            if (missing.modName){
                assign("rxDll",
                       RxODE::rxCompile(model, extraC = extraC, debug = debug,
                                        calcJac=calcJac, calcSens=calcSens,
                                        collapseModel=collapseModel),
                       envir=.(.env));
            } else {
                assign("rxDll",
                       RxODE::rxCompile(model, dir=mdir, extraC = extraC,
                                        debug = debug, modName = modName,
                                        calcJac=calcJac, calcSens=calcSens,
                                        collapseModel=collapseModel),
                       envir=.(.env));
            }
        });
    }));
    .env$compile();
    .env$get.modelVars <- eval(bquote(function(){
        with(.(.env), {
            .mv <- RxODE::rxModelVars(rxDll);
            .ret <- RxODE::rxModelVars(rxDll)[c("params", "state", "lhs")];
            .p <- .ret["params"];
            .ini <- names(.mv$.ini)
            .init <- RxODE::rxInit(rxDll);
            .ret$params <- .ret$params[!(.ret$params %in% names(.init))]
            class(.ret) <- "list"
            return(.ret);
        })
    }))
    .env$state <- RxODE::rxModelVars(.env$rxDll)$state;
    .env$lhs <- RxODE::rxModelVars(.env$rxDll)$lhs;
    .env$params <- RxODE::rxModelVars(.env$rxDll)$params;
    .env$version <- RxODE::rxVersion()["version"];
    .env$solve <- eval(bquote(function(..., matrix=TRUE, object=NULL){
        RxODE::rxSolve(object=get("rxDll", envir=.(.env)), ..., matrix=matrix);
    }))
    .env$.c <- with(.env, function(...){.C(...)});
    .env$dll <- new.env(parent=baseenv())
    .env$dll$.c <- .env$.c;
    .env$assignPtr <- eval(bquote(function(){
        RxODE::rxAssignPtr(get("rxDll", envir=.(.env)));
    }))
    .env$run <- .env$solve;
    .env$modName <- modName;
    .env$model <- model; # actual model code
    ## cmpMgr = cmpMgr,
    .env$dynLoad <- eval(bquote(function(force = FALSE){
        rx <- .(.env);
        class(rx) <- "RxODE";
        RxODE::rxDynLoad(rx);
    }));
    .env$load <- .env$dynLoad;
    .env$dynUnload <- eval(bquote(function(){
        rx <- .(.env);
        class(rx) <- "RxODE";
        RxODE::rxDynUnload(rx);
    }));
    .env$unload <- .env$dynUnload;
    .env$isValid <- eval(bquote(function(){
        return(file.exists(RxODE::rxDll(get("rxDll", envir=.(.env)))));
    }))
    .env$isLoaded <- eval(bquote(function(){
        rx <- .(.env);
        class(rx) <- "RxODE";
        RxODE::rxIsLoaded(rx);
    }))
    .env$delete <- eval(bquote(function(){
        rx <- .(.env);
        class(rx) <- "RxODE";
        RxODE::rxDelete(rx);
    }))
    .env$parse <- with(.env, function(){
        stop("$parse is no longer supported");
    })
    .env$get.index <- eval(bquote(function(s){
        return(rxState(get("rxDll", envir=.(.env)), s));
    }))
    .mv <- RxODE::rxModelVars(.env$rxDll);
    .env$lib.name <- .mv$trans["lib.name"]
    tmp <- list(dllfile=RxODE::rxDll(.env$rxDll),
                ode_solver=as.vector(.mv$trans["ode_solver"]),
                ode_solver_ptr=as.vector(.mv$trans["ode_solver_ptr"]),
                prefix=as.vector(.mv$trans["prefix"]),
                model=model,
                isValid=eval(bquote(function(){with(.(.env), isValid())})),
                parse = eval(bquote(function(){with(.(.env), parse())})),
                compile = eval(bquote(function(){with(.(.env), compile())})),
                dynLoad = eval(bquote(function(){with(.(.env), dynLoad())})),
                dynUnload = eval(bquote(function(){with(.(.env), dynUnload())})),
                modelDir = .env$mdir, # model directory
                get.modelVars = eval(bquote(function(){with(.(.env), get.modelVars())})),
                delete = eval(bquote(function(){with(.(.env), delete())})),
                get.index = eval(bquote(function(...){with(.(.env), get.index(...))})),
                extraC        = .env$extraC,
                .rxDll         = .env$rxDll,
                rxDll=eval(bquote(function(){with(.(.env), return(rxDll))})));
    tmp <- list2env(tmp, parent=.env);
    class(tmp) <- "RxCompilationManager"
    .env$cmpMgr <- tmp;
    .env$calcJac <- (length(.mv$dfdy) > 0);
    .env$calcSens <- (length(.mv$sens) > 0)
    assign(.mv$trans["prefix"], .env, getFromNamespace(".rxModels", "RxODE"));
    class(.env) <- "RxODE"
    RxODE::rxForget();
    return(.env);
}

##' Get model properties without compiling it.
##'
##' @param model RxODE specification
##' @inheritParams RxODE
##' @return RxODE trans list
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxGetModel <- memoise::memoise(function(model, calcSens=FALSE, calcJac=FALSE, collapseModel=FALSE){
    if (is(substitute(model), "call")){
        model <- model;
    }
    if (is(substitute(model), "{")){
        model <- deparse(substitute(model))
        if (model[1] == "{"){
            model <- model[-1];
            model <- model[-length(model)];
        }
        model <- paste(model, collapse="\n");
    } else if (is(model, "function") || is(model, "call")){
        model <- deparse(body(model));
        if (model[1] == "{"){
            model <- model[-1];
            model <- model[-length(model)];
        }
        model <- paste(model, collapse="\n");
    } else if (is(model, "character")){
        if (file.exists(model)){
            ret$useModelName <- TRUE;
        }
    } else if (is(model, "name")){
        model <- eval(model);
    } else {
        stop(sprintf("Can't figure out how to handle the model argument (%s).", class(model)));
    }
    .parseModel <- tempfile("parseModel");
    .cFile <- tempfile("cFile");
    on.exit({unlink(.parseModel); unlink(.cFile)});
    sink(.parseModel);
    cat(model);
    cat("\n");
    sink();
    return(rxTrans(.parseModel, .cFile, calcSens=calcSens, calcJac=calcJac, collapseModel=collapseModel, modVars=TRUE));
})

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
    .args <- rev(as.list(match.call())[-1]);
    names(.args) <- c("obj", "solvedObject");
    return(do.call("rxChain2", .args, envir = parent.frame(1)));
}

##' @rdname rxChain
##' @export
'+.solveRxDll' <- function(obj1, obj2){
    return(RxODE::rxChain(obj1, obj2))
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
    .args <- as.list(match.call());
    stop(sprintf("Do not know how to add %s to RxODE solved object %s.", toString(.args[[2]]), toString(.args[[3]])))
}

##' @rdname rxChain2
##' @export
rxChain2.EventTable <- function(obj, solvedObject){
    .args <- rev(as.list(match.call())[-1]);
    names(.args) <- c("object", "events");
    return(do.call("rxSolve", .args, envir = parent.frame(1)));
}

.isLatex <- function() {
    if (!("knitr" %in% loadedNamespaces())) return(FALSE)
    get("is_latex_output", asNamespace("knitr"))()
}

.useUtf <- function() {
    opt <- getOption("cli.unicode", NULL)
    if (! is.null(opt)) {
        isTRUE(opt)
    } else {
        l10n_info()$`UTF-8` && !.isLatex()
    }
}
.getBound <- function(x, parent=parent.frame(2)){
    bound <- do.call("c", lapply(ls(parent), function(cur){
                              if (identical(parent[[cur]], x)){
                                  return(cur)
                              }
                              return(NULL);
                          }));
    if (length(bound) > 1) bound <- bound[1];
    if (length(bound) == 0){
        bound  <- ""
    }
    return(bound)
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
    .bound <- .getBound(x, parent.frame(2));
    .valid <- x$isValid()
    .msg2 <- "";
    if (!.valid){
        .msg <- crayon::red$bold("invalid");
        .msg2 <- paste0(" re-create with ", crayon::blue("RxODE::"), crayon::yellow("RxODE"))
        .ico <- crayon::red(cli::symbol$cross)
    } else {
        .loaded <- x$isLoaded();
        if (.loaded){
            .msg <- crayon::green$bold("ready")
            .ico <- crayon::green(cli::symbol$tick)
        } else{
            .msg <- crayon::yellow$bold("unloaded");
            .ico <- crayon::yellow(cli::symbol$warning);
            .msg2 <- paste0(" reload with ", crayon::blue("RxODE::"), crayon::yellow("rxLoad"));
        }
    }
    if (.useUtf()){
        .ico <- paste0(.ico, " ");
    } else {
        .ico <- ""
    }
    .dll <- basename(RxODE::rxDll(x));
    .dll <- substr(.dll, 1, nchar(.dll) - nchar(.Platform$dynlib.ext) - nchar(.Platform$r_arch) - 1)
    message(paste0(crayon::bold("RxODE "), as.vector(RxODE::rxVersion()["version"]), " model named ",
                   crayon::yellow$bold(.dll), " model (", .ico, .msg,
                   .msg2, ")."))
    if (!any(names(list(...)) == "rxSuppress") && .valid){
        .cur <- RxODE::rxState(x);
        if (length(.cur) > 0)
            message(paste0(crayon::yellow(.bound), crayon::blue$bold("$state"), ": ", paste(.cur, collapse=", ")))
        .cur <- RxODE::rxParams(x);
        if (length(.cur) > 0)
            message(paste0(crayon::yellow(.bound), crayon::blue$bold("$params"), ": ", paste(.cur, collapse=", ")))
        .cur <- RxODE::rxLhs(x);
        if (length(.cur) > 0)
            message(paste0(crayon::yellow(.bound), crayon::blue$bold("$lhs"), ": ", paste(.cur, collapse=", ")))
    }
    invisible(x)
}

##'@export
print.rxModelVars <- function(x, ...)
{
    .bound <- .getBound(x, parent.frame(2));
    message("RxODE model variables (see str to see all variables)");
    .cur <- x$state;
    if (length(.cur) > 0)
        message(paste0(crayon::yellow(.bound), crayon::blue$bold("$state"), ": ", paste(.cur, collapse=", ")))
    .cur <- x$params;
    if (length(.cur) > 0)
        message(paste0(crayon::yellow(.bound), crayon::blue$bold("$params"), ": ", paste(.cur, collapse=", ")))
    .cur <- x$lhs;
    if (length(.cur) > 0)
        message(paste0(crayon::yellow(.bound), crayon::blue$bold("$lhs"), ": ", paste(.cur, collapse=", ")))
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
##'
##' @export
coef.RxODE <- function(object,
                       ...){
    .ret <- RxODE::rxModelVars(object)[c("params", "state", "ini", "sens", "fn.ini")];
    .ret$RxODE <- object;
    class(.ret) <- "rxCoef";
    return(.ret);
}

##' @rdname coef.RxODE
##' @export
coef.RxCompilationManager <- function(...){
    coef.RxODE(...);
}

##' @rdname coef.RxODE
##' @export
coef.solveRxODE <- function(object, ...){
    .env <- attr(object, ".env");
    .rxode <- .env$env$out;
    .lst <- c(list(params=.env$params, inits=.env$inits), .env$extra.args)
    ret <- list();
    ret$params <- .lst$params;
    ret$state <- .lst$inits;
    ret$sens <- RxODE::rxModelVars(.rxode)["sens"];
    ret$state <- ret$state[regexpr(getFromNamespace("regSens", "RxODE"), names(ret$state)) == -1]
    ret$RxODE <- .rxode;
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
    .rxDllObj <- x$RxODE;
    if (length(rxParams(.rxDllObj)) > 0){
        cat(cli::rule(left="User supplied parameters:"), "\n");
        print(RxODE::rxInits(.rxDllObj, c(), RxODE::rxParams(.rxDllObj), NA, TRUE))
        cat(cli::rule(left="User initial conditions:"), "\n");
        .tmp <- RxODE::rxInits(.rxDllObj, c(), RxODE::rxState(.rxDllObj), 0, TRUE);
        if (length(x$sens) > 0){
            .tmp <- .tmp[regexpr(getFromNamespace("regSens", "RxODE"), names(.tmp)) == -1];
        }
        .tmp <- .tmp[!(names(.tmp) %in% x$fn.ini)];
        print(.tmp);
    }
    if (length(x$fn.ini) > 0){
        cat(cli::rule(left="Parameter-based initial conditions:"), "\n");
        print(x$fn.ini);
    }
    cat(cli::rule(left="Compartments:"), "\n");
    .tmp <- RxODE::rxState(.rxDllObj);
    if (length(.tmp) > 0){
        names(.tmp) <- paste0("cmt=", 1:length(.tmp));
        if (length(x$sens) > 0){
            .tmp1 <- .tmp[regexpr(getFromNamespace("regSens", "RxODE"), .tmp) == -1];
            print(.tmp1);
            cli::rule(left="Sensitivities:")
            .tmp2 <- gsub(getFromNamespace("regSens", "RxODE"), "d/dt(d(\\1)/d(\\2))",
                          .tmp[regexpr(regSens, .tmp) != -1]);
            print(.tmp2);
        } else {
            print(.tmp);
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
    if (length(RxODE::rxInits(rxDllObj)) > 0){
        cat("\nDefault parameter values:\n")
        print(RxODE::rxInits(rxDllObj));
    }
    cat("\nCompartments:\n");
    tmp <- RxODE::rxState(rxDllObj);
    names(tmp) <- paste0("cmt=", 1:length(tmp));
    sens <- tmp[regexpr(getFromNamespace("regSens", "RxODE"), tmp) != -1];
    if (length(sens) > 0){
        tmp <- tmp[regexpr(getFromNamespace("regSens", "RxODE"), tmp) == -1]
    }
    print(tmp);
    if (length(sens) > 0){
        sens <- gsub(getFromNamespace("regSens", "RxODE"), "d/dt(d(\\1)/d(\\2))", sens);
        cat("\nSensitivities:\n");
        print(sens)
    }
    return(invisible());
}

.rxPrefix <- function(model,          # Model or file name of model
                      modName = NULL, # Model name, overrides calculated model name.
                      calcJac=NULL,
                      calcSens=NULL,
                      collapseModel=FALSE,
                      ...){
    ## rxPrefix returns a prefix for a model.
    if (!is.null(modName)){
        .modelPrefix <- sprintf("%s_", gsub("\\W", "_", modName));
    } else if (file.exists(model)){
        .modelPrefix <- sprintf("%s_", gsub("\\W", "_", gsub("[.].*$", "", basename(model))));
    } else {
        .parseModel <- tempfile("parseModel2");
        .cFile <- tempfile("cFile2");
        on.exit({unlink(.parseModel); unlink(.cFile)});
        sink(.parseModel);
        cat(model);
        cat("\n");
        sink();
        .trans <- RxODE::rxTrans(.parseModel, .cFile, calcSens=calcSens, calcJac=calcJac, collapseModel=collapseModel)
        .modelPrefix <- sprintf("rx_%s_", .trans["parsed_md5"]);
    }
    .modelPrefix <- sprintf("%s%s_", .modelPrefix, .Platform$r_arch);
    return(.modelPrefix);
}

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
    if (is(model, "character")){
        if (length(model) == 1){
            if (file.exists(model)){
                .ret <- suppressWarnings({readLines(model)});
                .mod <- paste(.ret, collapse = "\n");
            } else {
                .ret <- model
                .mod <- model;
            }
        } else {
            .ret <- model;
            .mod <- paste(.ret, collapse="\n");
        }
        if (is(extraC, "character")){
            if (file.exists(extraC)){
                .ret <- c(.ret, gsub(rex::rex(or(any_spaces, any_newlines)), "", readLines(extraC), perl = TRUE));
            }
        }
        .tmp <- c(RxODE.syntax.assign, RxODE.syntax.star.pow, RxODE.syntax.require.semicolon, RxODE.syntax.allow.dots,
                  RxODE.syntax.allow.ini0, RxODE.syntax.allow.ini, RxODE.calculate.jacobian,
                  RxODE.calculate.sensitivity, calcJac, calcSens, collapseModel);
        .ret <- c(.ret, .tmp);
        .tmp <- getLoadedDLLs()$RxODE;
        class(.tmp) <- "list";
        ## new RxODE DLLs gives different digests.
        .ret <- c(.ret, digest::digest(.tmp$path, file=TRUE, algo="md5"));
        ## Add version and github repository information
        .ret <- c(.ret, RxODE::rxVersion());
        return(list(text = .mod,
                    digest = digest::digest(.ret, algo="md5")));
    } else {
        RxODE::rxModelVars(model)$md5;
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
    .mv <- RxODE::rxModelVars(model)
    if (modVars){
        return(.mv);
    } else {
        return(c(.mv$trans, .mv$md5));
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
        modelPrefix <- .rxPrefix(model, modName, calcSens, calcJac, collapseModel);
    }
    if (file.exists(model)){
        if (missing(md5)){
            md5 <- RxODE::rxMd5(model, extraC, calcJac, calcSens, collapseModel)$digest;
        }
    } else {
        stop("This only translates a file (currently; Try rxCompile).");
    }
    .parseModel <- tempfile("parseModel3");
    .out3 <- tempfile("out3a");
    on.exit(unlink(.parseModel));
    RxODE::rxReq("dparser");
    .ret <- try(.Call(trans, model, model, cFile, extraC, modelPrefix, md5, .parseModel, .out3, PACKAGE="RxODE"));
    if (inherits(.ret, "try-error")){
        message("Model")
        message(suppressWarnings(readLines(model)))
        stop("Cannot Create RxODE model");
    }
    md5 <- c(file_md5 = md5, parsed_md5 = rxMd5(c(.ret$model["normModel"],
                                                  .ret$ini,
                                                  .ret$state,
                                                  .ret$params,
                                                  .ret$lhs), extraC, calcJac, calcSens, collapseModel)$digest);
    .ret$md5 <- md5
    ## dparser::dpReload();
    ## rxReload()
    if (file.exists(cFile)){
        .ret$md5 <- md5
        if (is(calcSens, "logical")){
            if (!calcSens){
                calcSens <- NULL;
            }
        }
        if (!is.null(calcSens)){
            if (length(rxState(.ret)) <= 0){
                stop("Sensitivities do not make sense for models without ODEs.")
            }
            .new <- rxSymPySensitivity(RxODE::rxModelVars(rxNorm(.ret)), calcSens=calcSens, calcJac=calcJac,
                                      collapseModel=collapseModel);
            .expandModel <- tempfile("expandModel");
            sink(.expandModel);
            cat(.new);
            cat("\n");
            sink()
            .ret <- try({.Call(trans, model, .expandModel, cFile, extraC, modelPrefix, md5, .parseModel, .out3)},
                        silent=TRUE);
            md5 <- c(file_md5 = md5, parsed_md5 = rxMd5(.ret$model["normModel"], extraC, calcJac, calcSens,
                                                        collapseModel)$digest)
            .ret$md5 <- md5
            if (inherits(.ret, "try-error")){
                RxODE::rxCat("\n", paste(readLines(.expandModel), collapse="\n"), "\n");
                stop("Syntax Errors")
            }
            ## dparser::dpReload();
            ## rxReload();
            unlink(.expandModel);
            .ret$md5 <- md5;
        } else if (calcJac){
            if (length(rxState(.ret)) <= 0){
                stop("Jacobians do not make sense for models without ODEs.")
            }
            .new <- .rxSymPyJacobian(RxODE::rxModelVars(rxNorm(.ret)));
            .expandModel <- tempfile("expandModel2");
            sink(.expandModel);
            cat(.new);
            cat("\n");
            sink()
            .ret <- .Call("trans", model, .expandModel, cFile, extraC, modelPrefix, md5, .parseModel, .out3,
                          PACKAGE="RxODE");
            md5 <- c(file_md5 = md5, parsed_md5 = rxMd5(.ret$model["normModel"], extraC, calcJac, calcSens,
                                                        collapseModel)$digest)
            .ret$md5 <- md5
            ## dparser::dpReload();
            ## rxReload();
            unlink(.expandModel);
            .ret$md5 <- md5;
        }
        if (modVars){
            return(.ret)
        } else {
            return(c(.ret$trans, .ret$md5));
        }
    } else {
        RxODE::rxCat(model);
        stop("Syntax Error (see above)");
    }
}

.rxTransMakevars <- function(rxProps, # rxTrans translation properties
                            rxDll, # Dll of file
                            cFile,
                            ## List of compile flags
                            compileFlags =c("parsed_md5", "ode_solver", "ode_solver_sexp",
                                            "ode_solver_ptr", "ode_solver_xptr", "inis",
                                            "model_vars", "calc_lhs", "calc_jac", "dydt", "dydt_liblsoda",
                                            "dydt_lsoda", "calc_jac_lsoda", "ode_solver_solvedata",
                                            "ode_solver_get_solvedata", "neq", "nlhs", "fix_inis"),
                            debug        = FALSE,
                            ...){
    ## rxTransCompileFlags returns a string for the compiler options
    .ctxt <- suppressWarnings({readLines(cFile)});
    .neededProps <- c("jac", compileFlags);
    if (all(.neededProps %in% names(rxProps))){
        .ret <- "";
        .tmp <- rxProps[compileFlags];
        for (x in c("parsed_md5", "ode_solver", "ode_solver_sexp", "ode_solver_ptr", "ode_solver_xptr",
                    "ode_solver_focei_outer", "model_vars")){
            .tmp[sprintf("%s_str", x)] <- sprintf("\"%s\"", .tmp[x]);
        }
        .tmp["lib_str"] <- sprintf("\"%s\"", gsub(.Platform$dynlib.ext, "", basename(rxDll)));
        for (.flg in names(.tmp)){
            .ctxt <- gsub(sprintf("__%s__", toupper(.flg)), .tmp[.flg], .ctxt);
        }
        .ctxt <- gsub("__R_INIT__", sprintf("R_init_%s", gsub(.Platform$dynlib.ext, "", basename(rxDll))), .ctxt)
        .ctxt <- gsub("__R_UNLOAD__", sprintf("R_unload_%s", gsub(.Platform$dynlib.ext, "", basename(rxDll))), .ctxt)
        .ctxt <- gsub("__TIMEID__", paste0(as.integer(Sys.time())), .ctxt)
        writeLines(.ctxt, cFile);
        .ret <- ""
        if (debug){
            .ret <- sprintf("%s -D__DEBUG__", .ret);
        }
        .ret <- sprintf("#RxODE Makevars\nPKG_CFLAGS=%s\nPKG_LIBS=$(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)", .ret);
        ## .ret <- paste(.ret, "-g");
        cat(.ret);
        return(.ret);
    } else {
        message(sprintf("Needed Variables: %s", paste(.neededProps, collapse=",")));
        print(setdiff(.neededProps, names(rxProps)));
        print(setdiff(names(rxProps), .neededProps));
        stop(sprintf("Cannot compile, only found %s.", paste(names(rxProps), collapse=",")));
    }
} # end function rxTransCompileFlags

##' @rdname rxIsLoaded
##' @export
rxDllLoaded <- rxIsLoaded
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
                                 dir=NULL,             # Directory
                                 prefix=NULL,     # Prefix
                                 extraC  = NULL,  # Extra C File.
                                 force   = FALSE, # Force compile
                                 modName = NULL,  # Model Name
                                 calcJac=NULL, # Calculate Jacobian
                                 calcSens=NULL, # Calculate Sensitivity
                                 collapseModel=FALSE,
                                 ...){
    ## rxCompile returns the DLL name that was created.
    if (is.null(dir)){
        if (RxODE.cache.directory != "."){
            dir <- RxODE.cache.directory;
        } else if (RxODE.tempfiles){
            dir <- .rxTempDir()
        } else {
            dir <- getwd();
        }
    } else {
        if (.Platform$OS.type == "windows"){
            dir <- suppressMessages(.normalizePath(dir, mustWork=FALSE));
        } else {
            dir <- file.path(getwd(), dir);
        }
    }
    if (!file.exists(dir))
        dir.create(dir, recursive = TRUE)
    if (is.null(prefix)){
        prefix <- .rxPrefix(model, modName, calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel);
    }
    .cFile <- file.path(dir, sprintf("%s.c", substr(prefix, 0, nchar(prefix)-1)));
    .cDllFile <- file.path(dir, sprintf("%s%s", substr(prefix, 0, nchar(prefix)-1), .Platform$dynlib.ext));
    if (file.exists(model)){
        .mFile <- .normalizePath(model);
    } else {
        .mFile <- sprintf("%s.rx", substr(.cFile, 0, nchar(.cFile)-2));
        sink(.mFile);
        cat(model);
        cat("\n");
        sink();
    }
    on.exit(unlink(.mFile));
    .md5 <- RxODE::rxMd5(.mFile, extraC, calcJac, calcSens, collapseModel);
    .allModVars <- NULL;
    .needCompile <- TRUE
    if (file.exists(.cDllFile)){
        try(dynLoad(.cDllFile), silent = TRUE);
        .modVars <- sprintf("%smodel_vars", prefix);
        if (is.loaded(.modVars)){
            .allModVars <- eval(parse(text = sprintf(".Call(\"%s\")", .modVars)), envir = .GlobalEnv)
            .modVars <- .allModVars$md5;
        }
        if (!any(names(.modVars) == "file_md5")){
            .needCompile <- FALSE;
        } else {
            if (.modVars["file_md5"] == .md5$digest){
                .needCompile <- FALSE;
            }
        }
    }
    if (force || .needCompile){
        .Makevars <- file.path(dir, "Makevars");
        if (file.exists(.Makevars)){
            if ("#RxODE Makevars" == readLines(.Makevars, 1)){
                unlink(.Makevars)
                on.exit({if (file.exists(.Makevars)){unlink(.Makevars)}}, add=TRUE);
            } else {
                file.rename(.Makevars, paste0(.Makevars, ".bakrx"));
                on.exit({
                    if (file.exists(.Makevars)){
                        unlink(.Makevars)
                    };
                    file.rename(paste0(.Makevars, ".bakrx"), .Makevars)
                }, add=TRUE)
            }
        } else {
            on.exit({if (file.exists(.Makevars)){unlink(.Makevars)}}, add=TRUE);
        }
        .trans <- rxTrans(.mFile, cFile = .cFile, md5 = .md5$digest, extraC = extraC, ..., modelPrefix = prefix,
                         calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel);
        if (file.exists(.cDllFile)){
            if (.modVars["parsed_md5"] == .trans["parsed_md5"]){
                RxODE::rxCat("Don't need to recompile, minimal change to model detected.\n");
                .needCompile <- FALSE;
            }
        }
        if (force || .needCompile){
            ## Setup Makevars
            .owd <- getwd();
            on.exit(setwd(.owd), add=TRUE);
            ## Now create C file
            .mv <- rxTrans(.mFile, cFile = .cFile, md5 = .md5$digest, extraC = extraC, ...,
                          modelPrefix = prefix, calcJac=calcJac, calcSens=calcSens, collapseModel=collapseModel,
                          modVars=TRUE);
            .j <- 0;
            .i <- 0;
            if (length(.mv$ini) > 0){
                .fixInis <- sprintf("double _theta[%d];\n%s\n", length(.mv$params),
                                   paste(sapply(.mv$params, function(x){
                                       if (!is.na(.mv$ini[x])){
                                           ret <- sprintf("  _theta[%d] = %.16f;", .i, as.vector(.mv$ini[x]));
                                           .i <<- .i + 1;
                                           return(ret)
                                       } else {
                                           ret <- sprintf("  _theta[%d] = theta[%d];", .i, .j);
                                           .i <<- .i + 1;
                                           .j <<- .j + 1;
                                           return(ret);
                                       }
                                   }), collapse="\n"))
            } else {
                .fixInis <- sprintf("double *_theta = theta;");
            }
            .trans <- c(.mv$trans, .mv$md5);
            .trans["fix_inis"] <- .fixInis;
            sink(.Makevars);
            cat(.rxTransMakevars(.trans, .cDllFile, .cFile, ...));
            sink();
            .sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths
            ## Change working directory
            setwd(dir);
            try(dyn.unload(.cDllFile), silent = TRUE);
            try(unlink(.cDllFile));
            .cmd <- sprintf("%s/bin/R CMD SHLIB %s",
                           Sys.getenv("R_HOME"), basename(.cFile));
            ## message(.cmd);
            do.call(.sh, list(.cmd, ignore.stdout=TRUE, ignore.stderr=TRUE));
            .tmp <- try(dynLoad(.cDllFile));
            if (inherits(.tmp, "try-error")){
                stop("Error loading model.")
            }
            .modVars <- sprintf("%smodel_vars", prefix);
            if (is.loaded(.modVars)){
                .allModVars <- eval(parse(text = sprintf(".Call(\"%s\")", .modVars)), envir = .GlobalEnv)
            }
        }
    }
    .call <- function(...){return(.Call(...))};
    .c <- function(...){return(.C(...))};
    .args <- list(model = model, dir = dir, prefix = prefix,
                 extraC = extraC, force = force, modName = modName,
                 ...);
    ret <- suppressWarnings({list(dll     = .cDllFile,
                                  c       = .cFile,
                                  model   = .allModVars$model["model"],
                                  extra   = extraC,
                                  modVars = .allModVars,
                                  .call   = .call,
                                  .c      = .c,
                                  args    = .args)});
    class(ret) <- "rxDll";
    return(ret);
}

##' @rdname rxCompile
##' @export
rxCompile.rxDll <- function(model, ...){
    .args <- as.list(match.call(expand.dots = TRUE));
    .rxDllArgs <- model$args;
    if (any(names(.rxDllArgs) == "dir")){
        .args$dir <- .rxDllArgs$dir;
    }
    if (any(names(.rxDllArgs) == "prefix")){
        .args$prefix <- .rxDllArgs$prefix;
    }
    if (any(names(.rxDllArgs) == "extraC")){
        .args$extraC <- .rxDllArgs$extraC;
    }
    if (any(names(.rxDllArgs) == "force")){
        .args$force <- .rxDllArgs$force;
    }
    if (any(names(.rxDllArgs) == "modName")){
        .args$modName <- .rxDllArgs$modName;
    }
    .args$model <- .rxDllArgs$model;
    return(do.call(getFromNamespace("rxCompile", "RxODE"), .args, envir = parent.frame(1)));
}

##' @export
print.rxC <- function(x, ...){
    message(sprintf("C File: %s  (summary for code)", x));
}

##' @export
summary.rxC <- function(object, ...){
    message(sprintf("//C File: %s", object));
    message("//");
    suppressWarnings(message(paste(readLines(object), collapse="\n")));
}

##' @rdname rxCompile
##' @export
rxCompile.RxODE <- function(model, ...){
    model$compile()
}

##' @rdname rxDynLoad
##' @export
rxLoad <- rxDynLoad

##' @rdname rxDynUnload
##' @export
rxUnload <- rxDynUnload

.rxConditionLst <- list();
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
    .key <- digest::digest(RxODE::rxNorm(obj, FALSE), algo="md5");
    if (!missing(condition) && is.null(condition)){
        condition <- FALSE;
    }
    if (is.null(condition)){
        return(getFromNamespace(".rxConditionLst", "RxODE")[[.key]]);
    } else if (any(condition == rxNorm(obj, TRUE))){
        .lst <- getFromNamespace(".rxConditionLst", "RxODE");
        .lst[[.key]] <- condition;
        assignInMyNamespace(".rxConditionLst", .lst);
        return(getFromNamespace(".rxConditionLst", "RxODE")[[.key]]);
    } else {
        .lst <- getFromNamespace(".rxConditionLst", "RxODE");
        .lst[[.key]] <- NULL;
        assignInMyNamespace(".rxConditionLst", .lst);
        return(getFromNamespace(".rxConditionLst", "RxODE")[[.key]]);
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
        .ret <- strsplit(rxNorm(obj, condition), "\n")[[1]];
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
            .ret <- .rxRmIni(.ret)
        }
        if (removeJac){
            .ret <- .rxRmJac(.ret)
        }
        if (removeSens){
            .ret <- .rxRmSens(.ret)
        }
        return(paste(.ret, collapse="\n"))
    } else {
        if (is(condition, "logical")){
            if (!condition){
                condition <- NULL;
            } else {
                .tmp <- RxODE::rxExpandIfElse(obj)
                return(names(.tmp))
            }
        } else if (is.null(condition)){
            condition <- RxODE::rxCondition(obj);
        }
        if (is.null(condition)){
            .tmp <- RxODE::rxModelVars(obj)$model["normModel"]
            names(.tmp) <- NULL;
            return(.tmp)
        } else {
            if (is(condition, "character")){
                .tmp <- RxODE::rxExpandIfElse(obj)[condition];
                names(.tmp) <- NULL;
                return(.tmp)
            } else {
                return(rxNorm(obj, FALSE));
            }
        }
    }
}

.rxModelVarsCharacter <- memoise::memoise(function(obj){
    if (length(obj) == 1){
        .cFile <- tempfile("cFile3");
        if (file.exists(obj)){
            .parseModel <- obj;
            on.exit({unlink(.cFile)});
        } else {
            .parseModel <- tempfile("parseModel4");
            sink(.parseModel);
            cat(paste(obj, collapse="\n"));
            sink()
            on.exit({unlink(.parseModel); unlink(.cFile)});
        }
        .ret <- rxTrans(.parseModel, .cFile, modVars=TRUE);
        return(.ret);
    } else {
        .rxModelVarsCharacter(paste(obj, collapse="\n"));
    }
})


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
    .args <- as.list(match.call(expand.dots = TRUE));
    if (any(names(.args) == "noprint")){
        .noprint <- .args$noprint;
    } else {
        .noprint <- FALSE;
    }
    if (!.noprint)
        print(object);
    cat(sprintf("DLL: %s\n", RxODE::rxDll(object)));
    cat(sprintf("Jacobian: %s\n",
                ifelse(RxODE::rxModelVars(object)$jac == "fulluser", "Full User Specified",
                       "Full Internally Calculated")));
    print(coef(object));
    if (length(RxODE::rxLhs(object)) > 0){
        cat("\nCalculated Variables:\n");
        print(RxODE::rxLhs(object));
    }
    .tmp <- as.vector(RxODE::rxModelVars(object)$model["model"])
    class(.tmp) <- "rxModelText"
    print(.tmp)
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
    .ret <- c(theta, eta)
    .ltheta <- length(theta);
    leta <- length(eta)
    if (.ltheta > 0){
        .ntheta <- paste0("THETA[", 1:length(theta), "]")
    } else  {
        .ntheta <- character()
        theta <- numeric()
    }
    if (leta > 0){
        .neta <- paste0("ETA[", 1:length(eta), "]");
    } else {
        .neta <- character();
        eta <- numeric();
    }
    .ret <- c(theta, eta);
    names(.ret) <- c(.ntheta, .neta)
    return(.ret)
}

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
    .tmp <- getLoadedDLLs()$RxODE
    class(.tmp) <- "list";
    dyn.unload(.tmp$path);
    .ret <- is.null(getLoadedDLLs()$RxODE)
    dynLoad(.tmp$path);
    .ret <- .ret && !is.null(getLoadedDLLs()$RxODE)
    return(.ret)
}

.rxModels <- new.env(parent = emptyenv())
##' Get the rxModels  information
##'@param env boolean that returns the environment where models are stored (TRUE), or the currently assigned RxODE model variables (FALSE).
##'@keywords internal
##'@export
rxModels_ <- function(env=TRUE){
    if (env){
        return(getFromNamespace(".rxModels", "RxODE"));
    } else {
        return(.Call(RxODE_get_mv, PACKAGE="RxODE"))
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
##' @author Matthew L. Fidler
##' @export
rxModelVars <- function(obj){
    .tmp <- try(obj, silent=TRUE);
    if (inherits(.tmp, "try-error")){
        obj <- as.character(substitute(obj));
    }
    rxModelVars_(obj);
}


##' Parameters specified by the model
##'
##' This return the model's parameters that are required to solve the
##' ODE system.
##'
##' @inheritParams rxModelVars
##'
##' @param constants is a boolean indicting if constants should be
##'     included in the list of parameters. Currently RxODE parses
##'     constants into variables in case you wish to change them
##'     without recompiling the RxODE model.
##'
##' @return a character vector listing the parameters in the model.
##'
##' @author Matthew L.Fidler
##' @export
rxParams <- function(obj, constants=TRUE){
    .ret <- rxParams_(obj)
    if (!constants){
        .init <- RxODE::rxInit(obj);
        .ret <- .ret[!(.ret %in% names(.init))]
    }
    return(.ret);
}

##' @rdname rxParams
##' @export
rxParam <- rxParams
