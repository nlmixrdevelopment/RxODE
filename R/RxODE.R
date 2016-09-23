rex::register_shortcuts("RxODE");
regIni <- rex::rex(or(group(one_of("_."),"0"),"0","(0)","[0]","{0}"),end);

.onLoad <- function(libname,pkgname){
    ## Setup RxODE.prefer.tbl
    op <- options();
    op.rx <- list(RxODE.prefer.tbl = FALSE);
    w <- !(names(op.rx) %in% names(op))
    if(any(w)) options(op.rx[w]);
}

#' Create an ODE-based model specification 
#'
#' Create a dynamic ODE-based model object suitably for translation
#' into fast C code
#'
#' @param model a string containing the set of ordinary differential
#'     equations (ODE) and other expressions defining the changes in
#'     the dynamic system (see also the \code{filename} argument).For
#'     details, see the sections \dQuote{Details} and
#'     \dQuote{\code{RxODE Syntax}} below.
#'
#' @param modName a string to be used as the model name. This string
#'     is used for naming various aspects of the computations,
#'     including generating C symbol names, dynamic libraries,
#'     etc. Therefore, it is necessary that \code{modName} consists of
#'     simple ASCII alphanumeric characters starting with a letter.
#'
#' @param wd character string with a working directory where to create
#'     a subdirectory according to \code{modName} (defaults to the
#'     current working directory \code{getwd()}. A subdirectoy named
#'     after the \dQuote{\code{modName.d}} will be created and
#'     populated with a C file, a dynamic loading library, plus
#'     various other working files.
#'
#' @param filename A file name or connection object where the
#'     ODE-based model specification resides. Only one of \code{model}
#'     or \code{filename} may be specified.
#'
#' @param do.compile logical specifying whether to carry out the
#'     parsing of the ODE into C code and its compilation. Default is
#'     \code{TRUE}
#' 
#' @param extraC a string indicating the path to a file with extra c
#'     functions needed for the compiled ODE model.
#'
#' @param debug is a boolean indicating if the executable should be
#'     compiled with verbose debugging information turned on.
#' 
#' @param ... any other arguments are passed to the function
#'     \code{\link{readLines}}, (e.g., encoding).
#'
#' The \dQuote{Rx} in the name \code{RxODE} is meant to suggest the
#' abbreviation \emph{Rx} for a medical prescription, and thus to
#' suggest the package emphasis on pharmacometrics modeling, including
#' pharmacokinetics (PK), pharmacodynamics (PD), disease progression,
#' drug-disease modeling, etc.
#'
#' The ODE-based model specification may be coded inside a character
#' string or in a text file, see Section \emph{RxODE Syntax} below for
#' coding details.  An internal \code{RxODE} compilation manager
#' object translates the ODE system into C, compiles it, and
#' dynamically loads the object code into the current R session.  The
#' call to \code{RxODE} produces an object of class \code{RxODE} which
#' consists of a list-like structure (closure) with various member
#' functions (see Section \emph{Value} below).
#'
#' For evaluating \code{RxODE} models, two types of inputs may be
#' provided: a required set of time points for querying the state of
#' the ODE system and an optional set of doses (input amounts).  These
#' inputs are combined into a single \emph{event table} object created
#' with the function \code{\link{eventTable}}.
#'
#' @section RxODE Syntax:
#'
#' An \code{RxODE} model specification consists of one or more
#' statements terminated by semi-colons, \sQuote{\code{;}}, and
#' optional comments (comments are delimited by \code{#} and an
#' end-of-line marker).  \strong{NB:} Comments are not allowed inside
#' statements.
#' 
#' A block of statements is a set of statements delimeted by curly
#' braces, \sQuote{\code{\{ ... \}}}. Statements can be either
#' assignments or conditional \code{if} statements. Assignment
#' statements can be: (1) \dQuote{simple} assignmets, where the left
#' hand is an identifier (i.e., variable), (2) special
#' \dQuote{time-derivative} assignments, where the left hand specifies
#' the change of that variable with respect to time e.g.,
#' \code{d/dt(depot)}, or (3) special \dQuote{jacobian} assignments,
#' where the left hand specifies the change of of the ODE with respect
#' to one of the parameters, e.g. \code{df(depot)/dy(kel)}.  The
#' \dQuote{jacobian} assignments are not required, and are only useful
#' for very stiff differential systems.  
#'
#' Expressions in assignment and \sQuote{\code{if}} statements can be
#' numeric or logical (no character expressions are currently
#' supported). Numeric expressions can include the following numeric
#' operators (\sQuote{\code{+}}, \sQuote{\code{-}}, \sQuote{\code{*}},
#' \sQuote{\code{/}}, \sQuote{\code{^}}), and those mathematical
#' functions defined in the C or the R math libraries (e.g.,
#' \code{fabs}, \code{exp}, \code{log}, \code{sin}).  (Notice that the
#' modulo operator \sQuote{\code{\%}} is currently not supported.)
#'
#' Identifiers in an \code{RxODE} model specification can refer to:
#' \itemize{
#'    \item state variables in the dynamic system (e.g., compartments in a
#'          pharmacokinetics/pharamcodynamics model);
#'    \item implied input variable, \code{t} (time), 
#'    \code{podo} (oral dose, for absorption models), and
#'    \code{tlast} (last time point);
#'    \item model parameters, (\code{ka} rate of absorption, \code{CL} 
#'        clearance, etc.);
#'    \item \code{pi}, for the constant pi. 
#'    \item others, as created by assignments as part of the model
#'          specification.
#' }
#'
#' Identifiers consists of case-sensitive alphanumeric characters,
#' plus the underscore \sQuote{_} character.  \strong{NB:} the dot
#' \sQuote{.} character is \strong{not} a valid character identifier.
#'
#' The values of these variables at pre-specified time points are
#' saved as part of the fitted/integrated/solved model (see
#' \code{\link{eventTable}}, in particular its member function
#' \code{add.sampling} that defines a set of time points at which to
#' capture a snapshot of the syste via the values of these variables).
#'
#' The ODE specification mini-language is parsed with the help of the
#' open source tool \emph{DParser}, Plevyak (2015).
#'
#' 
#' @return An object (closure) of class \dQuote{\code{RxODE}} (see Chambers and Temple Lang (2001))
#'      consisting of the following list of strings and functions:
#' 
#'     \item{modName}{the name of the model (a copy of the input argument).}
#'     \item{model}{a character string holding the source model specification.}
#'     \item{get.modelVars}{a function that returns a list with 3 character 
#'         vectors, \code{params}, \code{state}, and \code{lhs} of variable names used in the model
#'         specification. These will be output when the model is computed (i.e., the ODE solved by integration).}
#' 
#'       \item{solve}{this function solves (integrates) the ODE. This
#'           is done by passing the code to \code{\link{rxSolve}}.
#'           This is as if you called \code{rxSolve(RxODEobject,...)},
#'           but returns a matrix instead of a rxSolve object.
#'       
#'           \code{params}: a numeric named vector with values for every parameter 
#'           in the ODE system; the names must correspond to the parameter 
#'           identifiers used in the ODE specification;
#'           
#'           \code{events}: an \code{EventTable} object describing the 
#'           input (e.g., doses) to the dynamic system and observation 
#'           sampling time points (see  \code{\link{eventTable}});
#'       
#'           \code{inits}: a vector of initial values of the state variables
#'           (e.g., amounts in each compartment), and the order in this vector
#'           must be the same as the state variables (e.g., PK/PD compartments);
#'           
#'       
#'           \code{stiff}: a logical (\code{TRUE} by default) indicating whether
#'           the ODE system is stifff or not.  
#'       
#'           For stiff ODE sytems (\code{stiff=TRUE}), \code{RxODE} uses
#'           the LSODA (Livermore Solver for Ordinary Differential Equations) 
#'           Fortran package, which implements an automatic method switching 
#'           for stiff and non-stiff problems along the integration interval, 
#'           authored by Hindmarsh and Petzold (2003).  
#'       
#'           For non-stiff systems (\code{stiff=FALSE}), \code{RxODE} uses DOP853,
#'           an explicit Runge-Kutta method of order 8(5,3) of Dormand and Prince
#'           as implemented in C by Hairer and Wanner (1993).
#'       
#'           \code{trans_abs}: a logical (\code{FALSE} by default) indicating
#'           whether to fit a transit absorption term 
#'           (TODO: need further documentation and example);
#'       
#'           \code{atol}: a numeric absolute tolerance (1e-08 by default);
#'       
#'           \code{rtol}: a numeric relative tolerance (1e-06 by default).e
#'       
#'           The output of \dQuote{solve} is a matrix with as many rows as there
#'           are sampled time points and as many columns as system variables 
#'           (as defined by the ODEs and additional assigments in the RxODE model 
#'               code).}
#'       
#'       \item{isValid}{a function that (naively) checks for model validity,
#'           namely that the C object code reflects the latest model 
#'           specification.}
#'       \item{version}{a string scaler with the version of the \code{RxODE}
#'           object (not the package).}
#'       \item{dynLoad}{a function with one \code{force=FALSE} argument
#'           that dynamically loads the object code if needed.}
#'       \item{dynUnload}{a function with no argument that unloads 
#'           the model object code.}
#'       \item{cmpMgr}{a \dQuote{compilation manager} object, see 
#'           \code{\link{rx.initCmpMgr}}.}
#'       \item{delete}{removes all created model files, including C and DDL files.
#'           The model object is no longer valid and should be removed, e.g., 
#'           \code{rm(m1)}.}
#'       \item{run}{deprecated, use \code{solve}.}
#'       \item{parse}{deprecated.}
#'       \item{compile}{deprecated.}
#'       \item{get.index}{deprecated.}
#'       \item{getObj}{internal (not user callable) function.}
#'
#' @references
#'
#' Chamber, J. M. and Temple Lang, D. (2001)
#' \emph{Object Oriented Programming in R}. 
#' R News, Vol. 1, No. 3, September 2001.
#' \url{https://cran.r-project.org/doc/Rnews/Rnews_2001-3.pdf}.
#'
#' Hindmarsh, A. C.
#' \emph{ODEPACK, A Systematized Collection of ODE Solvers}.
#' Scientific Computing, R. S. Stepleman et al. (Eds.),
#' North-Holland, Amsterdam, 1983, pp. 55-64.
#'
#' Petzold, L. R.
#' \emph{Automatic Selection of Methods for Solving Stiff and Nonstiff
#' Systems of Ordinary Differential Equations}.
#' Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
#'
#' Hairer, E., Norsett, S. P., and Wanner, G.
#' \emph{Solving ordinary differential equations I, nonstiff problems}.
#' 2nd edition, Springer Series in Computational Mathematics,
#' Springer-Verlag (1993).
#'
#' Plevyek, J.
#' \emph{Dparser}, \url{http://dparser.sourceforge.net}. Web. 12 Oct. 2015.
#'
#' @author Melissa Hallow, Wenping Wang, David James and Matthew Fidler
#'
#' @seealso \code{\link{eventTable}}
#'
#' @examples
#' # Step 1 - Create a model specification
#' ode <- "
#'    # A 4-compartment model, 3 PK and a PD (effect) compartment
#'    # (notice state variable names 'depot', 'centr', 'peri', 'eff')
#' 
#'    C2 = centr/V2;
#'    C3 = peri/V3;
#'    d/dt(depot) =-KA*depot;
#'    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
#'    d/dt(peri)  =                    Q*C2 - Q*C3;
#'    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
#' "
#' 
#' m1 <- RxODE(model = ode, modName = "m1")
#' print(m1)
#' 
#' # Step 2 - Create the model input as an EventTable,
#' # including dosing and observation (sampling) events
#' 
#' # QD (once daily) dosing for 5 days.
#' 
#' qd <- eventTable(amount.units="ug", time.units = "hours")
#' qd$add.dosing(dose=10000, nbr.doses=5, dosing.interval = 24)
#' 
#' # Sample the system hourly during the first day, every 8 hours
#' # then after
#' 
#' qd$add.sampling(0:24)
#' qd$add.sampling(seq(from = 24+8, to = 5*24, by = 8))
#' 
#' # Step 3 - set starting parameter estimates and initial
#' # values of the state
#' 
#' theta <- 
#'     c(KA=.291, CL=18.6, 
#'       V2=40.2, Q=10.5, V3=297.0,
#'       Kin=1.0, Kout=1.0, EC50=200.0)
#' 
#' # init state variable
#' inits <- c(0, 0, 0, 1)      
#' 
#' # Step 4 - Fit the model to the data
#' 
#' qd.cp <- m1$solve(theta, events = qd, inits)
#' 
#' head(qd.cp)
#' 
#' @keywords models nonlinear
#' @concept Nonlinear regression
#' @concept ODE models
#' @concept Ordinary differential equations
#' @concept Pharmacokinetics (PK)
#' @concept Pharmacodynamics (PD)
#' @useDynLib RxODE trans
#' @export
RxODE <-
    function(model, modName = basename(wd), wd = getwd(),
             filename = NULL, do.compile = NULL, extraC=NULL,
             debug = FALSE,...)
{
    if(!missing(model) && !missing(filename))
        stop("must specify exactly one of 'model' or 'filename'")
    if (missing(model) && !missing(filename))
        model <- filename;
    ## RxODE compilation manager (location of parsed code, generated C,  shared libs, etc.)

    cmpMgr <- rx.initCmpMgr(model, modName, wd,  extraC, debug, missing(modName));
    ## NB: the set of model variables (modelVars) is only available 
    ## after parsing, thus it needs to be dynamically computed in cmpMgr
    get.modelVars <- cmpMgr$get.modelVars

    .version <- "0.6"          # object version
    .last.solve.args <- NULL   # to be populated by solve()

    solve <- function(...){
        .last.solve.args <<- as.list(match.call(expand.dots = TRUE));
        rx <- cmpMgr$rxDll();
        return(as.matrix(rxSolve(rx,...)));
    }
    force <- FALSE
    if (class(do.compile) == "logical"){
        if (do.compile)
            force <- TRUE;
    }
    if (is.null(do.compile)){
        do.compile <- TRUE
    }
    if (do.compile){
        cmpMgr$compile(force);
        ## Add for backward compatibility.
        cmpMgr$dllfile <- rxDll(cmpMgr$rxDll());
        cmpMgr$ode_solver <- rxTrans(cmpMgr$rxDll())["ode_solver"];
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
            getObj = function(obj) get(obj, envir = environment(solve)))
    class(out) <- "RxODE"
   out
}

RxODE.quoteVar <- function(var,word=TRUE){
    ret <- c();
    for (i in 1:nchar(var)){
        ret <- c(ret,sprintf("one_of(\"%s\")",substr(var,i,i)));
    }
    ret <- parse(text=sprintf(ifelse(word,"rex::rex(boundary,%s,boundary)",
                                     "rex::rex(%s)"),
                              paste(ret,collapse=",")));
    ret <- eval(ret)
    return(ret);
}

#' Predict an RxODE object
#'
#' \code{predict} solves the odinary differential equations specified by
#'  RxODE object
#'
#' @param object An RxODE object
#' @param ... Solve arguments sent to \code{rxSolve}.  See
#'     \code{\link{rxSolve}}.
#' 
#' @export 
predict.RxODE <- function(object,...){
    rxSolve(object,...);
}

#' Solve RxODE objects
#'
#' @param ... Additional arguments sent to \code{rxSolve}
#'
#' @seealso \code{\link{rxSolve}}
#' @export
solve.RxODE <- function(...){
    rxSolve(...)
}
#' @rdname solve.RxODE
#' @export 
solve.RxCompilationManager <- function(...){
    rxSolve(...)
}
#' @rdname solve.RxODE
#' @export 
solve.solveRxDll <- function(...){
    rxSolve(...)
}
#' @rdname solve.RxODE
#' @export 
solve.character <- function(...){
    rxSolve(...)
}

#' @rdname solve.RxODE
#' @export
solve.rxDll <- function(...){
    rxSolve(...);
}

#' Add item to solved system of equations
#'
#' @param obj1 Solved object
#'
#' @param obj2 New object to be added/piped/chained to solved object.
#'
#' When \code{newObject} is an event table, return a new solved object
#' with the new event table.
#' @name rxChain

#' @rdname rxChain
#' @export
'+.solveRxDll' <- function(obj1, obj2) {
    args <- rev(as.list(match.call())[-1]);
    names(args) <- c("obj","solvedObject");
    return(do.call("rxChain2",args,envir = parent.frame(1)));
}

#' @rdname rxChain
#' @export
'%>%.solveRxDll' <- '+.solveRxDll';

#' Second command in chaining commands
#'
#' This is s3 method is called internally with \code{+} and \code{\%>\%} operators.
#'
#' @param obj the object being added/chained/piped to the solved object
#' @param solvedObject the solved object
#' @param envir
#' 
#' @keywords internal
#' @export
rxChain2 <- function(obj,solvedObject){
    UseMethod("rxChain2")
}

#' @rdname rxChain2
#' @export
rxChain2.default <- function(obj,solvedObject){
    args <- as.list(match.call());
    stop(sprintf("Do not know how to add %s to RxODE solved object %s",toString(args[[2]]), toString(args[[3]])))
}

#' @rdname rxChain2
#' @export
rxChain2.EventTable <- function(obj,solvedObject){
    args <- rev(as.list(match.call())[-1]);
    names(args) <- c("object","events");
    return(do.call("rxSolve",args, envir = parent.frame(1)));
}
 
#' Print information about the RxODE object.
#'
#' This prints the model name and its status for being able to be solved
#'
#' @param x An rxode object
#' @param ... Ignored parameters
#' @export print.RxODE
print.RxODE <-
    function(x, ...)
{
    valid <- x$cmpMgr$isValid()
    if(!valid){
        .msg <- "invalid object, needs to be re-created"
    } else {
        .ready <- x$cmpMgr$getObj(".compiled")
        .msg <- if(.ready) "ready to run" else "needs compilation"
    }
    cat(sprintf('RxODE model named "%s" (%s)\n', x$modName, .msg))
    invisible(x)
}

#' Print expanded information about the RxODE object.
#'
#' This prints the expanded information about the RxODE object.
#'
#' @param object RxODE object
#' @param ... Ignored parameters
#' @export summary.RxODE
summary.RxODE <- function(object, ...)
{
    print.RxODE(object);
    summary.rxDll(object$cmpMgr$rxDll(),noprint=TRUE)
    invisible(object)
}

#' @rdname summary.RxODE
#' @export
summary.RxCompilationManager <- function(object, ...)
{
    print.RxCompilationManager(object);
    summary.rxDll(object$rxDll(),noprint=TRUE)
    invisible(object);
}



#' Return the RxODE coefficients
#'
#' This returns the parameters , state variables
#'
#' @param object is an RxODE object
#' @param ... ignored arguments
#'
#' @return a rxCoef object with the following
#'
#' \item{params}{ is a list of strings for parameters for the RxODE object}
#' \item{state}{ is a list of strings for the names of each state in
#'     the RxODE object.}
#' \item{ini}{ is the model specified default values for the
#'     parameters.}
#' \item{RxODE}{ is the referring RxODE object}
#' @export 
coef.RxODE <- function(object,
                       ...){
    ret <- rxModelVars(object)[c("params","state","ini")];
    ret$RxODE <- object;
    class(ret) <- "rxCoef";
    return(ret);
}

#' @rdname coef.RxODE
#' @export 
coef.RxCompilationManager <- function(...){
    coef.RxODE(...);
}

#' @rdname coef.RxODE
#' @export 
coef.solveRxDll <- function(...){
    coef.RxODE(...);
}

#' @rdname coef.RxODE
#' @export 
coef.rxDll <- function(...){
    coef.RxODE(...);
};

#' Print the rxCoef object
#'
#' This prints out the user supplied arguments for rxCoef object
#'
#' @param x rxCoef object
#' 
#' @keywords internal
#' @export print.rxCoef
print.rxCoef <- function(x,...){
    rxDllObj <- x$RxODE;
    if (length(rxParams(rxDllObj)) > 0){
        cat("\nUser Supplied Parameters:\n");
        print(rxInits(rxDllObj,c(),rxParams(rxDllObj),NA,TRUE))
        cat("\nUser Initial Conditions:\n");
        print(rxInits(rxDllObj,c(),rxState(rxDllObj),0,TRUE))
    }
    if (length(rxInits(rxDllObj)) > 0){
        cat("\nDefault parameter values:\n")
        print(x$ini);
    }
    cat("\nCompartents:\n");
    tmp <- rxState(rxDllObj);
    names(tmp) <- paste0("cmt=",1:length(tmp));
    print(tmp);
    return(invisible());
}

igraph <- function(obj,...){
    ## In case anyone else wants to use the method...
    UseMethod("igraph");
}

rxSub <- function(regexp, # regular expression with named groups
                  what,   # replacement with \g{named_group} replacement
                  text,   # Text to replace
                  ...){
    ## rxSub returns a replaced character string, like gsub with named groups
    reg <- regexpr(regexp,text,perl=TRUE);
    if (reg != -1){
        n <- apply(rbind(attr(reg,"capture.start"),attr(reg,"capture.length")),2,function(x){
            return(substr(text,x[1],x[1]+x[2]-1))
        })
        ## Support \g{} backreference in replacement
        names(n) <- paste0("\\g{",attr(reg,"capture.names"),"}");
        m <- what;
        for (i in 1:length(n)){
            m <- gsub(names(n)[i],n[i],m,fixed=TRUE);
        }
        text <- paste0(substr(text,0,reg-1),m,
                      substr(text,reg+attr(reg,"match.length"),
                             nchar(text)));
        return(text);
    } else{
        return(text);
    }
} # end function rxSub

idrInfo <- function(cmt="eff",text="+Kin-Kout*(1-centr/V2/(EC50+centr/V2))*eff"){
    text <- sigmoidInfo(text);
    id <- rex::rex(one_of("_","a":"z","A":"Z"),any_of("_","a":"z","A":"Z","0":"9"));
    resp <- rex::rex("Hill(",capture(except_some_of(",",newline),name="emax"),## Emax
                     ",",capture(except_some_of(",",newline),name="cp"),## Cp
                     ",",capture(except_some_of(",",newline),name="e50"),## e50
                     ",",capture(except_some_of(",",newline),name="gamma"),## gamma
                     ")")
    kout <- rex::rex(capture(id,name="kout"),any_spaces);
    kin <- rex::rex(capture(id,name="kin"),any_spaces);
    R <- rex::rex(capture(cmt,name="R"),any_spaces);
    inb <- rex::rex("(",any_spaces,or("1","1.0"),any_spaces,"-",any_spaces,resp,any_spaces,")",any_spaces);
    stim <- rex::rex("(",any_spaces,or("1","1.0"),any_spaces,"+",any_spaces,resp,any_spaces,")",any_spaces);
    stim2 <- rex::rex("(",any_spaces,resp,any_spaces,"+",any_spaces,or("1","1.0"),any_spaces,")",any_spaces);
    regs <- list(
        ## IDR1 kin*(1-sig)-kout*R
        c(rex::rex(kin,"*",any_spaces,inb,"-",any_spaces,R,"*",any_spaces,kout),
          "IDR1"),
        c(rex::rex(kin,"*",any_spaces,inb,"-",any_spaces,kout,"*",any_spaces,R),
          "IDR1"),
        c(rex::rex(inb,"*",any_spaces,kin,"-",any_spaces,R,"*",any_spaces,kout),
          "IDR1"),
        c(rex::rex(inb,"*",any_spaces,kin,"-",any_spaces,kout,"*",any_spaces,R),
          "IDR1"),
        ## IDR2
        c(rex::rex(kin,"-",any_spaces,kout,"*",any_spaces,inb,"*",any_spaces,R),
          "IDR2"),
        c(rex::rex(kin,"-",any_spaces,kout,"*",any_spaces,R,"*",any_spaces,inb),
          "IDR2"),
        c(rex::rex(kin,"-",any_spaces,inb,"*",any_spaces,kout,"*",any_spaces,R),
          "IDR2"),
        c(rex::rex(kin,"-",any_spaces,inb,"*",any_spaces,R,"*",any_spaces,kout),
          "IDR2"),
        c(rex::rex(kin,"-",any_spaces,R,"*",any_spaces,kout,"*",any_spaces,inb),
          "IDR2"),
        c(rex::rex(kin,"-",any_spaces,R,"*",any_spaces,inb,"*",any_spaces,kout),
          "IDR2"),
        ## IDR3
        c(rex::rex(kin,"*",any_spaces,stim,"-",any_spaces,kout,"*",any_spaces,R),
          "IDR3"),
        c(rex::rex(stim,"*",any_spaces,kin,"-",any_spaces,kout,"*",any_spaces,R),
          "IDR3"),
        c(rex::rex(kin,"*",any_spaces,stim,"-",any_spaces,R,"*",any_spaces,kout),
          "IDR3"),
        c(rex::rex(stim,"*",any_spaces,kin,"-",any_spaces,R,"*",any_spaces,kout),
          "IDR3"),
        c(rex::rex(kin,"*",any_spaces,stim2,"-",any_spaces,kout,"*",any_spaces,R),
          "IDR3"),
        c(rex::rex(stim2,"*",any_spaces,kin,"-",any_spaces,kout,"*",any_spaces,R),
          "IDR3"),
        c(rex::rex(kin,"*",any_spaces,stim2,"-",any_spaces,R,"*",any_spaces,kout),
          "IDR3"),
        c(rex::rex(stim2,"*",any_spaces,kin,"-",any_spaces,R,"*",any_spaces,kout),
          "IDR3"),
        ## IDR4
        c(rex::rex(kin,"-",any_spaces,kout,"*",any_spaces,stim,"*",any_spaces,R),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,kout,"*",any_spaces,R,"*",any_spaces,stim),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,stim,"*",any_spaces,kout,"*",any_spaces,R),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,stim,"*",any_spaces,R,"*",any_spaces,kout),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,R,"*",any_spaces,stim,"*",any_spaces,kout),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,R,"*",any_spaces,kout,"*",any_spaces,stim),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,kout,"*",any_spaces,stim2,"*",any_spaces,R),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,kout,"*",any_spaces,R,"*",any_spaces,stim2),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,stim2,"*",any_spaces,kout,"*",any_spaces,R),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,stim2,"*",any_spaces,R,"*",any_spaces,kout),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,R,"*",any_spaces,stim2,"*",any_spaces,kout),
          "IDR4"),
        c(rex::rex(kin,"-",any_spaces,R,"*",any_spaces,kout,"*",any_spaces,stim2),
          "IDR4")
        );
    ret <- text
    for (i in 1:length(regs)){
        ret <- rxSub(regs[[i]][1],sprintf("%s(\"\\g{kin}\",\"\\g{kout}\",\"\\g{e50}\",\"\\g{R}\",\"\\g{cp}\",edgeList,nodes)",regs[[i]][2]),ret);
    }
    ret <- gsub(rex::rex(start,"+"),"",ret);
    return(ret);
}

sigmoidInfo <- function(text = "+Kin-Kout*(1-centr/V2/(EC50+centr/V2))*eff"){
    ## Emax models
    id <- rex::rex(one_of("_","a":"z","A":"Z"),any_of("_","a":"z","A":"Z","0":"9"));

    cp <- rex::rex(or(id, ## Cp
                      group(id,any_spaces,"/",any_spaces,id), ## AMT/V
                      group("(",any_spaces,id,any_spaces,"/",any_spaces,id,any_spaces,")") ## (AMT/V)
                      ))
    regs <- list(
        ## Emax*C/(C+E50)
        c(rex::rex(capture(id,name="emax"),any_spaces,"*",any_spaces,capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")"),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## Emax*C/(E50+C)
        c(rex::rex(capture(id,name="emax"),any_spaces,"*",any_spaces,capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")"),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),

        ## C*Emax/(C+E50)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",any_spaces,capture(id,name="emax"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")"),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C*Emax/(E50+C)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",any_spaces,capture(id,name="emax"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")"),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C/(C+E50)*Emax
        c(rex::rex(capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C/(E50+C)*Emax
        c(rex::rex(capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),

        ## (C/(C+E50))*Emax
        c(rex::rex("(",any_spaces,capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,")",any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        
        ## (C/(E50+C))*Emax
        c(rex::rex("(",any_spaces,capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,")",any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C/(C+E50)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"/",
                            any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")"),"Hill(1,\\g{cp},\\g{e50},1)"),
        ## C/(E50+C)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"/",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")"),"Hill(1,\\g{cp},\\g{e50},1)"),
        ## C*(E50+c)^-1 variants
        ## Emax*C*(C+E50)
        c(rex::rex(capture(id,name="emax"),any_spaces,"*",any_spaces,capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## Emax*C*(E50+C)
        c(rex::rex(capture(id,name="emax"),any_spaces,"*",any_spaces,capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),

        ## C*Emax*(C+E50)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",any_spaces,capture(id,name="emax"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C*Emax*(E50+C)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",any_spaces,capture(id,name="emax"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C*(C+E50)*Emax
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C*(E50+C)*Emax
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)"),any_spaces,"*",
                   any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),

        ## (C*(C+E50))*Emax
        c(rex::rex("(",any_spaces,capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)"),
                   any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        
        ## (C*(E50+C))*Emax
        c(rex::rex("(",any_spaces,capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)"),
                   any_spaces,"*",any_spaces,capture(id,name="emax")),"Hill(\\g{emax},\\g{cp},\\g{e50},1)"),
        ## C*(C+E50)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture_group("cp"),any_spaces,"+",any_spaces,capture(id,name="e50"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)")),"Hill(1,\\g{cp},\\g{e50},1)"),
        ## C*(E50+C)
        c(rex::rex(capture(cp,name="cp"),any_spaces,"*",
                   any_spaces,"(",any_spaces, capture(id,name="e50"),any_spaces,"+",any_spaces,capture_group("cp"),
                   any_spaces,")",any_spaces,or("**","^"),any_spaces,or("-1","(-1)","-1.0","(-1.0)")),"Hill(1,\\g{cp},\\g{e50},1)")
        ## FIXME: Emax/(1+EC50/Cp) variants
        ## FIXME: -- log linerized forms?
        ## FIXME: Hill Variants
        ## FIXME: Hodgkin?
        ## FIXME: Douglas?
        ## FIXME: Gompertz?
        );
    ret <- text
    for (i in 1:length(regs)){
        ret <- rxSub(regs[[i]][1],regs[[i]][2],ret);
    }
    return(ret);
}

nodeInfo <- function(x,       # RxODE normalized model
                     modVars, # model Vars
                     ...){
    ## nodeInfo returns a list containing edgeList, biList and nodes
    
    ##  The nodes object is the compartments in the model that will be
    ##  drawn.  Compartments that should be hidden are prefixed with a
    ##  "."

    ## The edgeList has a list of nodes and the directions they
    ## connect like c("A","B") is an arrow from compartent A to B.

    ## The names of the edgeList is the label applied to the arrow.
    ##
    
    ## Move any of the equals logical operators to ~ operators so they
    ## wont be split...
    mod0 <- gsub(rex::rex(capture(one_of(">!<")),"="),
                 "\\1~",x,perl=TRUE);
    mod0 <- gsub(rex::rex("=="),"~~",mod0,perl=TRUE);
    ## Make sure there are no spaces around remaining = sign
    mod0 <- gsub(rex::rex(any_spaces,"=",any_spaces),"=",mod0,perl=TRUE)
    ## remove the newlines
    mod0 <- gsub(rex::rex(newlines),"",mod0,perl=TRUE)
    ## Remove any if statements
    mod0 <- gsub(rex::rex(any_spaces,"if",any_spaces,
                          one_of("("),except_any_of(")"),one_of(")"),
                          any_spaces,one_of("{"),
                          except_any_of("}"),one_of("}")),
                 "",mod0,perl=TRUE);
    ## Remove any else statements
    mod0 <- gsub(rex::rex(any_spaces,"else",any_spaces,one_of("{"),except_any_of("}"),one_of("}")),
                 "",mod0,perl=TRUE);
    ## Split by semicolon, and =.
    mod0 <- strsplit(strsplit(mod0,rex::rex(some_of(";")))[[1]],
                     rex::rex(one_of("=")));
    mod <- eval(parse(text=sprintf("c(%s)",paste(unlist(lapply(mod0,function(x){
        if (length(x) == 2){
            x1 <- gsub(rex::rex(start,any_spaces,capture(anything),any_spaces,end),
                       "\\1",x[1]);
            x2 <- gsub(rex::rex(start,any_spaces,capture(anything),any_spaces,end),
                       "\\1",x[2]);
            return(sprintf("\"%s\"=\"%s\"",x1,x2));
        } else {
            return(NULL);
        }
    })),collapse=","))))
    
    ## Expand  - abc*() and - ()*abc expressions
    expandFactor <- function(var = "F*KA*depot-C2*(CL+Q)+Q*C3" ){
        ## Expand any simple factored expressions.
        pr <- rex::rex(anything,
                       capture(one_of("-+"),any_spaces,
                               or(
                                   ## - abc * (a+b+c+d+e...)
                                   group(except_some_of("*"),
                                         any_spaces,
                                         one_of("*"),
                                         any_spaces,
                                         one_of("("),
                                         except_any_of("()"),
                                         one_of(")")),
                                   ## - (a+b+c+d+e) * abc
                                   group(one_of("("),
                                         except_any_of("()"),
                                         one_of(")"),
                                         any_spaces,
                                         one_of("*"),
                                         or(except_some_of("*-+;"),
                                            newline)))),
                       anything);
        neg <- rex::rex(start,any_spaces,one_of("-+"))
        ## Add initial + if needed
        if (regexpr(neg,var,perl=TRUE) == -1){
            var <- sprintf("+%s",var);
        }
        while(regexpr(pr,var,perl=TRUE) != -1){
            ## Get expression to expand
            totalExpr <- gsub(pr,"\\1",var,perl=TRUE);
            v <- gsub(rex::rex(any_of("-+"),any_spaces),"",
                      gsub(rex::rex(any_of("*"),
                                    any_spaces,
                                    one_of("("),except_any_of("()"),one_of(")"),
                                    any_spaces,
                                    any_of("*")),
                           "",totalExpr,perl=TRUE),perl=TRUE);
            ## Is this an negative expression?
            negVar <- FALSE
            if (regexpr(neg,totalExpr,perl=TRUE) != -1){
                negVar <- TRUE;
            }
            ## Capture the parenthetical expression
            p <- gsub(rex::rex(anything,
                               one_of("("),capture(except_any_of("()")),one_of(")"),
                               anything),"\\1",totalExpr,perl=TRUE)
            ## Add a positive to the beginning of the expression, if needed.
            if (regexpr(rex::rex(start,any_spaces,one_of("-")),
                       p,perl=TRUE) == -1){
                p <- paste0("+",p);
            }
            ## Protect negative expressions for later use
            p <- gsub(rex::rex(one_of("-")),"~~~~",p,perl=TRUE);
            ## Expand the positive expression
            p <- gsub(rex::rex(one_of("+")),sprintf("%s%s*",ifelse(negVar,"-","+"),v),p,perl=TRUE)
            ## Expand the negative expression
            p <- gsub(rex::rex(n_times(one_of("~"),4)),sprintf("%s%s*",ifelse(negVar,"+","-"),v),p,perl=TRUE)
            ## replace in original
            var <- gsub(RxODE.quoteVar(totalExpr,FALSE),p,var,perl=TRUE);
            print(var);
        }
        ## Replace initial +
        var <- gsub(rex::rex(start,any_spaces,one_of("+"),any_spaces),"",var);
        return(var)
    }
    ## Expand factored expressions
    for (i in 1:length(mod)){
        mod[i] <- expandFactor(mod[i]);
    }
    ## now replace defined variables in overall expressions.
    for (i in 1:length(mod)){
        if (i > 1){
            for (j in seq(1,i-1)){
                mod[i] <- gsub(RxODE.quoteVar(names(mod)[j]),
                               mod[j],mod[i],perl=TRUE);
            }
        }
    }
    ## Get the variables from the parser
    
    ## Subset to the ODE equations.
    de <- c()
    for (v in modVars$state){
        de <- c(de,mod[sprintf("d/dt(%s)",v)])
    }
    w <- which(regexpr(rex::rex(start,any_spaces,one_of("-")),de,perl=TRUE) == -1);
    de[w] <- paste0("+",de[w]);
    names(de) <- modVars$state;

    fullDe <- de;
    
    ## Currently only parses simple expressions.
    ## FIXME: Add michelis-menton / Hill kinetics
    de <- de[which(regexpr("[(]",de) == -1)];
    parDe <- de;

    sortDe <- function(x =c("centr/V2*CL","centr/V2*Q")){
        ## sortDe returns
        ## First protect the / operators.
        x <- unlist(lapply(strsplit(gsub(rex::rex(any_spaces,one_of("/"),any_spaces),
                                  "/zzzzzzz",x[x != ""],perl=TRUE),
                             rex::rex(any_spaces,one_of("/*"),any_spaces),
                             perl=TRUE),
                    function(x){
                        x <- paste(gsub(rex::rex(one_of("*"),n_times("z",7)),"/",paste0("*",sort(x))),collapse="");
                        return(x)
                    }));
        return(x);
    } # end function sortDe
    
    ## Find the negative exressions in the ODEs
    negDe <- gsub(rex::rex(one_of("+"),any_spaces,except_any_of("+-")),"",de,perl=TRUE);
    negDe <- negDe[negDe != ""];
    negDe <- lapply(strsplit(negDe,rex::rex(any_spaces,one_of("-"),any_spaces)),sortDe);
    
    
    ## Find the positive expressions in the ODEs
    de <- gsub(rex::rex(one_of("-"),any_spaces,except_any_of("+-")),"",de,perl=TRUE) 
    nde <- names(de)
    de <- de[de != ""];
    de <- lapply(strsplit(de,rex::rex(any_spaces,one_of("+"),any_spaces)),sortDe);
    for (v in names(negDe)){
        negDe[[v]] <- gsub(rex::rex(any_spaces,one_of("*"),any_spaces,RxODE.quoteVar(v)),"",negDe[[v]])
    }
    ## Look for connections between compartments from the positive direction
    nodes <- c();
    edges <- NULL;
    edgeList <- list();
    biList <- list();
    for (n in names(de)){
        for (n2 in nde){
            reg <- rex::rex(any_spaces,one_of("*"),any_spaces,RxODE.quoteVar(n2))
            w <- which(regexpr(reg,de[[n]],perl=TRUE) != -1);
            if (length(w) == 1){
                nodes <- c(nodes,n2,n);
                var <- gsub(reg,"",de[[n]][w],perl=TRUE);
                for (e in edgeList){
                    if (e[1] == n & e[2] == n2){
                        biList[[length(biList) + 1]] <- c(n,n2);
                    }
                }
                edgeList[[var]] <- c(n2,n);
                ## If the corresponding negative direction exists,
                ## remove it from negDe.  The negDe should have overall body clearances.
                noFvar <- gsub(rex::rex(any_spaces,one_of("*"),any_spaces,one_of("Ff"),any_alnums),"",var,perl=TRUE)
                if (any(names(negDe) == n2)){
                    negDe[[n2]] <- negDe[[n2]][negDe[[n2]] != var];
                    negDe[[n2]] <- negDe[[n2]][negDe[[n2]] != noFvar];

                }
            }
        }
    }
    ## Assign overall body clerances
    nodes <- unique(nodes);
    outn <- 1;
    if (length(negDe) > 0){
        for (i in 1:length(negDe)){
            if (length(negDe[[i]]) == 1){
                outNode <- sprintf(".out%s",outn)
                edgeList[[negDe[[i]]]] <- c(names(negDe)[i],outNode);
                nodes <- c(nodes,outNode);
                outn <- outn +1;
            }
        }
    }
    reg <- rex::rex(capture(anything),any_spaces,one_of("*"),any_spaces,capture(one_of("Ff"),any_alnums),capture(anything));
    w <- which(regexpr(reg,names(edgeList),perl=TRUE) != -1);
    if (length(w) > 0){
        names(edgeList)[w] <- gsub(reg,"\\1\\3\n\\2",names(edgeList)[w],perl=TRUE);
    }
    names(edgeList) <- gsub(rex::rex(any_spaces,one_of("*"),any_spaces),"",names(edgeList),perl=TRUE)
    ## Now see if we can recognize kin/kout from indirect response models
    idr <- sapply(1:length(fullDe),function(x){
        tmp <- idrInfo(names(fullDe)[x],fullDe[x]);
        if (regexpr(rex::rex(start,"IDR"),tmp) != -1){
            return(tmp)
        } else {
            return("")
        }
    })
    names(idr) <- names(fullDe);
    idr <- idr[idr != ""];
    w <- which(names(fullDe) %in% names(idr));        
    fReg <- sprintf("^.*(%s).*$",paste(names(fullDe)[-w],collapse="|"));
    IDR1 <- function(kin,kout,e50,ef,cp,edgeList,nodes){
        kout <- sprintf("%s\n",kout);
        kin  <- sprintf("%s\n\u25BC %s",kin,e50);
        idr <- "Indirect\nEffect (I)";
        idrFrom <- gsub(fReg,"\\1",cp);
        edgeList[[idr]] <- c(idrFrom,ef);
        edgeList[[kin]] <- c(".Kin",ef);
        edgeList[[kout]] <- c(ef,".Kout");
        nodes <- c(nodes,ef,".Kin",".Kout");
        return(list(edgeList,nodes))
    }
    IDR2 <- function(kin,kout,e50,ef,cp,edgeList,nodes){
        kout <- sprintf("%s\n\u25BC %s",kout,e50);
        kin  <- sprintf("%s\n",kin);
        idr <- "Indirect\nEffect (II)";
        idrFrom <- gsub(fReg,"\\1",cp);
        edgeList[[idr]] <- c(idrFrom,ef);
        edgeList[[kin]] <- c(".Kin",ef);
        edgeList[[kout]] <- c(ef,".Kout");
        nodes <- c(nodes,ef,".Kin",".Kout");
        return(list(edgeList,nodes))
    }
    IDR3 <- function(kin,kout,e50,ef,cp,edgeList,nodes){
        kout <- sprintf("%s\n",kout);
        kin  <- sprintf("%s\n\u25B3 %s",kin,e50);
        idr <- "Indirect\nEffect (III)";
        idrFrom <- gsub(fReg,"\\1",cp);
        edgeList[[idr]] <- c(idrFrom,ef);
        edgeList[[kin]] <- c(".Kin",ef);
        edgeList[[kout]] <- c(ef,".Kout");
        nodes <- c(nodes,ef,".Kin",".Kout");
        return(list(edgeList,nodes))
    }
    IDR4 <- function(kin,kout,e50,ef,cp,edgeList,nodes){
        kout <- sprintf("%s\n\u25B3 %s",kout,e50);
        kin  <- sprintf("%s\n",kin);
        idr <- "Indirect\nEffect (IV)";
        idrFrom <- gsub(fReg,"\\1",cp);
        edgeList[[idr]] <- c(idrFrom,ef);
        edgeList[[kin]] <- c(".Kin",ef);
        edgeList[[kout]] <- c(ef,".Kout");
        nodes <- c(nodes,ef,".Kin",".Kout");
        return(list(edgeList,nodes))
    }    
    if (length(idr) > 0){
        tmp <- eval(parse(text=idr));
        edgeList <- tmp[[1]];
        nodes <- tmp[[2]];
        print(edgeList);
        print(nodes);
    }
    return(list(nodes    = nodes,
                edgeList = edgeList,
                biList   = biList));
} # end function nodeInfo

igraph.rxDll <- function(x,                                   #  object
                         shape      = c("square","circle","csquare","rectangle","crectangle","vrectangle","sphere","none"),
                         size       = 30,                     # Size of square
                         colors     = c("accent1"="#0460A9"), # Colors
                         fillColor  = "accent1",              # Color to fill 
                         family     = "sans",                 # font Family
                         font       = 2,                      # Font (1: plain, 2: bold, 3: italic, 4: bold/italic, 5:symbol)
                         labelColor = "white",                # Label Color
                         lineColor  = "accent1",              # Line Color
                         shapeEnd   = c("none","sphere","circle","square","csquare","rectangle","crectangle","vrectangle"),
                         sizeEnd    = 10,                     # Size of End
                         tk         = FALSE,                  # For tkplot; transparancy not supported.
                         ...){
    ## igraph.rxDll returns igraph object from rxDll
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package igraph needed for this function to work. Please install it.",
             call. = FALSE)
    }
    with(nodeInfo(rxModelVars(x)$model["normModel"],rxModelVars(x)),{
        ret <- eval(parse(text=sprintf("igraph::graph_from_literal(%s);",paste(unlist(lapply(edgeList,function(x){sprintf("\"%s\" -+ \"%s\"",x[1],x[2])})),collapse=","))));
        if (length(shape) > 1){
            shape <- shape[1];
        }
        if (length(shapeEnd) > 1){
            shapeEnd <- shapeEnd[1];
        }
        getColor <- function(x){
            if (any(names(colors) == x)){
                return(colors[x]);
            } else {
                return(x);
            }
        }
        ret <- igraph::set_vertex_attr(ret,"shape",value=shape);
        ret <- igraph::set_vertex_attr(ret,"size",value=size);
        ret <- igraph::set_vertex_attr(ret,"color",value=getColor(fillColor));
        ret <- igraph::set_vertex_attr(ret,"label.family",value=family)
        ret <- igraph::set_vertex_attr(ret,"label.font",value=font);
        ret <- igraph::set_vertex_attr(ret,"label.color",value=getColor(labelColor));
        for (n in nodes){
            if (substring(n,0,1) == "."){
                if (!tk){
                    ret <- igraph::set_vertex_attr(ret,"label.color",igraph::V(ret)[n],value="transparent");
                    ret <- igraph::set_vertex_attr(ret,"shape",igraph::V(ret)[n],value=shapeEnd);
                    ret <- igraph::set_vertex_attr(ret,"size",igraph::V(ret)[n],value=sizeEnd);
                } else{
                    ret <- igraph::set_vertex_attr(ret,"size",igraph::V(ret)[n],value=sizeEnd);
                }
            }
        }
        ret <- igraph::set_edge_attr(ret,"color",value=getColor(lineColor))
        ret <- igraph::set_edge_attr(ret,"curved",value=FALSE)
        ret <- igraph::set_edge_attr(ret,"label.font",value=font)
        ## Set bidirectional to curved.
        eval(parse(text=paste(lapply(biList,function(x){sprintf("ret <- igraph::set_edge_attr(ret,\"curved\",igraph::E(ret)[\"%s\" %%--%% \"%s\"],value=TRUE)",x[1],x[2])}),collapse=";")));
        ## add the labels
        for (l in names(edgeList)){
            v <- edgeList[[l]];
            eg <- igraph::E(ret)[v[1] %->% v[2]];
            ret <- igraph::set_edge_attr(ret,"label",eg,value=l);
        }
        return(ret);
    });
} # end function igraph.rxDll

#' Plot a model digram for simple ODE models
#'
#' This plots a line diagram for the current system of differential
#' equations.
#'
#' @param RxODEobj RxODE object
#' @param family Font family to use
#' @param interactive If true, use Tcl/tk to plot the diagram and
#'     allow movement of the compartments
#' @param shape Box shape, by default it is a square
#' @param size size of square for compartments
#' @param colors is a named list of colors.  If any of the names match
#'     the colors specified by this function, it will use the color in
#'     this list.  For example \code{accent1} would match the default
#'     color in this list.
#' @param fillColor Color of the fill
#' @param font Font (1: plain, 2: bold, 3: italic, 4: bold/italic,
#'     5:symbol)
#' @param labelColor Color of the label
#' @param lineColor Color of the line
#' @param shapeEnd End shape
#' @param sizeEnd Size of end
#' @param ... Other arguments passed to igraph's plot method.
#'
#' Currently this method supports simple elimination type models and
#' indirect-effect models.  Michelis Menton models are not yet
#' supported.
#'
#' @export
rxPlot <- function(RxODEobj,
                   family      ="sans",
                   interactive = FALSE,
                   shape       = c("square","circle","csquare","rectangle","crectangle","vrectangle","sphere","none"),
                   size        = 30,                     # Size of square
                   colors      = c("accent1"="#0460A9"), # Colors
                   fillColor   = "accent1",              # Color to fill 
                   font        = 2,                      # Font (1: plain, 2: bold, 3: italic, 4: bold/italic, 5:symbol)
                   labelColor  = "white",                # Label Color
                   lineColor   = "accent1",              # Line Color
                   shapeEnd    = c("none","sphere","circle","square","csquare","rectangle","crectangle","vrectangle"),
                   sizeEnd     = 10,
                   ...){
    ## rxPlot returns nothing, but plots a diagram
    if (class(RxODEobj) == "RxODE"){
        x <- RxODEobj$cmpMgr$rxDll()
    } else if (class(RxODEobj) == "RxCompilationManager") {
        x <- RxODEobj$rxDll()
    } else if (class(RxODEobj) == "rxDll"){
        x <- RxODEobj;
    } else {
        stop("The RxODEobj is not a supported RxODE object")
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package igraph needed for this function to work. Please install it.",
             call. = FALSE)
    }
    ig <- igraph(x,family=family,tk=interactive,
                 shape = shape, size = size, colors = colors, fillColor = fillColor, font = font,
                 labelColor = labelColor, lineColor = lineColor, shapeEnd = shapeEnd, sizeEnd = sizeEnd);
    layout <- -igraph::layout.grid(ig);
    
    if (!interactive){
        plot(ig,edge.label.family=family, layout = layout, ...);
    } else {
        igraph::tkplot(ig,edge.label.family=family, layout = layout, ...);
    }
    invisible(NULL)
} # end function rxPlot

#' Plot a model digram for simple ODE models
#'
#' This plots a line diagram for the current system of differential
#' equations.
#'
#' @param x is an RxODE object
#' @param ... are other arguments
#'
#' This function is an alias for the \code{\link{rxPlot}} function.
#'
#' @seealso \code{\link{rxPlot}},\code{\link{RxODE}}.
#' @export
plot.RxODE <- function(x,...){
    rxPlot(x$cmpMgr$rxDll(),...);
}
#' @rdname plot.RxODE
#' @export
plot.rxDll <- function(x,...){
    rxPlot(x,...);
}
#' @rdname plot.RxODE
#' @export
plot.RxCompilationManager <- function(x,...){
    rxPlot(x,...);
}

#' A compilation manager for RxODE models
#'
#' This function parses, compiles, links, and loads the shared object
#' that implements an RxODE model.
#'
#' The function parses and compiles (if needed) the \code{RxODE}
#' model specified in the string \code{model} into a dynamic link
#' library (DLL on Windows) or a shared object (\code{*.so} on
#' Unix-like systems).
#' 
#' It then dynamically loads this code into the current R
#' session. (Models previously parsed and compiled in previous R
#' sessions only need to be dynamically loaded into the current R session.)
#'
#' @param model a string containing the either the source or the file
#'     name of the \code{RxODE} model
#' @param modName a string with a model identifier (this string will
#'     be used to construct filenames and C symbols, therefore it must
#'     be suitable for creating those identifiers); Model
#'     specification or model file
#' @param wd a string with a file path to a working directory where to
#'     create various C and object files.
#' @param extraC a string with a file path to an extra C file to be
#'     included in the model.  This is useful for adding your own C
#'     functions.
#' @param debug a boolean specifying if debugging is used for the
#'     compiled code by adding -D__DEBUG__ to the compile-time
#'     options.
#' @param mmod A boolean telling if the modName from \code{RxODE} was
#'     missing.  This affects how the model is created and used.
#' @return  An object (closure) with the following member functions:
#' \item{parse}{
#'     this function parses (translates) the ODE-based model 
#'     specification and generates a C file to implement the model.}
#' \item{compile}{
#'     compiles the generated C file for the ODE system
#'     and dynamically loads the machine code in the shared object.}
#' \item{dynLoad}{
#'     if needed, dynamically loads the dynamic library
#'     produced in the \code{compile()} step.  Note that the shared
#'     object persists across R sessions, thus the \code{dynLoad} needs
#'     to be issued as needed.}
#' \item{dynUnload}{
#'     this function unloads the previously dynamically loaded
#'     model object code.  Mainly for internal use.}
#' \item{ode_solver}{
#'     a string with the name of the C symbol for this model solver.}
#' \item{dllfile}{
#'     a string with the name of the dynamic link (or shared object) file.}
#' \item{get.modelVars}{
#'     function that returns a list with 3 character
#'     vectors, \code{params}, \code{state}, and \code{lhs} of variable 
#'     names (identifiers) used in the model specification. 
#'     These will be output when the model is computed (i.e., the ODE solved).}
#' \item{isValid}{
#'     a function that (naively) checks for model validity,
#'     namely that the C object code reflects the latest model 
#'     specification.}
#' \item{get.index}{
#'     helper function to extract the index of one or
#'     more system variables (state, parameter, or other).}
#' \item{getObj}{
#'     internal (not user callable) function.}
#' @examples
#' \dontrun{
#'   cmpMgt <- rx.initCmpMgr(model, "tst1", wd = ".")
#' }
#' @keywords internal models ODE
#' @concept ordinary differential equations
#' @seealso \code{\link{RxODE}}
#' @export rx.initCmpMgr
rx.initCmpMgr <-
    function(model, modName, wd, extraC = NULL,debug = TRUE,mmod = FALSE)
{
    ## Initialize the RxODE compilation manager.  This is a stub
    ## function for backward compatability.
    .model <- model;
    .mmod <- mmod;
    .modName <- modName;
    .wd <- wd;
    .parsed <- FALSE;
    .compiled <- FALSE;
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
        if(!do.it)
            return(invisible(.parsed));
        .parsed <<- TRUE;
        .compiled <<- FALSE;
        invisible(.parsed);
    }
    
    compile <- function(force = FALSE){
        do.it <- force|| !.compiled;
        if(!do.it)
            return(invisible(.compiled));
        lwd <- getwd();
        if (!file.exists(.wd))
            dir.create(.wd, recursive = TRUE)
        if (!file.exists(.wd))
            setwd(.wd);
        on.exit(setwd(lwd));
        if(.mmod){
            .rxDll <<- rxCompile(.model, extraC = .extraC, debug = .debug);

        } else {
            .rxDll <<- rxCompile(.model, .mdir, extraC = .extraC, debug = .debug,
                                 modName = .modName);
        }
        if (class(.rxDll) == "rxDll"){
            .compiled <<- TRUE;
        }
        invisible(.compiled);
    }
   
    dynLoad <- function(force = FALSE){
        ## NB: we may need to reload (e.g., when we re-start R and
        ## re-instantiate the RxODE object from a save.image.
        return(rxLoad(.rxDll));
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
            return(rxState(.rxDll,s));
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
            get.modelVars = function() rxModelVars(.rxDll)[c("params","state","lhs")],
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

#' @rdname print.RxODE
#' @export
"print.RxCompilationManager" <-
function(x, ...)
{
    modName <- x$getObj(".modName")
    cat(sprintf("RxCompilationManager for RxODE model '%s'\n", modName))
   invisible(x)
}

rxShortPath <- function(p, # Path
                        ...){
    ## rxShortPath returns a short path name under windows
    return(gsub("\\\\", "/", utils::shortPathName(p)))
} # end function rxShortPath

rxDvode <- file.path(system.file("common", package = "RxODE"), "call_dvode.c");
rxInclude <- system.file("include", package = "RxODE");
rxLibs <- system.file(file.path("libs", .Platform$r_arch), package = "RxODE");
if (!file.exists(rxInclude)){
    ## Testing instead of running
    rxInclude <- file.path(getwd(),"src","ode");
    if (file.exists(rxInclude)){
        rxDvode <- file.path(getwd(),"inst","common","call_dvode.c");
        rxLibs <- rxInclude;
    } else {
        stop("Required files not found in the package.")
    }
}
if (.Platform$OS.type=="windows"){
    rxDvode <- rxShortPath(rxDvode);
    rxInclude <- rxShortPath(rxInclude);
    rxLibs <- rxShortPath(rxLibs);
}


rxPrefix <- function(model,          # Model or file name of model
                     modName = NULL, # Model name, overrides calculated model name.
                     ...){
    ## rxPrefix returns a prefix for a model.
    if (!is.null(modName)){
        modelPrefix <- sprintf("%s_",gsub("\\W", "_", modName));
    } else if (file.exists(model)){
        modelPrefix <- sprintf("%s_",gsub("\\W", "_", gsub("[.].*$","",base::basename(model))));
    } else {
        parseModel <- tempfile();
        cFile <- tempfile();
        on.exit({unlink(parseModel); unlink(cFile)});
        sink(parseModel);
        cat(model);
        cat("\n");
        sink();
        trans <- rxTrans(parseModel,cFile)
        modelPrefix <- sprintf("rx_%s_",trans["parsed_md5"]);
    }
    modelPrefix <- sprintf("%s%s_",modelPrefix,.Platform$r_arch);
    return(modelPrefix);
} # end function rxPrefix

#' Return the md5 of an RxODE object or file
#'
#' This md5 is based on the model and possibly the extra c code
#' supplied for the model
#' 
#' @param model either a string representing the path to a model file
#'     or a RxODE object
#' 
#' @param extraC The extra C file used to help create the md5.  If the
#'     model is an RxODE object, this argument is ignored.
#'
#' @param ... ignored arguments
#'
#' @return If this is a RxODE object, return a named list:
#'
#' \item{\code{file_md5}}{ is the model's file's md5}
#' 
#' \item{\code{parsed_md5}}{ is the parsed model's file's md5.}
#'
#' Otherwise return the md5 based on the arguments provided
#'
#' @keywords internal
#' @export rxMd5
rxMd5 <- function(model,         # Model File
                  extraC = NULL, # Extra C
                  ...){
    ## rxMd5 returns MD5 of model file.
    ## digest(file=TRUE) includes file times, so it doesn't work for this needs.
    if (class(model) == "character"){
        if (file.exists(model)){
            ret <- readLines(model);
            mod <- paste(ret,collapse="\n");
        } else {
            stop("Requires model to be a file.");
        }
        if (class(extraC) == "character"){
            if (file.exists(extraC)){
                ret <- c(ret,gsub(rex::rex(or(any_spaces,any_newlines)),"",readLines(extraC),perl=TRUE));
            }
        }
        return(list(text=mod,
                    digest=digest::digest(ret)));
    } else {
        rxModelVars(model)$md5;
    }
} # end function rxMd5
#' Translate the model to C code if needed
#'
#' This function translates the model to C code, if needed
#'
#' 
#' @param model This can be either a string specifying a file name for
#'     the RxODE code, or an RxODE family of objects
#'
#' @param cFile The C file where the code should be output
#'
#' @param extraC Extra c code to include in the model.  This can be
#'     useful to specify functions in the model.  These C functions
#'     should usually take \code{double} precision arguments, and
#'     return \code{double} precision values.
#'
#' @param modelPrefix Prefix of the model functions that will be
#'     compiled to make sure that multiple RxODE objects can coexist
#'     in the same R session.
#'
#' @param md5 Is the md5 of the model before parsing, and is used to
#'     embed the md5 into dll, and then provide for functions like
#'     \code{\link{rxModelVars}}.
#' 
#' @param modName is a string specifying the model name.  This string
#'     is used to generate the model's dll file.  If unspecified, and
#'     the model does not come from the file, the model dll name is
#'     based on the parsed md5.
#'
#' @param modVars returns the model variables instead of the named
#'     vector of translated properties.
#'
#' @param ... Ignored parameters.
#'
#' @return a named vector of translated model properties
#'       including what type of jacobian is specified, the \code{C} function prefixes,
#'       as well as the \code{C} functions names to be called through the compiled model.
#' @seealso \code{\link{RxODE}}, \code{\link{rxCompile}}.
#' @export
rxTrans <- function(model,
                    cFile       = sprintf("%s.c",gsub("[.][^.]*$","",model)), # C file output
                    extraC      = NULL,                                       # Extra C file(s)
                    modelPrefix = "",                                         # Model Prefix
                    md5         = "",                                         # Md5 of model
                    modName     = NULL,                                       # Model name for dll
                    modVars     = FALSE,                                      # Return modVars
                    ...){
    ## rxTrans returns a list of compiled properties
    if (class(model) == "character"){
        if (missing(modelPrefix)){
            modelPrefix <- rxPrefix(model,modName);
        }
        if (file.exists(model)){
            if (missing(md5)){
                md5 <- rxMd5(model,extraC)$digest;
            }
        } else {
            stop("This only translates a file (currently; Try rxCompile).");
        }
        parseModel <- tempfile();
        on.exit(unlink(parseModel));
        ret <- .Call(trans, model, cFile, extraC, modelPrefix,md5,parseModel);
        if (file.exists(cFile)){
            ret$md5 <- c(file_md5=md5,parsed_md5 = rxMd5(parseModel,extraC)$digest);
            if (modVars){
                return(ret)
            } else {
                return(c(ret$trans,ret$md5));
            }
        } else {
            stop("Syntax Error (see above)");
        }
        
    } else {
        mv <- rxModelVars(model)
        if (modVars){
            return(mv);
        } else {
            return(c(mv$trans,mv$md5));
        }
    }
} # end function rxTrans

rxTransMakevars <- function(rxProps,                                                                              # rxTrans translation properties
                            compileFlags =c("parsed_md5","ode_solver","model_vars","calc_lhs","calc_jac","dydt"), # List of compile flags
                            debug        = FALSE,                                                                 # Debug compile?
                            ...){
    ## rxTransCompileFlags returns a string for the compiler options
    neededProps <- c("jac",compileFlags);
    if (all(neededProps %in% names(rxProps))){
        if (rxProps["jac"] == "fulluser"){
            ret <- " -D__JT__=1 -D__MF__=21";
        } else if (rxProps["jac"] == "fullint"){
            ret <- " -D__JT__=2 -D__MF__=22";
        }
        tmp <- rxProps[compileFlags];
        tmp["parsed_md5_str"] <- sprintf("\"\\\"%s\\\"\"",tmp["parsed_md5"]);
        ret <- paste(c(ret,sprintf("-D__%s__=%s",toupper(names(tmp)),tmp)),collapse=" ");
        if (debug){
            ret <- sprintf("%s -D__DEBUG__",ret);
        }
        ret <- sprintf("PKG_CPPFLAGS=-I%s %s\nPKG_LIBS=-L%s -lodeaux $(BLAS_LIBS) $(FLIBS)",rxInclude,ret,rxLibs);
        cat(ret);
        return(ret);
    } else {
        stop("Cannot figure out what needs to be specified in the compiler.")
    }
} # end function rxTransCompileFlags

#' Compile a model if needed
#'
#' This is the compilation workhorse creating the RxODE model dll
#' files.
#'
#' @param model This can be either a string specifying a file name for
#'     the RxODE code, a string representing the model specification,
#'     or an RxODE family of objects to recompile if needed.
#'
#' @param dir This is the model directory where the C file will be
#'     stored for compiling.
#' 
#'     If unspecified, the C code is stored in a temporary directory,
#'     then the model is compiled and moved to the current directory.
#'     Afterwards the C code is removed.
#'
#'     If specified, the C code is stored in the specified directory
#'     and then compiled in that directory.  The C code is not removed
#'     after the dll is created in the same directory.  This can be
#'     useful to debug the c-code outputs.
#'
#' @param prefix is a string indicating the prefix to use in the C
#'     based functions.  If missing, it is calculated based on file
#'     name, or md5 of parsed model.
#'
#' @param extraC Extra c code to include in the model.  This can be
#'     useful to specify functions in the model.  These C functions
#'     should usually take \code{double} precision arguments, and
#'     return \code{double} precision values.
#'
#' @param force is a boolean stating if the (re)compile should be
#'     forced if RxODE detects that the models are the same as already
#'     generated.
#'
#' @param modName is a string specifying the model name.  This string
#'     is used to generate the model's dll file.  If unspecified, and
#'     the model does not come from the file, the model dll name is
#'     based on the parsed md5.
#' @param ... Other arguments sent to the \code{\link{rxTrans}} function.
#'
#' @return A rxdll object that has the following slots
#' 
#' \item{dll}{dll path}
#' \item{model}{model specification}
#' \item{.c}{A function to call C code in the correct context from the dll
#'          using the \code{\link{.C}} function.}
#' \item{.call}{A function to call C code in the correct context from the dll
#'          using the \code{\link{.Call}} function.}
#' \item{args}{A list of the arguments used to create the rxDll object.}
#' @seealso \code{\link{RxODE}}
#' @export
rxCompile <-  function(model,           # Model
                       dir,             # Directory
                       prefix,          # Prefix
                       extraC  = NULL,  # Extra C File.
                       force   = FALSE, # Force compile
                       modName = NULL,  # Model Name 
                       ...){
    ## rxCompile returns the dll name that was created.
    if (class(model) == "character"){
        dllCopy <- FALSE;
        if (missing(dir)){
            dir <- tempfile();
            dllCopy <-  TRUE;
            on.exit(unlink(dir, recursive = TRUE))
        }
        if (missing(prefix)){
            prefix <- rxPrefix(model, modName);
        }
        if(!file.exists(dir))
            dir.create(dir, recursive = TRUE)
        cFile <- file.path(dir,sprintf("%s.c",substr(prefix,0,nchar(prefix)-1)));
        cDllFile <- file.path(dir,sprintf("%s%s",substr(prefix,0,nchar(prefix)-1), .Platform$dynlib.ext));
        if (dllCopy){
            finalDll <- file.path(getwd(),basename(cDllFile));
        } else {
            finalDll <-  cDllFile;
        }
        if (!file.exists(model)){
            mFile <- sprintf("%s.rx",substr(cFile,0,nchar(cFile)-2));
            sink(mFile);
            cat(model);
            cat("\n");
            sink();
        } else {
            mFile <- model;
        }
        md5 <- rxMd5(mFile,extraC);
        allModVars <- NULL;
        needCompile <- TRUE
        if (file.exists(finalDll)){
            try(dyn.load(finalDll, local = FALSE), silent=TRUE);
            modVars <- sprintf("%smodel_vars",prefix);
            if (is.loaded(modVars)){
                allModVars <- eval(parse(text=sprintf(".Call(\"%s\")",modVars)),envir=.GlobalEnv)
                modVars <- allModVars$md5;
            }
            if (modVars["file_md5"] == md5$digest){
                needCompile <- FALSE;
            }
        }
        if (force || needCompile){
            dvODE <- file.path(dir,"call_dvode.c");
            dvODEo <- file.path(dir,"call_dvode.o");
            file.copy(rxDvode, dvODE);
            Makevars <- file.path(dir,"Makevars");
            trans <- rxTrans(mFile,cFile = cFile,md5=md5$digest,extraC=extraC,...,modelPrefix=prefix);
            if (file.exists(finalDll)){
                if (modVars["parsed_md5"] == trans["parsed_md5"]){
                    cat("Don't need to recompile, minimal change to model detected.\n");
                    needCompile <- FALSE;
                }
            }
            if (force || needCompile){
                ## Setup Makevars
                if (file.exists(Makevars)){
                    unlink(Makevars);
                }
                if (file.exists(dvODEo)){
                    unlink(dvODEo);
                }
                sink(Makevars);
                cat(rxTransMakevars(trans,...));
                sink();
                sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths    
                ## Change working directory
                owd <- getwd();
                on.exit(setwd(owd));
                setwd(dir);
                try(dyn.unload(finalDll), silent=TRUE);
                try(unlink(finalDll));
                cmd <- sprintf("%s/bin/R CMD SHLIB %s %s", 
                               Sys.getenv("R_HOME"), base::basename(cFile),
                               base::basename(dvODEo));
                cat(sprintf("%s\n",cmd));
                rc <- try(do.call(sh, list(cmd)), silent = FALSE)
                if(inherits(rc, "try-error"))
                    stop(sprintf("error compiling %s", cFile));
                if (dllCopy){
                    file.copy(cDllFile,finalDll);
                }
                try(dyn.load(finalDll, local = FALSE), silent=TRUE);
                modVars <- sprintf("%smodel_vars",prefix);
                if (is.loaded(modVars)){
                    allModVars <- eval(parse(text=sprintf(".Call(\"%s\")",modVars)),envir=.GlobalEnv)
                }
            }
        }
        .c <- function(...){return(.C(...))};
        .call <- function(...){return(.Call(...))};
        args <- list(model = model, dir = dir, prefix = prefix,
                     extraC = extraC, force = force, modName = modName,
                     ...);
        ret <- list(dll     = finalDll,
                    model   = md5$text,
                    extra   = extraC,
                    modVars = allModVars,
                    .c      = .c,
                    .call   = .call,
                    args    = args);
        class(ret) <- "rxDll";
        return(ret);
    } else if (class(model) == "rxDll") {
        args <- model$args;
        if (!missing(dir)){
            args$dir <- dir;
        }
        if (!missing(prefix)){
            args$prefix <- prefix;
        }
        if (!missing(extraC)){
            args$extraC <- extraC;
        }
        if (!missing(force)){
            args$force <- force;
        }
        if (!missing(modName)){
            args$modName <- modName;
        }
        ret <- rxCompile(model  = args$model,  dir     = args$dir,
                         prefix = args$prefix, extraC  = args$extra,
                         force  = args$force,  modName = args$modName,
                         ...);
        return(ret);
    }
} # end function rxCompile

#' Determine if the rxDll is loaded or not.
#'
#' @param x is a RxODE family of objects
#'
#' @param retry is a flag to retry to load if the function can't
#'     determine if the object is loaded or not...
#'
#' @return a boolean stating if the dll is loaded
#'
#' @export
rxDllLoaded <- function(x,retry = TRUE){
    m <- rxTrans(x);
    if (any(names(m) == "ode_solver")){
        return(is.loaded(m["ode_solver"]));
    } else if (retry) {
        m <- rxCompile(x,force=FALSE);
        return(rxDllLoaded(m,retry = FALSE))
    } else {
        stop("Can't figure out if the object...");
    }
}

#' Return the dll associated with the RxODE object
#'
#' This will return the dynamic load library or shared object used to
#' run the C code for RxODE.
#'
#' @param obj A RxODE family of objects or a character string of the
#'     model specification or location of a file with a model
#'     specification.
#'
#' @return a path of the library
#'
#' @keywords internal
#' @export
rxDll <- function(obj,...){
    UseMethod("rxDll");
}

#' @rdname rxDll
#' @export
rxDll.character <- function(obj,...){
    return(rxDll(rxCompile(obj,...)))
}

#' @rdname rxDll
#' @export
rxDll.rxDll <- function(obj,...){
    return(obj$dll)
}

#' @rdname rxDll
#' @export
rxDll.RxODE <- function(obj,...){
    return(rxDll(obj$cmpMgr$rxDll()))
}

#' Load the dll for the object
#'
#' This loads the dll into the current R session to allow C functions
#' to be called in R.
#'
#' @param obj a RxODE family of objects
#'
#' @export
rxLoad <- function(obj){
    if (!(rxDllLoaded(obj))){
        dll <- rxDll(obj);
        rc <- try(dyn.load(dll), silent = TRUE)
        if(inherits(rc, "try-error"))
            stop(sprintf("error loading dll file %s", dll));
    }
    return(invisible());
}

#' Unload the dll for the object
#'
#' This unloads the dll in the R session so that the dll can be
#' deleted.  All the c functions will no longer be accessible.
#'
#' @param obj a RxODE family of objects
#'
#' @export
rxUnload <- function(obj){
    if ((rxDllLoaded(obj))){
        dll <- rxDll(obj);
        rc <- try(dyn.unload(dll), silent = TRUE)
        if(inherits(rc, "try-error"))
            stop(sprintf("error unloading dll file %s", dll));
    }
    return(invisible());
}

#' Delete the dll for the model
#'
#' This function deletes the dll, but doesn't delete the model
#' information in the object.
#'
#' @param obj RxODE family of objects
#'
#' @return a boolean stating if the operation was successful.
#' 
#' @export
rxDelete <- function(obj){
    dll <- rxDll(obj);
    rxUnload(obj)
    unlink(dll);
    return(file.exists(dll));
}

#' Parameters specified by the model
#'
#' This return the model's parameters that are required to solve the
#' ODE system.
#'
#' @param obj is a RxODE family of objects
#' @param ... Ignored arguments
#'
#' @return a character vector listing the parameters in the model.
#'
#' @export
rxParams <- function(obj,...){
    return(rxModelVars(obj)$params);
}

#' @rdname rxParams
#' @export
rxParam <- rxParams

#' State variables
#'
#' This returns the model's compartments or states.
#'
#' @param obj is a RxODE family of objects
#'
#' @param state is a string indicating the state or compartment that
#'     you would like to lookup.
#'
#' @param ... Ignored arguments
#'
#' @return If state is missing, return a character vector of all the states.
#'
#' If state is a string, return the compartment number of the named state.
#'
#' @seealso \code{\link{RxODE}}
#' 
#' @export
rxState <- function(obj,state,...){
    if (missing(state)){
        return(rxModelVars(obj)$state);
    } else {
        objState <- rxState(obj);
        if (length(objState) == 1)
            warning("only one state variable should be input", immediate = TRUE);
        w <- which(objState == state)
        if (length(w) != 1){
            stop(sprintf("Cannot locate state \"%s\"",state));
        }
        return(w);
    }
}
#' Left handed Variables
#'
#' This returns the model calculated variables
#'
#' @param obj RxODE family of objects
#' @param ... ignored arguments
#'
#' @return a character vector listing the calculated parameters
#' @seealso \code{\link{RxODE}}
#' 
#' @export
rxLhs <- function(obj,...){
    return(rxModelVars(obj)$lhs);
}

#' All model variables for a RxODE object
#'
#' Return all the known model variables for a specified rxode object
#'
#' These items are only calculated after compilation; they are
#' built-into the RxODE compiled dll.
#'
#' @param obj RxODE family of objects
#'
#' @return A list of RxODE model properties including:
#'
#' \item{params}{ a character vector of names of the model parameters}
#' \item{lhs}{ a character vector of the names of the model calculated parameters}
#' \item{state}{ a character vector of the compartments in RxODE object}
#' \item{trans}{ a named vector of translated model properties
#'       including what type of jacobian is specified, the \code{C} function prefixes,
#'       as well as the \code{C} functions names to be called through the compiled model.}
#' \item{md5}{a named vector that gives the digest of the model (\code{file_md5}) and the parsed model (\code{parsed_md5})}
#' \item{model}{ a named vector giving the input model (\code{model}),
#'    normalized model (no comments and standard syntax for parsing, \code{normModel}),
#'    and interim code that is used to generate the final C file \code{parseModel}}
#'
#' @keywords internal
#' @export
rxModelVars <- function(obj,...){
    UseMethod("rxModelVars");
}

#' @rdname rxModelVars
#' @export
rxModelVars.rxDll <- function(obj,...){
    return(obj$modVars)
}

#' @rdname rxModelVars
#' @export
rxModelVars.RxCompilationManager <- function(obj,...){
    return(rxModelVars.rxDll(obj$rxDll()))
}

#' @rdname rxModelVars
#' @export
rxModelVars.RxODE <- function(obj,...){
    return(rxModelVars.rxDll(obj$cmpMgr$rxDll()))
}

#' @rdname rxModelVars
#' @export
rxModelVars.solveRxDll <- function(obj,...){
    lst <- attr(obj,"solveRxDll");
    return(rxModelVars.rxDll(lst$object));
}

#' Print rxDll object
#'
#' This tells if the rxDll is loaded, ready and/or deleted.
#'
#' @keywords internal
#' @export print.rxDll
print.rxDll <- function(x,...){
    if (file.exists(x$dll)){
        cat(sprintf("RxODE dll named \"%s\"",basename(x$dll)));
        if (rxDllLoaded(x)){
            cat(" is loaded and ready to use.\n");
        } else {
            cat(" is not loaded now.\n");
        }
    } else {
        cat(sprintf("RxODE dll named \"%s\" has been deleted.\n",basename(x$dll)));
    }
    invisible(x);
}

#' Summary of rxDll object
#'
#' This gives expanded information about the rxDll object
#'
#' @param object RxDll object
#'
#' @param ... Other arguments.  Includes \code{noprint}, which is a
#'     logical telling if the object should print the rxDll object
#'     first. By default this is FALSE
#' 
#' 
#' @keywords internal
#' @export summary.rxDll
summary.rxDll <- function(object,...){
    args <- as.list(match.call(expand.dots = TRUE));
    if (any(names(args) == "noprint")){
        noprint <- args$noprint;
    } else {
        noprint <- FALSE;
    }
    if (!noprint)
        print(object);
    cat(sprintf("dll: %s\n",rxDll(object)));
    cat(sprintf("Jacobian: %s\n",ifelse(rxModelVars(object)$jac == "fulluser","Full User Specified","Full Internally Caluclated")));
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

#' Initial Values and State values for a RxODE object
#'
#' Returns the initial values of the rxDll object
#'
#' @param rxDllObj rxDll, RxODE, or named vector representing default
#'     initial arguments
#'
#' @param vec If supplied, named vector for the model.
#'
#' @param req Required names, and the required order for the ODE solver
#'
#' @param default a number or NA representing the default value for
#'     parameters missing in \code{vec}, but required in \code{req}.
#'
#' @param noerror is a boolean specifying if an error should be thrown
#'     for missing parameter values when \code{default} = \code{NA}
#'
#' @keywords internal
#' @export
rxInits <- function(rxDllObj,        # rxDll object
                    vec,             # Supplied parameters
                    req,             # Required names, and order
                    default = 0,     # Default value; If NA, then throw error if doesn't exist
                    noerror = FALSE, # no error thrown for NA
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
    
    ini <- rxModelVars(rxDllObj)$ini;
    if (!missing(vec)){
        if (length(vec) > 0){
            nv <- names(vec)
            nr <- names(ini)
            if (is.null(nv)){
                if (!missing(req) && length(req) == length(vec)){
                    warning(sprintf("Assumed order of inputs: %s",paste(req,collapse=", ")))
                    return(vec)
                } else {
                    stop(sprintf("Length mismatch\nreq: c(%s)\nvec: c(%s)\n%s",paste(req,collapse=", "),paste(vec,collapse=", "),rxModelVars(rxDllObj)))
                }
            } else {
                shared <- vec[nv %in% nr];
                new <- vec[!(nv %in% nr)];
                old <- ini[!(nr %in% nv)];
                vec <- c(shared,old,new);
            }
        } else {
            vec <- ini;
        }
    } else {
        vec <- ini;
    }
    ret <- vec;
    if (!missing(req)){
        nv <- names(vec);
        if (!is.null(nv)){
            ret <- ret[nv %in% req];
            missing <- req[!(req %in%  names(ret))];
            if (is.na(default) && length(missing) > 0 && !noerror){
                stop(sprintf("Missing the following parameter(s): %s.",paste(missing,collapse=", ")))
            }
            if (length(missing) > 0){
                ret[missing] <- default;
                if (!noerror)
                    warning(sprintf("Assiged %s to %s.",paste(missing,collapse=", "),default))
            }
            ret <- ret[req];
        }
    }
    return(ret);
    
} # end function rxInits

#' @rdname rxInits
#' @export
rxInit <- rxInits;

#' Solves a ODE equation
#'
#' This uses RxODE family of objects, file, or model specification to
#' solve a ODE system.
#'
#' @param object is a either a RxODE family of objects, or a file-name
#'     with a RxODE model specification, or a string with a RxODE
#'     model specification.
#'
#' @param params a numeric named vector with values for every
#'     parameter in the ODE system; the names must correspond to the
#'     parameter identifiers used in the ODE specification;
#'
#' @param events an \code{EventTable} object describing the input
#'     (e.g., doses) to the dynamic system and observation sampling
#'     time points (see \code{\link{eventTable}});
#'
#' @param inits a vector of initial values of the state variables
#'     (e.g., amounts in each compartment), and the order in this
#'     vector must be the same as the state variables (e.g., PK/PD
#'     compartments);
#'
#' @param covs a matrix or dataframe the same number of rows as the
#'     sampling points defined in the events \code{EventTable}.  This
#'     is for time-varying covariates.
#'
#' @param covs_interpolation specifies the interpolation method for
#'     time-varying covariates. When solving ODEs it often samples
#'     times outside the sampling time specified in \code{events}.
#'     When this happens, the time varying covariates are
#'     interpolated.  Currently this can be \code{"Linear"}
#'     interpolation (the default), which interpolates the covariate
#'     by solving the line between the observed covariates and
#'     extrapolating the new covariate value. The other possibility is
#'     \code{"LOCF"}, or Last observation carried forward.  In this
#'     approach, the last observation of the covariate is considered
#'     the current value of the covariate.
#'
#' @param stiff a logical (\code{TRUE} by default) indicating whether
#'     the ODE system is stifff or not.
#' 
#'     For stiff ODE sytems (\code{stiff=TRUE}), \code{RxODE} uses the
#'     LSODA (Livermore Solver for Ordinary Differential Equations)
#'     Fortran package, which implements an automatic method switching
#'     for stiff and non-stiff problems along the integration
#'     interval, authored by Hindmarsh and Petzold (2003).
#'       
#'     For non-stiff systems (\code{stiff=FALSE}), \code{RxODE} uses
#'     DOP853, an explicit Runge-Kutta method of order 8(5,3) of
#'     Dormand and Prince as implemented in C by Hairer and Wanner
#'     (1993).
#'
#' @param transit_abs boolean indicating if this is a transit
#'     compartment absorption
#'
#' @param atol a numeric absolute tolerance (1e-08 by default);
#'
#' @param rtol a numeric relative tolerance (1e-06 by default).
#'
#' @param maxsteps maximum number of (internally defined) steps allowed
#'     during one call to the solver. (5000 by default)
#'
#' @param hmin The minimum absolute step size allowed. The default
#'     value is 0.
#'
#' @param hmax The maximum absolute step size allowed.  The default
#'     checks for the maximum difference in times in your sampling and
#'     events, and uses this value.  The value 0 is equivalent to
#'     infinite maximum absolute step size.
#'
#' @param maxordn The maximum order to be allowed for the nonstiff
#'     (Adams) method.  The default is 12.  It can be between 1 and
#'     12.
#'
#' @param maxords The maximum order to be allowed for the stiff (BDF)
#'     method.  The default value is 5.  This can be between 1 and 5.
#'
#' @param ... Other arguments including scaling factors for each
#'     compartment.  This includes S#=numeric will scale a compartment
#'     # by a dividing the compartment amount by the scale factor,
#'     like NONMEM.
#'
#' @return An \dQuote{rxSolve} solve object that stores the solved
#'     value in a matrix with as many rows as there are sampled time
#'     points and as many columns as system variables (as defined by
#'     the ODEs and additional assignments in the RxODE model code).
#'     It also stores information about the call to allow dynmaic
#'     updating of the solved object.
#'
#'     The operations for the object are simialar to a data-frame, but
#'     expand the \code{$} and \code{[[""]]} access operators and
#'     assignment operators to resolve based on different parameter
#'     values, initial conditions, solver parameters, or events (by
#'     updaing the \code{time} variable).
#'
#'     You can call the \code{\link{EventTable}} methods on the solved
#'     object to update the event table and resolve the system of
#'     equations.  % Should be able to use roxygen templates...
#' 
#' @references
#'
#' Hindmarsh, A. C.
#' \emph{ODEPACK, A Systematized Collection of ODE Solvers}.
#' Scientific Computing, R. S. Stepleman et al. (Eds.),
#' North-Holland, Amsterdam, 1983, pp. 55-64.
#'
#' Petzold, L. R.
#' \emph{Automatic Selection of Methods for Solving Stiff and Nonstiff
#' Systems of Ordinary Differential Equations}.
#' Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
#'
#' Hairer, E., Norsett, S. P., and Wanner, G.
#' \emph{Solving ordinary differential equations I, nonstiff problems}.
#' 2nd edition, Springer Series in Computational Mathematics,
#' Springer-Verlag (1993).
#' 
#' 
#' @seealso \code{\link{RxODE}}
#' @export
rxSolve <- function(object,                      # RxODE object
                    params,                      # Parameter 
                    events,                      # Events
                    inits              = NULL,   # Initial Events
                    covs               = NULL,   # Covariates
                    stiff              = TRUE,   # Is the system stiff
                    transit_abs        = FALSE,  # Transit compartment absorption?
                    atol               = 1.0e-6, # Absoltue Tolerance for LSODA solver
                    rtol               = 1.0e-6, # Relative Tolerance for LSODA solver
                    maxsteps           = 5000,   # Maximum number of steps
                    hmin               = 0,      # Hmin
                    hmax               = NULL,   # Hmax
                    hini               = 0,      # Hini
                    maxordn            = 12,     # maxordn
                    maxords            = 5,      # maxords
                    ...,
                    covs_interpolation = c("Linear","LOCF")
                    ) {
    ## rxSolve returns
    UseMethod("rxSolve");
} # end function rxSolve

#' @rdname rxSolve
#' @export
rxSolve.solveRxDll <- function(object,params, events, inits, covs, stiff, transit_abs, atol, rtol, maxsteps, hmin, hmax, hini, maxordn, maxords, ...,
                               covs_interpolation= c("Linear","LOCF")){
    call <- as.list(match.call(expand.dots = TRUE));
    lst <- attr(object,"solveRxDll");
    lst <- lst[names(lst) != "matrix"];
    lst <- lst[names(lst) != "object"];
    for (n in names(call)[c(-1,-2)]){
        if (n == "params"){
            for (n2 in names(eval(call$params))){
                lst$params[n2] <-  eval(call$params)[n2];
            }
        } else if(n == "inits"){
            for (n2 in names(eval(call$inits))){
                lst$inits[n2] <- eval(call$inits)[n2]
            }
        } else {
            lst[[n]] <- call[[n]];
        }
    }
    lst$object <- object$object;
    return(do.call(rxSolve.rxDll,lst,envir=parent.frame(1)));
}

#' @rdname rxSolve
#' @export
rxSolve.RxODE <- function(object,params,events,inits = NULL,covs = NULL, stiff = TRUE, transit_abs = FALSE,atol = 1.0e-8, rtol = 1.0e-6, maxsteps=5000,hmin=0,hmax=NULL,hini=0, maxordn=12, maxords=5,...,
                          covs_interpolation = c("Linear","LOCF")){
    rxSolve.rxDll(object$cmpMgr$rxDll(),params,events,inits,covs,stiff, transit_abs,atol,rtol,maxsteps, hmin, hmax, hini, maxordn, maxords,..., covs_interpolation = covs_interpolation)
}
#' @rdname rxSolve
#' @export
rxSolve.RxCompilationManager <- function(object,params,events,inits = NULL, covs = NULL,stiff = TRUE, transit_abs = FALSE,atol = 1.0e-8, rtol = 1.0e-6,maxsteps=5000,hmin=0,hmax=NULL,hini=0, maxordn=12, maxords=5,...,
                                         covs_interpolation = c("Linear","LOCF")){
    rxSolve.rxDll(object$rxDll(),params,events,inits,covs,stiff, transit_abs,atol,rtol,maxsteps,hmin,hmax,hini, maxordn, maxords,...,
                  covs_interpolation = covs_interpolation);
}
#' @rdname rxSolve
#' @export
rxSolve.character <- function(object,params,events,inits = NULL, covs = NULL,stiff = TRUE, transit_abs = FALSE,atol = 1.0e-8, rtol = 1.0e-6,maxsteps=5000,hmin=0,hmax=NULL,hini=0, maxordn=12, maxords=5,...,
                              covs_interpolation = c("Linear","LOCF")){
    rxSolve.rxDll(rxCompile(object),params,events,inits,covs,stiff, transit_abs,atol,rtol,maxsteps,hmin,hmax,hini, maxordn, maxords,...,
                  covs_interpolation = covs_interpolation);
}
#' @rdname rxSolve
#' @export
rxSolve.rxDll <- function(object,params,events,inits = NULL, covs = NULL,stiff = TRUE, transit_abs = FALSE,atol = 1.0e-8, rtol = 1.0e-6,maxsteps=5000,hmin=0,hmax=NULL,hini=0, maxordn=12, maxords=5,...,
                          covs_interpolation = c("Linear","LOCF")){
    ## rxSolve.rxDll returns a solved object
    if (missing(events) && class(params) == "EventTable"){
        events <- params;
        params <- c();
    }
    last.solve.args <-
        list(params = params, events = events$copy(),
             inits = inits, covs = covs, stiff = stiff, 
             transit_abs = transit_abs, atol = atol, rtol = rtol, maxsteps = maxsteps,
             hmin = hmin, hmax = hmax, hini = hini, maxordn = maxordn, maxords = maxords,
             covs_interpolation = covs_interpolation, ...)    ;

    event.table <- events$get.EventTable()
    if (!is.numeric(maxordn)) 
        stop("`maxordn' must be numeric")
    if (maxordn < 1 || maxordn > 12) 
        stop("`maxordn' must be >1 and <=12")
    if (!is.numeric(maxords)) 
        stop("`maxords' must be numeric")
    if (maxords < 1 || maxords > 5) 
        stop("`maxords' must be >1 and <=5")
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
        if(hini < 0)
            stop("`hini' must be a non-negative value")
    } else {
        hini <- 0;
    }
    ## preserve input arguments. 
    inits <- rxInits(object,inits,rxState(object),0);
    params <- rxInits(object,params,rxParams(object),NA,!is.null(covs));
    if (!is.null(covs)){
        cov <- as.matrix(covs);
        pcov <- sapply(dimnames(cov)[[2]],function(x){
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
    s <- as.list(match.call(expand.dots = TRUE)) 
    wh <- grep(pattern="S\\d+$", names(s))[1]
    ## HACK: fishing scaling variables "S1 S2 S3 ..." from params call
    ## to solve(). Maybe define a "scale=c(central=7.6, ...)" argument
    ## similar to "params="?
    scaler.ix <- 0
    if (!is.na(wh)) {
        if (s[[wh]] %in% names(params)) {
            scaler <- params[s[[wh]]]
            scaler.ix <- as.numeric(substring(names(s)[wh], 2))
        } else {
            warning(paste("scaler variable not found:", s[[wh]]))
        }
    }
    state_vars <- rxState(object);
    neq  <- length(state_vars);
    lhs_vars <- rxLhs(object);
    nlhs <- length(lhs_vars);

    ntime <- dim(event.table)[1];
    ret   <- rep(0.0, ntime*neq);
    lhs   <- rep(0.0, ntime*nlhs);
    rc <- as.integer(0)  # return code 0 (success) or IDID in call_dvode.c
    if (is.null(inits))
        inits <- rep(0.0, neq)

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
        warning(sprintf("The initial conditions are at t=%s instead of t=0.",event.table$time[1]))
    }
    counts <- rep(0,3);
    xx <- object$.c(rxTrans(object)["ode_solver"],
                    as.integer(neq),
                    as.double(params),
                    as.double(event.table$time),
                    as.integer(event.table$evid),
                    length(event.table$time),
                    as.double(inits),
                    as.double(event.table$amt[event.table$evid>0]),
                    as.double(ret),
                    as.double(atol),
                    as.double(rtol),
                    as.integer(maxsteps),
                    as.integer(stiff),
                    as.integer(transit_abs),
                    as.integer(nlhs),
                    as.double(lhs),
                    ## parameter covariates
                    as.integer(pcov),
                    as.double(cov),
                    as.integer(n_cov),
                    as.integer(isLocf),
                    ## Solver options
                    as.double(hini), ## Mabye accept H0?
                    as.double(hmin),
                    as.double(hmax),
                    as.integer(maxordn),
                    as.integer(maxords),
                    ## Counts
                    as.integer(counts),
                    ## Return Code
                    rc
                    );
    counts <- xx[[length(xx)-1]];
    rc <- xx[[length(xx)]]
    if(rc!=0)
        stop(sprintf("could not solve ODE, IDID=%d (see further messages)", rc))
    x <- cbind(
        matrix(xx[[8]], ncol=neq, byrow=T),
        if(nlhs) matrix(xx[[15]], ncol=nlhs, byrow=T) else NULL,
        if(n_cov) as.matrix(covs) else NULL
    )
    colnames(x) <- c(state_vars, lhs_vars, covnames)

    if (scaler.ix) {
        x[, scaler.ix] <- x[, scaler.ix]/scaler
    }
    
    ret <- cbind(time=event.table$time, x)[events$get.obs.rec(),];
    ## Ensure the objects have names
    names(inits) <- rxState(object);
    names(params) <- rxParams(object);
    lst <- last.solve.args;
    lst[["inits"]] <- inits;
    lst[["params"]] <-  params;
    lst[["object"]] <- object;
    lst[["matrix"]] <- ret;
    names(counts) <- c("solver","dadt","user_jac");
    lst[["counts"]] <-  counts;
    names(ret) <- dimnames(ret)[[2]]; ## For compatability with tidyr::spread
    length(ret) <- length(dimnames(ret)[[2]])
    class(ret) <- c("solveRxDll");
    attr(ret,"solveRxDll") <- lst;
    return(ret);
} # end function rxSolve.rxDll
 
#' @export
print.solveRxDll <- function(x,...){
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
    lst <- attr(x,"solveRxDll")
    cat("Solved RxODE object\n");
    cat(sprintf("Dll: %s\n\n",rxDll(lst$object)))
    cat("Parameters:\n")
    is.dplyr <- requireNamespace("dplyr", quietly = TRUE);
    w <- which((names(lst$params) %in% names(as.data.frame(x))))
    if (length(w) > 0){
        print(lst$params[-w]);
        cat("\n\nFirst Part of Time Varying Covariates:\n");
        d <- as.data.frame(lst$covs)[,names(lst$params)[w]];
        if (length(w) == 1){
            d <- data.frame(d=d);
            names(d) <- names(lst$params)[w];
        }
        if (!is.dplyr){
            print(head(d),n=n);
        } else {
            print(dplyr::as.tbl(d),n=n,width=width);
        }
    }  else {
        print(lst$params);
    }
    cat("\n\nInitial Conditions:\n")
    print(lst$inits);
    cat("\n\nFirst part of data:\n")
    if (!is.dplyr){
        print(head(as.matrix(x)),n=n);
    } else {
        print(dplyr::as.tbl(x),n=n,width=width);
    }
}

#' @export
summary.solveRxDll <- function(object,...){
    lst <- attr(object,"solveRxDll")
    cat("Solved RxODE object\n");
    cat(sprintf("Dll: %s\n\n",rxDll(lst$object)))
    cat("Model:\n");
    cat(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;\n");
    cat(rxModelVars(object)$model["model"]);
    cat(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;\n");
    cat("Parameters:\n")
    w <- which((names(lst$params) %in% names(as.data.frame(object))))
    if (length(w) > 0){
        print(lst$params[-w]);
        cat("\n\nSummary of Time Varying Covariates:\n");
        d <- as.data.frame(lst$covs)[,names(lst$params)[w]];
        if (length(w) == 1){
            d <- data.frame(d=d);
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


#' @export
as.matrix.solveRxDll <- function(x,...){
    lst <- attr(x,"solveRxDll")
    return(lst$matrix);
}


#' @export
as.data.frame.solveRxDll <- function(x,row.names = NULL, optional = FALSE, ...,
                                     stringsAsFactors = default.stringsAsFactors()){
    return(as.data.frame(as.matrix(x),row.names = row.names, optional = optional,...,
                         stringAsFactors = stringAsFactors));
}


#' @export
#' @importFrom utils head 
head.solveRxDll <- function(x, n = 6L, ...){
    return(utils::head.matrix(as.matrix(x),n = n, ...));
}


#' @export
#' @importFrom utils tail
tail.solveRxDll <- function(x, n = 6L, addrownums = TRUE,...){
    return(utils::tail.matrix(as.matrix(x),n = n, addrownums = addrownums, ...));
}

#' Convert Solved RxODE object to tbl
#' @param x Solved RxDll object
#' @param ... other arguments (ignored)
#' @export as.tbl.solveRxDll
as.tbl.solveRxDll <- function(x,...){
    return(dplyr::as.tbl(as.data.frame(x)));
}

solveRxDll_updateEventTable <- function(obj,name,...){
    cat("Update with new event specification.\n");
    tmp <- attr(obj,"solveRxDll");
    events <- tmp$events;
    events[[name]](...)
    tmp <- rxSolve.solveRxDll(obj,events = eval(events));
    lst <- attr(tmp,"solveRxDll");
    lst$events <- events;
    call <- as.list(match.call(expand.dots = TRUE));
    eval(parse(text=sprintf("attr(%s,\"solveRxDll\") <<- lst",toString(call$obj))));
    invisible()
}

accessComp <- function(obj,arg){
    lst <- obj;
    class(lst) <- "list";
    if (any(names(obj) == arg)){
        return(lst[[arg]]);
    } else {
        if (any(rxState(obj) == gsub(regIni,"",arg))){
            arg <- gsub(regIni,"",arg);
            ret <- rxInits(obj)[arg];
            if (is.na(ret)){
                ret <- NA;
                names(ret) <- arg;
                return(ret)
            } else {
                return(ret)
            }
        } else if (any(rxParams(obj) == arg)){
            ret <- rxInits(obj)[arg];
            if (is.na(ret)){
                ret <- NA;
                names(ret) <- arg;
                return(ret)
            } else {
                return(ret)
            }
        } else {
            return(NULL);
        }
    }
}

#' @export
"[[.RxODE" <- function(obj,arg){
    accessComp(obj,arg)
}

#' @export
"[[.RxCompilationManager" <- function(obj,arg){
    accessComp(obj,arg)
}

#' @export
"[[.rxDll" <- function(obj,arg){
    accessComp(obj,arg)
}


#'@export
"$.RxODE" <- function(obj,arg){
    accessComp(obj,arg)
}

#' @export
"$.RxCompilationManager" <- function(obj,arg){
    accessComp(obj,arg)
}

#' @export
"$.rxDll" <- function(obj,arg){
    accessComp(obj,arg)
}

rxTbl <- function(x,msg){
    if (getOption("RxODE.prefer.tbl",TRUE) && class(x) == "data.frame" && requireNamespace("dplyr", quietly = TRUE)){
        if (!missing(msg)){
            cat(sprintf("Change solved object to dplyr's tbl for %s\n",msg));
        }
        return(dplyr::as.tbl(x))
    } else {
        if (!missing(msg)){
            cat(sprintf("Change solved object to data.frame for %s\n",msg))
        }
        return(x)
    }
}

asTbl <- function(obj){
    if (getOption("RxODE.prefer.tbl",TRUE) && requireNamespace("dplyr", quietly = TRUE)){
        return(dplyr::as.tbl(as.data.frame(obj)));
    } else {
        return(as.data.frame(obj))
    }
}

#' @export
"$.solveRxDll" <-  function(obj,arg, exact = TRUE){
    m <- as.data.frame(obj);
    ret <- m[[arg, exact = exact]];
    if (is.null(ret) & class(arg) == "character"){
        if (arg == "t"){
            return(m[["time"]]);
        } else {
            tmp <- attr(obj,"solveRxDll");
            if (isTRUE(exact)){
                w <- which(names(tmp) == arg);
            } else {
                w <- which(regexpr(rex::rex(start,arg),names(tmp)) != -1)
                if (length(w) == 1 && !any(names(tmp) == arg) && is.na(exact)){
                    warning(sprintf("partial match of '%s' to '%s'",arg,names(tmp)[w]));
                }
            }
            if (length(w) == 1){
                return(tmp[[w]]);
            }
            if (any(names(tmp$param) == arg)){
                return(tmp$param[arg]);
            }
            if (any(names(tmp$init) == gsub(regIni,"",arg))){
                arg <- gsub(regIni,"",arg);
                return(tmp$init[arg]);
            }
            if (any(arg == names(tmp$events))){
                if (substr(arg,0,4) == "get."){
                    return(tmp$events[[arg]]);
                } else {
                    call <- as.list(match.call(expand.dots = TRUE));
                    return(eval(parse(text=sprintf("function(...){solveRxDll_updateEventTable(%s,\"%s\",...)}",toString(call$obj),toString(call$arg)))))
                }
            }
            return(NULL);
        }
    } else {
        return(rxTbl(ret));
    }
}


#' @export
"[.solveRxDll" <- function(x,i,j,drop){
    df <- as.data.frame(x);
    if (!missing(i) && missing(j) && missing(drop)){
        df <- df[i,];
    } else if (missing(i) && !missing(j) && missing(drop)){
        df <- df[,j];
    } else if (!missing(i) && !missing(j) && missing(drop)){
        df <- df[i,j];
    } else if (missing(i) && missing(j) && missing(drop)){
        df <- df[];
    } else if (!missing(i) && missing(j) && !missing(drop)){
        df <- df[i, drop = drop];
    } else if (missing(i) && !missing(j) && !missing(drop)){
        df <- df[,j, drop = drop];
    } else if (!missing(i) && !missing(j) && !missing(drop)){
        df <- df[i,j, drop = drop];
    } else if (missing(i) && missing(j) && !missing(drop)){
        df <- df[drop = drop];
    }
    return(rxTbl(df))
}

#' @export
"[[.solveRxDll" <- function(obj,arg,exact=TRUE,internal = FALSE){
    if (internal){
        tmp <- attr(obj,"solveRxDll");
        return(tmp[[arg,exact=exact]]);
    } else {
        "$.solveRxDll"(obj,arg,exact = exact);
    }
}

#' $ Assign for RxODE solved objects
#' 
#' Assign objects by argumnt as obj$arg <- value
#'
#' This also works as obj[[arg]] <- value
#' 
#' @param obj solveRxDll object
#' @param arg Dollar sign name
#' @param value assigned Value.
#' @seealso \code{\link{rxSolve}}
#' @keywords internal
#' @export
"$<-.solveRxDll" <- function(obj,arg,value){
    if (arg == "t"){
        arg <- "time";
    }
    lst <- attr(obj,"solveRxDll")
    iarg <- gsub(regIni,"",arg);
    if (arg == "time"){
        if (class(value) == "EventTable"){
            cat("Update event table and solved object.\n");
            return(rxSolve(obj,events = eventTable))
        } else if (class(value) == "data.frame"){
        } else if (class(value) == "numeric"){
            cat("Updating sampling times in the event table updating object.\n");
            eventTable <- lst$events$copy();
            eventTable$clear.sampling();
            eventTable$add.sampling(value);
            return(rxSolve(obj,events = eventTable));
        }
    } else if (arg != iarg && any(rxState(obj$object) == iarg) && length(value) == 1){
        cat("Updating object with new initial conditions.\n");
        inits <- c(value);
        names(inits) <- gsub(regIni,"",arg);
        return(rxSolve(obj,inits = inits));
    } else if (any(rxParams(obj$object) == arg)){
        cat("Updating object with new paramter values.\n");
        if (length(value) == 1){
            covs <- as.data.frame(lst$covs);
            if (any(names(covs) == arg)){
                cat(sprintf("Changing time varying covariate \"%s\" to a simple parameter value %s\n",arg,value));
                ncovs <- names(covs);
                covs <- as.data.frame(covs[,names(covs) != arg]);
                names(covs) <- ncovs[ncovs != arg];
                params <- c(value);
                names(params) <- arg;
                return(rxSolve.solveRxDll(obj,params = params,covs = covs));
            } else {
                params <- c(value);
                names(params) <- arg;
                return(rxSolve.solveRxDll(obj,params = params));
            }
        } else if (length(value) == length(lst$events$get.sampling()$time)){
            cat(sprintf("Changing simple parameter \"%s\" to a time-varying covariate.\n",arg));
            covs <- as.data.frame(lst$covs);
            covs[[arg]] <- value;
            return(rxSolve.solveRxDll(obj,covs = covs));
        }
    } else if (arg == "params"){
        cat("Updating object with new paramter values.\n");
        return(rxSolve.solveRxDll(obj,params = value));
    } else if (arg == "inits"){
        cat("Updating object with new initial conditions.\n");
        return(rxSolve.solveRxDll(obj,inits = value));
    } else if (any(arg==names(lst))){
        args <- list(obj,value);
        names(args) <- c("object",arg);
        cat(sprintf("Updating object with new solving argument %s=%s.\n",arg,value))
        return(do.call("rxSolve",args, envir = parent.frame(1)))
    } else {
        df <- as.data.frame(obj);
        df <- "$<-.data.frame"(df,arg,value);
        obj <- rxTbl(df,"assignment");
    }
    return(obj);
}


#' Assign solved objects using the [[]] syntax
#' @param obj solved object
#' @param arg element of solved object
#' @param value value assumed
#' @seealso \code{\link{rxSolve}}
#' @keywords internal
#' @export
"[[<-.solveRxDll" <- function(obj,arg,value){
    return("$<-.solveRxDll"(obj,arg,value=value))
}

#' Assign solved objects using the [] syntax
#' @param obj solved object
#' @param arg1 first argument
#' @param arg2 second argument
#' @param value value assumed
#' @keywords internal
#' @export
"[<-.solveRxDll" <- function(obj,arg1,arg2,value){
    if (any(arg2 == c(1,"t","time")) && missing(arg1)){
        obj$time <- value
        return(obj);
    }
    lst <- attr(obj,"solveRxDll")
    df <- as.data.frame(obj);
    if (missing(arg1) & missing(arg2)){
        df[] <- value;
    } else if (missing(arg1)){
        df[,arg2] <-  value;
    } else if (missing(arg2)){
        df[arg1,] <-  value;
    } else {
        d3[arg1,arg2] <- value;
    }
    return(rxTbl(df,"assignment"))
}

#' Assign rownames to rxSolve object
#'
#' row.names(x) <- value;
#'
#' @param x rxode object
#' @param value value assigned.
#' 
#' @keywords internal
#' @export
"row.names<-.solveRxDll" <- function(x,value){
    lst <- attr(x,"solveRxDll")
    df <- as.data.frame(x);
    m <- as.matrix("row.names<-.data.frame"(df,value));
    lst$matrix <- m;
    attr(x,"solveRxDll") <- lst;
    return(x);
}

#' Get row names of rxSolve object
#'
#' @param x rxSovle object
#' @param ... ignored arguments
#' @keywords internal
#' @export
row.names.solveRxDll <- function(x,...){
    lst <- attr(x,"solveRxDll")
    df <- as.data.frame(x);
    return(row.names(df));
}

#' @export
by.solveRxDll <- function(data, INDICES, FUN, ..., simplify = TRUE){
    by(as.data.frame(data),INDICES, FUN, ...,simplify = simplify);
}

#' @export
#' @importFrom stats aggregate
aggregate.solveRxDll <- function(x, by, FUN, ..., simplify = TRUE){
    aggregate.data.frame(as.data.frame(x), by, FUN, ... , simplify = simplify);
}

#' @export
anyDuplicated.solveRxDll <- function(x, incomparables = FALSE,
                                     fromLast = FALSE, ...){
    anyDuplicated(as.data.frame(x),incomparables, fromLast, ...);
}

#' Get dimensions of rxSolve object
#'
#' @param x rxSolve object
#' 
#' @keywords internal
#' @export
dim.solveRxDll <- function(x){
    return(dim(as.matrix(x)));
}

#' Assign dimensions of rxSolve object
#' 
#' dim(x) <- value
#' 
#' @param x rxSolve object
#' @param value dimensions
#' @keywords internal
#' @export
"dim<-.solveRxDll" <- function(x,value){
    lst <- attr(x,"solveRxDll")
    m <- as.matrix(x);
    dim(m) <- value;
    lst$matrix <- m;
    attr(x,"solveRxDll") <- lst;
    return(x);
}

#' Get dimension names for rxSolve object
#'
#' @param x rxSolve object
#'
#' @keywords internal
#' @export
dimnames.solveRxDll <- function(x){
    return(dimnames(as.matrix(x)));
}

#' Assign dimension names for rxSolve object
#'
#' dimnames(x) <- value
#' 
#' @param x rxSolve object
#' @param value dimension names assigned
#' @keywords internal
#' @export
"dimnames<-.solveRxDll" <- function(x,value){
    lst <- attr(x,"solveRxDll")
    m <- as.matrix(x);
    dimnames(m) <- value;
    lst$matrix <- m;
    attr(x,"solveRxDll") <- lst;
    return(x);
}


#' @export
droplevels.solveRxDll <- function(x, except, ...){
    droplevels(as.data.frame(x),except,...);
}


#' @export
duplicated.solveRxDll <- function(x, incomparables = FALSE,
                                  fromLast = FALSE, nmax = NA, ...){
    duplicated(as.data.frame(x),incomparables, fromLast, nmax,...);
}


#' Edit data frame of rxSolve
#'
#' @param name name of the object
#' @keywords internal
#' @export
edit.solveRxDll <- function(name, ...){
    edit(as.data.frame(name),...);
}


#' @export
is.na.solveRxDll <- function(x){
    is.na(as.data.frame(x));
}


#' @export
# Should this method be exported???
Math.solveRxDll <- function(x,...){
    Math(as.data.frame(x),...);
}


#' @export
rowsum.solveRxDll <- function(x, group, reorder = TRUE, na.rm = FALSE, ...){
    rowsum.data.frame(as.data.frame(x),group, reorder, na.rm,...);
}


#' @export
split.solveRxDll <- function(x, f, drop = FALSE, ...){
    split(as.data.frame(x),f,drop,...);
}


#' @export
"split<-.solveRxDll" <- function(x, f, drop = FALSE, ...,value){
    "split<-"(as.data.frame(x),f,drop,...,value=value);
}
 
#' @export
subset.solveRxDll <- function(x, subset, select, drop = FALSE, ...){
    subset.data.frame(as.data.frame(x),subset,select,drop,...)
}


#' @export
#' @importFrom utils stack
stack.solveRxDll <- function(x, select, ...){
    stack(as.data.frame(x),select,...)
}


#' @export
t.solveRxDll <- function(x){
    t(as.matrix(x))
}


#' @export
#' @importFrom utils unstack
unstack.solveRxDll <- function(x, form, ...){
    unstack(as.data.frame(x),form,...)
}


#' @export
unique.solveRxDll <- function(x,incomparables = FALSE, fromLast = FALSE,...){
    unique(as.data.frame(x),incomparables,fromLast,...);
}


#' @export
#'
within.solveRxDll <- function(data,expr,...){
    within(as.data.frame(data),expr,...)
}

with.solveRxDll <- function(data,expr,...){
    with(as.data.frame(data),expr,...)
}

## FIXME rbind, cbind could be possible...

#' @rdname cbind.solveRxDll
#' @export
rbind.solveRxDll <- function(...){
    stop("rbind is unsupported.  First convert to a data.frame with as.data.frame(x).")
}
#' cbind/rbind solveRxDll
#'
#' Cbind/rbind is disabled for RxOde solved objects.  Use as.data.frame(x)
#' to use cbind/rbind.
#'
#' @param ... ignored parameters
#' 
#' @keywords internal
#' @export
cbind.solveRxDll <- function(...){
    stop("cbind is unsupported.  First convert to a data.frame with as.data.frame(x).")
}

## Might work -- merge

#' dyplr filter_ support
#'
#' @param .data solveRxDll object that filter_ is being applied to.
#'
#' @param ... arguments to filter_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
filter_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(filter_(.data=asTbl(.data),...))
    } else {
        return(filter_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr slice_ support
#'
#' @param .data solveRxDll object that slice_ is being applied to.
#'
#' @param ... arguments to slice_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
#' @export
slice_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(slice_(.data=asTbl(.data),...))
    } else {
        return(slice_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr arrange_ support
#'
#' @param .data solveRxDll object that arrange_ is being applied to.
#'
#' @param ... arguments to arrange_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
arrange_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(arrange_(.data=asTbl(.data),...))
    } else {
        return(arrange_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr select_ support
#'
#' @param .data solveRxDll object that select_ is being applied to.
#'
#' @param ... arguments to select_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
select_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(select_(.data=asTbl(.data),...))
    } else {
        return(select_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr rename_ support
#'
#' @param .data solveRxDll object that rename_ is being applied to.
#'
#' @param ... arguments to rename_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
rename_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(rename_(.data=asTbl(.data),...))
    } else {
        return(rename_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr distinct_ support
#'
#' @param .data solveRxDll object that distinct_ is being applied to.
#'
#' @param ... arguments to distinct_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
distinct_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(distinct_(.data=asTbl(.data),...))
    } else {
        return(distinct_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr mutate_ support
#'
#' @param .data solveRxDll object that mutate_ is being applied to.
#'
#' @param ... arguments to mutate_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
mutate_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(mutate_(.data=asTbl(.data),...))
    } else {
        return(mutate_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr transmute_ support
#'
#' @param .data solveRxDll object that transmute_ is being applied to.
#'
#' @param ... arguments to transmute_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
transmute_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(transmute_(.data=asTbl(.data),...))
    } else {
        return(transmute_(.data=asTbl(.data),...,.dots = .dots))
    }
}
#' dyplr summarise_ support
#'
#' @param .data solveRxDll object that summarise_ is being applied to.
#'
#' @param ... arguments to summarise_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
summarise_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(summarise_(.data=asTbl(.data),...))
    } else {
        return(summarise_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr arrange_ support
#'
#' @param .data solveRxDll object that arrange_ is being applied to.
#'
#' @param ... arguments to arrange_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
arrange_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(arrange_(.data=asTbl(.data),...))
    } else {
        return(arrange_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr rename_ support
#'
#' @param .data solveRxDll object that rename_ is being applied to.
#'
#' @param ... arguments to rename_ method
#'
#' @param .dots dplyr .dots argument
#' 
#' @keywords internal
#' @export
rename_.solveRxDll <- function(.data, ... , .dots){
    if (missing(.dots)){
        return(rename_(.data=asTbl(.data),...))
    } else {
        return(rename_(.data=asTbl(.data),...,.dots = .dots))
    }
}

#' dyplr group_by_ support
#'
#' @param .data solveRxDll object that groupt_by_ is being applied to.
#'
#' @param ... arguments to groupt_by_ method
#'
#' @param .dots dplyr .dots argument
#'
#' @param add dpylr \code{add} argument.  By default, when \code{add =
#'     FALSE}, \code{group_by} will override existing groups.  To
#'     instead add to the existing groups, use \code{add = TRUE}.
#' 
#' @keywords internal
#' @export
group_by_.solveRxDll <- function(.data, ..., .dots, add = FALSE){
    if (missing(.dots)){
        return(group_by_(.data=asTbl(.data),..., add = add));
    } else {
        return(group_by_(.data=asTbl(.data),...,.dots = .dots, add = add));
    }
}

## I'm not sure  these functions make sense for a solved rxDll object..

#' dyplr support of sample_n and sample_frac
#'
#' sample_n samples n from the solved RxODE solved object.
#'
#' sample_frac samples a fraction of the number rows of the solved
#' RxODE object.
#'
#' @param tbl solved RxODE object
#' @param size size of sampled object
#' @param replace Sample with or without replacement?
#' @param weight Sampling weights. This expression is evaluated in the
#'     context of the data frame. It must return a vector of
#'     non-negative numbers the same length as the input. Weights are
#'     automatically standardised to sum to 1.
#' @param .env Environment in which to look for non-data names used in
#'     \code{weight}. Non-default settings for experts only.
#' @keywords internal
#' @export
sample_n.solveRxDll <- function(tbl,size,replace = FALSE, weight = NULL, .env = parent.frame()){
    if (missing(weight)){
        return(sample_n(dplyr::as.tbl(tbl),size = size, replace = replace, .env = .env));
    } else {
        return(sample_n(dplyr::as.tbl(tbl),size = size, replace = replace, weight = weight, .env = .env));
    }
}

#' @rdname sample_n.solveRxDll
#' @export
sample_frac.solveRxDll <- function(tbl, size = 1, replace = FALSE, weight = NULL, .env = parent.frame()){
    if (!missing(weight)){
        return(sample_frac(dplyr::as.tbl(tbl),size = size, replace = replace, weight = weight, .env = .env));
    } else {
        return(sample_frac(dplyr::as.tbl(tbl),size = size, replace = replace, .env = .env));
    }
}

#' tidyr support of gather_ method
#'
#' @param data Solved RxODE object
#' @param key_col,value_col Strings giving names of key and value
#'     columns to create.
#' @param gather_cols Character vector giving column names to be
#'     gathered into pair of key-value columns.
#' @param na.rm If \code{TRUE}, will remove rows from output where the
#'     value column in \code{NA}.
#' @param convert If \code{TRUE} will automatically run
#'     \code{type.convert} on the key column. This is useful if the
#'     column names are actually numeric, integer, or logical.
#' @param factor_key If \code{FALSE}, the default, the key values will
#'     be stored as a character vector. If \code{TRUE}, will be stored
#'     as a factor, which preserves the original ordering of the
#'     columns.
#' @keywords internal
#' @export
gather_.solveRxDll <- function(data, key_col, value_col, gather_cols, na.rm = FALSE, 
                               convert = FALSE, factor_key = FALSE){
    return(gather_(dplyr::as.tbl(data),key_col = key_col, value_col = value_col,
                   gather_cols = gather_cols, na.rm = na.rm, 
                   convert = convert, factor_key = factor_key));
}

#' tidyr separate support
#' 
#' @param data RxODE solved object
#' @param col name of column to split, as string
#' @param into Names of new variables to create as character vector.
#' @param sep Separator between columns.  See tidyr for more
#'     information.
#' @param remove If \code{TRUE}, remove input column from output data
#'     frame.
#' @param convert If 'TRUE', will run 'type.convert' with 'as.is =
#'     TRUE' on new columns.
#' @param extra If \code{sep} is a character vector, this controls
#'     what happens when there are too many pieces.  See tidyr
#'     separate for more details.
#' @param fill If \code{sep} is a character vector, this controls what
#'     happens when there are not enough pieces.  See tidyr separate
#'     for more details.
#' @param ... Ignored
#' @keywords internal
#' @export
separate_.solveRxDll <- function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE, 
                                 convert = FALSE, extra = "warn", fill = "warn", ...){
    return(separate_(data = dplyr::as.tbl(data), col =col, into = into, sep = sep, remove = remove, 
                     convert = convert, extra = extra, fill = fill, ...));
}

#' tidyr unite support
#'
#' @param data solved RxODE object
#'
#' @param col Name of new column as string
#'
#' @param from Names of existing columns as character vector
#'
#' @param sep Separator to use between values
#'
#' @param remove If \code{TRUE}, remove input columns from output data
#'     frame.
#' 
#' @keywords internal
#' @export
unite_.solveRxDll <- function(data, col, from, sep = "_", remove = TRUE){
    return(unite_(data = dplyr::as.tbl(data), col = col, from = from, sep = sep, remove = remove));
}

#' tidyr's spread_ support
#'
#' @param data RxODE solved data.
#'
#' @param key_col,value_col Strings giving names of key and value cols.
#'
#' @param fill If set, missing values will be replaced with this
#'     value. See spread for more details.
#' @param convert If \code{TRUE}, \code{type.convert} with \code{asis
#'     = TRUE} will be run on each of the new columns.  See spread for
#'     more details.
#' @param drop If \code{FALSE}, will keep factor levels that don't
#'     appear in the data, filling in missing combinations with
#'     \code{fill}.
#' @keywords internal
#' @export
spread_.solveRxDll <- function(data, key_col, value_col, fill = NA, convert = FALSE, drop = TRUE){
    return(spread_(data = dplyr::as.tbl(data), key_col = key_col, value_col = value_col,
                   fill = fill, convert = convert, drop = drop));
}

#' Cleanup anonymous dlls
#'
#' This cleans up any dlls created by text files
#'
#' @param wd What directory should be cleand
#'
#' This cleans up all files named rx-*.dll and associated files as
#' well as call_dvode.o and associated files
#'
#' @return TRUE if successful
#' #' 
#' @export
#'
rxClean <- function(wd = getwd()){
    owd <- getwd();
    setwd(wd);
    on.exit(setwd(owd));
    pat <- "^(Makevars|(rx.*|call_dvode)[.](o|dll|s[ol]|c|rx))$"
    files <- list.files(pattern=pat);
    for (f in files){
        try(dyn.unload(f),silent=TRUE);
        unlink(f);
    }
    return(length(list.files(pattern=pat)) == 0);
}

