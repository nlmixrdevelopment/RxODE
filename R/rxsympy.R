regIdentifier1 <- rex::rex(one_of("a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."))
regIdentifier2 <- rex::rex(at_least(".",1), one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."));
regIdentifier <- rex::rex(or(regIdentifier1, regIdentifier2));
regSens <- rex::rex("rx__sens_", capture(regIdentifier), "_",  capture(regIdentifier), "__");
regToSens1 <- rex::rex( capture(regIdentifier), or("_", ".", ":"),  capture(regIdentifier));
regToSens2 <- rex::rex( "d/dt(d", capture(regIdentifier), "/d",  capture(regIdentifier), ")");
regJac <- rex::rex( "df(", capture(regIdentifier), ")/dy(",  capture(regIdentifier), ")");
regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                         group(any_of("0":"9"), ".", some_of("0":"9"))),
                      at_most(group(or("E", "e"), at_most(or("+", "-"), 1), some_of("0":"9")), 1));
regFloat2 <- rex::rex(some_of("0":"9"), or("E", "e"), at_most(or("-", "+"), 1), some_of("0":"9"));
regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
regNum <- rex::rex(at_most("-", 1), or(regDecimalint, regFloat1, regFloat2))
regDDt <- rex::rex(start, "rx__d_dt_", capture(anything), "__", end);
regDfDy <- rex::rex(start, "rx__df_", capture(anything), "_dy_", capture(anything), "__", end);

## Currently unsupported sympy functions
sympy.bad <- c(
    "Chi",
    "Ci",
    "DiracDelta",
    "E1",
    "Eijk",
    "Heaviside",
    "KroneckerDelta",
    "LeviCivita",
    "Li",
    "Shi",
    "Si",
    "Ynm",
    "Ynm_c",
    "Znm",
    "airyai",
    "airyaiprime",
    "airybi",
    "airybiprime",
    "assoc_laguerre",
    "assoc_legendre",
    "bspline_basis",
    "bspline_basis_set",
    "chebyshevt",
    "chebyshevt_poly",
    "chebyshevt_root",
    "chebyshevu",
    "chebyshevu_poly",
    "chebyshevu_root",
    "dirichlet_eta",
    "ei",
    "elliptic_k",
    "erf2",
    "erf2inv",
    "erfcinv",
    "erfi",
    "erfinv",
    "expint",
    "fresnelc",
    "fresnels",
    "gegenbauer",
    "gegenbauer",
    "gegenbauer",
    "gegenbauer_poly",
    "hankel1",
    "hankel2",
    "hermite",
    "hermite_poly",
    "hyper",
    "jacobi",
    "jacobi_normalized",
    "jacobi_poly",
    "jn",
    "jn_zeros",
    "laguerre",
    "laguerre_poly",
    "legendre",
    "legendre_poly",
    "lerchphi",
    "li",
    "loggamma",
    "lowergamma",
    "mathieuc",
    "mathieucprime",
    "mathieus",
    "mathieusprime",
    "meijerg",
    "polylog",
    "uppergamma",
    "yn",
    "zeta")

regBadSymPy <- rex::rex(capture(or(sympy.bad)), "(")

## Start DSL based on http://adv-r.had.co.nz/dsl.html
unaryOp <- function(left, right) {
    force(left)
    force(right)
    function(e1) {
        paste0(left, e1, right)
    }
}

binaryOp <- function(sep) {
    force(sep)
    function(e1, e2) {
        if (missing(e2)){
            paste0(gsub(" ", "", sep), e1)
        } else {
            paste0(e1, sep, e2)
        }
    }
}

binaryOp2 <- function(sep) {
    force(sep)
    function(e1, e2) {
        paste0("(", e1, ")", sep, "(", e2, ")")
    }
}

divOp <- function(){
    function(e1, e2){
        if (e1 == "d" && grepl(rex::rex(start, "__dt__"), e2)){
            paste0("rx__d_dt_",gsub(rex::rex(start, "__dt__"), "", e2));
        } else if (grepl(rex::rex(start, "__df_", anything, "_", end), e1) && grepl(rex::rex(start, "_dy_",anything, "__", end), e2)){
            paste0("rx", substring(e1, 0, nchar(e1) - 1), e2)
        } else {
            paste0(e1, " / ", e2)
        }
    }
}

besselOp <- function(type){
    force(type)
    function(z, nu, df){
        base <- paste0("bessel", type, "(", nu, ",", z, ")");
        if (missing(df)){
            if (any(type == c("i", "k"))){
                stop(sprintf("bessel_%s requires 3 arguments.", type))
            }
            return(base)
        } else if (!any(type == c("i", "k"))){
            stop(sprintf("bessel_%s requires 2 arguments.", type))
        }
        if (df == 1){
            return(base);
        } else if (df == 2) {
            if (type == "i"){
                return(paste0("exp(-", z, ") * ", base));
            } else {
                return(paste0("exp(", z, ") * ", base));
            }
        } else {
            stop(sprintf("The third argument for bessel_%s must be 1 (unscaled) or 2 (scaled).", type))
        }
    }
}

besselOp2 <- function(type){
    force(type);
    function(nu, z){
        paste0("bessel_", type, "(", z,", ", nu, ifelse(any(type == c("i", "k")), ", 1.0", ""), ")");
    }
}

functionOp <- function(fn){
    force(fn);
    function(...){
        paste0(fn, "(", paste(unlist(list(...)), collapse=", "), ")")
    }
}

functionOp2 <- function(fn, end){
    force(fn);
    force(end);
    function(...){
        paste0(fn, paste(unlist(list(...)), collapse=", "), end)
    }
}

functionBrewx <- function(brew){
    force(brew);
    function(x){
        file <- tempfile()
        on.exit(unlink(file));
        sink(file)
        brew::brew(text=brew);
        sink();
        return(readLines(file));
    }
}

functionBrewxy <- function(brew){
    force(brew);
    function(x, y){
        file <- tempfile()
        on.exit(unlink(file));
        sink(file)
        brew::brew(text=brew);
        sink();
        return(readLines(file));
    }
}

sympyRxFEnv <- new.env(parent = emptyenv())
rxSympyFEnv <- new.env(parent = emptyenv())
for (op in c("+", "-", "*")){
    rxSympyFEnv[[op]] <- binaryOp(paste0(" ", op, " "));
    sympyRxFEnv[[op]] <- binaryOp(paste0(" ", op, " "));
}
rxSympyFEnv$"/" <- divOp();
sympyRxFEnv$"/" <- binaryOp(" / ");
rxSympyFEnv$"^" <- binaryOp("**")
rxSympyFEnv$"**" <- binaryOp("**")
sympyRxFEnv$"**" <- binaryOp("^")
sympyRxFEnv$"^" <- binaryOp("^")
rxSympyFEnv$dt <- function(e1){
    paste0("__dt__", e1, "__");
}

rxSympyFEnv$df <- function(e1){
    paste0("__df_", e1, "_")
}

rxSympyFEnv$dy <- function(e1){
    paste0("_dy_", e1, "__")
}

## rx -> sympy
rxSympyFEnv$bessel_i <- besselOp("i");
rxSympyFEnv$bessel_j <- besselOp("j");
rxSympyFEnv$bessel_k <- besselOp("k");
rxSympyFEnv$bessel_y <- besselOp("y");

## sympy -> rx
sympyRxFEnv$besseli <- besselOp2("i");
sympyRxFEnv$besselj <- besselOp2("j");
sympyRxFEnv$besselk <- besselOp2("k");
sympyRxFEnv$bessely <- besselOp2("y");


## rxSympyFEnv$"[" <- binaryOp("_")
##
rxSympyFEnv$"(" <- unaryOp("(", ")")
sympyRxFEnv$"(" <- unaryOp("(", ")")

## pow -> **
rxSympyFEnv$pow <- binaryOp2("**");
sympyRxFEnv$pow <- binaryOp2("^");

## gammafn -> gamma
rxSympyFEnv$gammafn <- functionOp("gamma")
rxSympyFEnv$lgammafn <- functionOp2("log(gamma(", "))")
rxSympyFEnv$lgamma <- functionOp2("log(gamma(", "))")
rxSympyFEnv$tetragamma <- functionOp2("psigamma(", ", 2)")
rxSympyFEnv$pentagamma <- functionOp2("psigamma(", ", 3)")
rxSympyFEnv$lbeta <- functionOp2("log(beta(", "))")
rxSympyFEnv$lgamma1p <- functionOp2("log(gamma((", ")+1))")
rxSympyFEnv$cospi <- functionOp2("cos(pi * (", "))")
rxSympyFEnv$sinpi <- functionOp2("sin(pi * (", "))")
rxSympyFEnv$tanpi <- functionOp2("tan(pi * (", "))")
rxSympyFEnv$logspace_add <- binaryOp2(" + ");
rxSympyFEnv$logspace_sub <- binaryOp2(" - ");
## Following R functions are not translated
## ftrunc, fround, fprec, fsign, sign, fmin2, fmax2, imin2, imax2, logspace_sum
## choose lchoose

rxSympyFEnv$R_pow <- binaryOp2("**");
rxSympyFEnv$R_pow_di <- binaryOp2("**");
rxSympyFEnv$log1p <- functionOp2("log(1 + (", "))");
rxSympyFEnv$log1pmx <- functionBrewx("(log(1 + (<%=x%>))-(<%=x%>))");
rxSympyFEnv$expm1 <- functionOp2("(exp(", ")-1)");

rxSympyFEnv$choose <- functionBrewxy("(factorial(<%=x%>)/(factorial(<%=y%>)*factorial((<%=x%>)-(<%=y%>))))");
rxSympyFEnv$lchoose <- functionBrewxy("(log(gamma((<%=x%>)+1))-log(gamma((<%=y%>)+1))-log(gamma((<%=x%>)-(<%=y%>)+1)))");

rxPrintOp <- function(op){
    force(op)
    function(...){
        txt <- do.call(op, list(...));
        eval(bquote(rxToSymPy(.(txt))))
    }
}
## rxSympyFEnv$sprintf <- rxPrintOp("sprintf")
## rxSympyFEnv$paste <- rxPrintOp("paste")
## rxSympyFEnv$paste0 <- rxPrintOp("paste0")

sympy.equiv.f <- c("abs", "acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
                   "cos", "cosh", "digamma", "erf", "erfc", "exp", "factorial",
                   "gamma", "log", "log10", "sin", "sinh", "sqrt", "tan",
                   "tanh", "trigamma");

for (f in sympy.equiv.f){
    rxSympyFEnv[[f]] <- functionOp(f);
    sympyRxFEnv[[f]] <- functionOp(f);
}

rxSympyFEnv$psigamma <- function(z, n){
    paste0("polygamma(", n, ", ", z, ")");
}
sympyRxFEnv$polygamma <- function(n, z){
    paste0("psigamma(", z, ", ", n, ")");
}

sympyTransit4 <- function(t, n, mtt, bio, podo="podo"){
    ktr <- paste0("((", n, " + 1) / (", mtt, "))");
    lktr <- paste0("(log((", n, ") + 1) - log(", mtt, "))");
    paste0("exp(log((", bio,") * (", podo, ")) + ", lktr, " + (",
           n, ") * ", "(", lktr," + log(", t, ")) - ",
           ktr," * (", t,") - log(gamma(1 + (", n, "))))");
}

rxSympyFEnv$transit <- function(n, mtt, bio){
    if (missing(bio))
        bio <- 1;
    sympyTransit4("t", n, mtt, bio);
}



allNames <- function(x) {
    if (is.atomic(x)) {
        character()
    } else if (is.name(x)) {
        as.character(x)
    } else if (is.call(x) || is.pairlist(x)) {
        children <- lapply(x[-1], allNames)
        unique(unlist(children))
    } else {
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
}

allCalls <- function(x) {
    if (is.atomic(x) || is.name(x)) {
        character()
    } else if (is.call(x)) {
        fname <- as.character(x[[1]])
        children <- lapply(x[-1], allCalls)
        unique(c(fname, unlist(children)))
    } else if (is.pairlist(x)) {
        unique(unlist(lapply(x[-1], allCalls), use.names = FALSE))
    } else {
        stop("Don't know how to handle type ", typeof(x), call. = FALSE)
    }
}

evalPrints <- function(x, envir=parent.frame()){
    if (is.atomic(x) || is.name(x)) {
        ## Leave unchanged
        return(x);
    } else if (is.call(x)) { # Call recurse_call recursively
        if (identical(x[[1]], quote(sprintf)) ||
            identical(x[[1]], quote(paste)) ||
            identical(x[[1]], quote(paste0))){
            txt <- sprintf("quote(%s)", eval(x, envir));
            txt <- eval(parse(text=txt))
            return(txt)
        } else {
            as.call(lapply(x, evalPrints, envir=envir));
        }
    } else if (is.pairlist(x)) {
        ## Call recurse_call recursively
        as.pairlist(lapply(x, evalPrints, envir=envir));
    } else { # User supplied incorrect input
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
}

unknownSympy <- function(op){
    force(op)
    function(...){
        stop(sprintf("RxODE doesn't know how to translate '%s' to sympy.", op));
    }
}

unknownRx <- function(op){
    force(op)
    function(...){
        stop(sprintf("RxODE doesn't know how to translate '%s' to a RxODE compatible function.", op));
    }
}

cloneEnv <- function(env, parent = parent.env(env)) {
    list2env(as.list(env), parent = parent)
}

sympyEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownSympy), calls)
    callEnv <- list2env(callList);
    rxSympyFEnv <- cloneEnv(rxSympyFEnv, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    n2[n2 == "time"] <- "t";
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSympyFEnv);
}

rxEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownRx), calls)
    callEnv <- list2env(callList);
    rxSympyFEnv <- cloneEnv(sympyRxFEnv, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- gsub(regDfDy, "df(\\1)/dy(\\2)", gsub(regDDt, "d/dt(\\1)", names));
    n2[n2 == "time"] <- "t";
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSympyFEnv);
}


## gamma -> gammafn
## polygamma -> polygamma(n, z) returns log(gamma(z)).diff(n + 1) = psigamma(z, n)
## trigamma -> trigamma

##' Converts model specification to/from a SymPy language
##'
##' @param x is either RxODE family of objects, character of RxODE
##'     expression or unquoted RxODE expression.
##' @return Code Lines for sympy/RxODE that are named by what
##'     variables are defined.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxToSymPy <- function(x, envir=parent.frame(1)) {
    if (class(substitute(x)) == "character"){
        if (length(x) == 1){
            txt <- strsplit(gsub(";", "\n", x), "\n+")[[1]];
            txt <- strsplit(txt, "=", txt);
            vars <- c();
            addNames <- TRUE;
            txt <- unlist(lapply(txt, function(x){
                var <- eval(parse(text=sprintf("rxToSymPy(%s)", x[1])));
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    eq <- eval(parse(text=sprintf("rxToSymPy(%s)", x[2])));
                    return(sprintf("%s = %s", var, eq));
                } else {
                    addNames <<- FALSE
                    return(var);
                }
            }));
            if (addNames){
                names(txt) <- vars;
            }
            return(txt);
        } else {
            return(eval(parse(text=sprintf("rxToSymPy(%s)", deparse(paste(x, collapse="\n"))))));
        }
    } else if (class(substitute(x)) == "name"){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (any(cls == c("list", "rxDll", "RxCompilationManager", "RxODE", "solveRxDll"))){
            ret <- strsplit(rxModelVars(x)$model["normModel"],"\n")[[1]];
            ret <- ret[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),any_spaces,"="), ret) == -1];
            ret <- ret[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),"(0)", any_spaces,"="), ret) == -1];
            ret <- eval(parse(text=sprintf("rxToSymPy(%s,envir=envir)", deparse(paste0(ret, collapse="\n")))), envir=envir);
            return(ret);
        } else if (cls == "character" && length(cls) == 1){
            return(eval(parse(text=sprintf("rxToSymPy(%s)", deparse(x)))));
        } else {
            expr <- evalPrints(substitute(x), envir=envir)
            return(eval(expr, sympyEnv(expr)))
        }
    } else {
        expr <- evalPrints(substitute(x), envir=envir)
        return(eval(expr, sympyEnv(expr)))
    }
}

##' @rdname rxToSymPy
##' @export
rxFromSymPy <- function(x, envir=parent.frame(1)) {
    if (class(substitute(x)) == "character"){
        if (length(x) == 1){
            txt <- strsplit(x, "\n+")[[1]];
            txt <- strsplit(txt, "=", txt);
            vars <- c();
            addNames <- TRUE;
            txt <- unlist(lapply(txt, function(x){
                var <- eval(parse(text=sprintf("rxFromSymPy(%s)", x[1])));
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    eq <- eval(parse(text=sprintf("rxFromSymPy(%s)", x[2])));
                    return(sprintf("%s = %s", var, eq));
                } else {
                    addNames <<- FALSE
                    return(var);
                }
            }));
            if (addNames){
                names(txt) <- vars;
            }
            return(txt);
        } else {
            return(eval(parse(text=sprintf("rxFromSymPy(%s)", deparse(paste(x, collapse="\n"))))));
        }
    } else if (class(substitute(x)) == "name"){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (cls == "character" && length(cls) == 1){
            return(eval(parse(text=sprintf("rxFromSymPy(%s)", deparse(x)))));
        } else {
            expr <- evalPrints(substitute(x), envir=envir)
            return(eval(expr, rxEnv(expr)))
        }
    } else {
        expr <- evalPrints(substitute(x), envir=envir)
        return(eval(expr, rxEnv(expr)))
    }
}
## Stop DSL

rxSymPy.vars <- c();
rxSymPy.sens <- NULL;
rxSymPy.jac <- NULL;
rxSymPy.jac2 <- NULL;
rxSymPy.model <- NULL;
rxSymPy.calcSens <- NULL;

##' Setup sympy variables
##'
##' This creates sympy variables for later evaulation in the CAS sympy.
##'
##' @param model RxODE family of objects
##' @return NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyVars <- function(model){
    if (!exists(".Jython", .GlobalEnv))
        rSymPy::sympyStart();
    if (!class(model) == "character"){
        vars <- c(rxParams(model), rxState(model),
                  "podo", "t", "time", "tlast");
    } else {
        vars <- model
    }
    known <- c(rxSymPy.vars, vars);
    assignInMyNamespace("rxSymPy.vars", known);
    .Jython$exec(sprintf("%s = symbols('%s')", paste(vars, collapse=", "), paste(vars, collapse=" ")));
    return(invisible());
}

##' Setup a sympy environment that sets up the specified RxODE model.
##'
##' @param model RxODE family of objects
##' @return a boolean indicating if the environment is cached.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPySetup <- function(model){
    setup <- rxToSymPy(model)
    cache <- identical(setup, rxSymPy.model);
    if (!cache){
        assignInMyNamespace("rxSymPy.model", setup);
        assignInMyNamespace("rxSymPy.jac", NULL);
        assignInMyNamespace("rxSymPy.calcSens", NULL);
        assignInMyNamespace("rxSymPy.sens", NULL);
        assignInMyNamespace("rxSymPy.jac2", NULL);
        rxSymPyVars(model)
        assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, names(setup)))
        for (line in setup){
            tmp <- line;
            names(tmp) <- NULL;
            .Jython$exec(tmp);
        }
    }
    return(invisible(cache));
}

rxSymPyClearCache <- function(){
    assignInMyNamespace("rxSymPy.model", NULL);
    assignInMyNamespace("rxSymPy.jac", NULL);
    assignInMyNamespace("rxSymPy.calcSens", NULL);
    assignInMyNamespace("rxSymPy.sens", NULL);
    assignInMyNamespace("rxSymPy.jac2", NULL);

}

##' Calculate df/dy derivatives
##'
##' @param model RxODE models
##' @param df is a string for the state in the df(.)/dy(.).  If
##'     missing and dy is missing, all the df(.)/dy(.) components are
##'     calulated according to the \code{vars} parameter below.
##' @param dy is a string for the state or varaible in the df(.)/dy(.).
##' @param vars is a boolean indicating if parameters will be included
##'     for the dy component in the df(.)/dy(.), instead of just state
##'     variables (required for sensitivity equations).
##' @return RxODE syntax lines for the df(.)/dy(.)
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxSymPyDfDy <- function(model, df, dy, vars=FALSE){
    if (missing(df) && missing(dy)){
        if (class(vars) == "logical"){
            if (vars){
                jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), rxParams(model)));
            } else  {
                jac <- expand.grid(s1=rxState(model), s2=rxState(model));
            }
        }
        jac <- with(jac, data.frame(jac,
                          rx=sprintf("df(%s)/dy(%s)", s1, s2),
                          sym=sprintf("__d_df_%s_dy_%s__", s1, s2)))
        cache <- rxSymPySetup(model);
        cache <- cache && !is.null(rxSymPy.jac);
        if (!cache){
            extraLines <- c();
            cat("Calculate Jacobian.");
            for (dfdy in jac$rx){
                if (!any(dfdy == rxDfdy(model))){
                    extraLines[length(extraLines) + 1] <- with(jac[jac$rx == dfdy, ], rxSymPyDfDy(NULL, s1, s2));
                    cat(".")
                }
            }
            cat("done!\n");
            assignInMyNamespace("rxSymPy.jac", extraLines);
        } else {
            extraLines <- rxSymPy.jac;
        }
        return(extraLines);
    } else {
        if (!is.null(model)){
            rxSymPySetup(model);
        }
        var1 <- rxToSymPy(sprintf("d/dt(%s)", df));
        var <- rxToSymPy(sprintf("df(%s)/dy(%s)", df, dy));
        if (grepl(rex::rex(dy), rSymPy::sympy(var1))){
            line <- rSymPy::sympy(sprintf("diff(%s,%s)", var1, dy));
            .Jython$exec(sprintf("%s=%s", var, line));
            ret <- sprintf("df(%s)/dy(%s) = %s", df, dy, rxFromSymPy(line));
        } else {
            .Jython$exec(sprintf("%s=0", var));
            ret <- sprintf("df(%s)/dy(%s) = 0;", df, dy);
        }
        return(ret);
    }
}
##' Calculate the full jacobain for a model
##'
##' This expand the model to caluclate the jacobian.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return RxODE syntax for model with jacobian specified.
##' @author Matthew L. Fidler
rxSymPyJacobian <- function(model){
    extraLines <- rxSymPyDfDy(model, vars=FALSE);
    extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, ";"), extraLines) == -1];
    model <- sprintf("%s\n%s", rxModelVars(model)$model["normModel"], paste(extraLines, collapse="\n"));
    rxSymPyClean()
    return(model);
}
##' Calculate the sensitivity equations for a model
##'
##' This expands the model to calculate sensitivities.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return Model syntax that includes the sensitivity parameters.
##' @author Matthew L. Fidler
##' @export
rxSymPySensitivity <- function(model, calcSens, calcJac=FALSE){
    cache <- rxSymPySetup(model);
    if (missing(calcSens)){
        calcSens <- rxParams(model);
    }
    if (class(calcSens) == "logical"){
        if (calcSens){
            calcSens <- rxParams(model);
        } else {
            stop("It is pointless to request a sensitivity calculation when calcSens=FALSE.")
        }
    }
    state <- rxState(model)
    extraLines <- rxSymPyDfDy(model, vars=TRUE);
    all.sens <- c();
    cache <- cache && identical(calcSens, rxSymPy.calcSens) && !is.null(rxSymPy.sens);
    if (!cache){
        cat("Calculate Sensitivities.");
        for (s1 in state){
            for (sns in calcSens){
                tmp <- c()
                vars <- c();
                for (s2 in state){
                    vars <- c(vars, sprintf("rx__sens_%s_%s__", s2, sns));
                    extra <- sprintf("df(%s)/dy(%s)*rx__sens_%s_%s__", s1, s2, s2, sns)
                    tmp <- c(tmp, extra);
                }
                tmp <- c(tmp, sprintf("df(%s)/dy(%s)", s1, sns));
                rxSymPyVars(vars);
                all.sens <- c(all.sens, vars);
                line <- rSymPy::sympy(rxToSymPy(paste(tmp, collapse=" + ")));
                var.rx <- sprintf("d/dt(rx__sens_%s_%s__)", s1, sns)
                var <- rxToSymPy(var.rx)
                .Jython$exec(sprintf("%s=%s", var, line));
                assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
                extraLines[length(extraLines) + 1] <- sprintf("%s=%s", var.rx, rxFromSymPy(line));
                cat(".")
            }
        }
        cat("\ndone!\n");
        assignInMyNamespace("rxSymPy.calcSens", calcSens);
        assignInMyNamespace("rxSymPy.sens", extraLines);
    } else {
        extraLines <- rxSymPy.sens;
    }
    if (calcJac){
        cache <- cache && !is.null(rxSymPy.jac2);
        if (!cache){
            cat("Expanding Jacobian for sensitivities.")
            jac2 <- expand.grid(s1=unique(all.sens), s2=unique(c(all.sens, rxState(model))))
            for (i in 1:length(jac2$s1)){
                extraLines[length(extraLines) + 1] <- rxSymPyDfDy(NULL, jac2$s1[i], jac2$s2[i]);
                cat(".");
                if (i %% 5 == 0){
                    cat(i);
                }
                if (i %% 50 == 0){
                    cat("\n");
                }
            }
            cat("\ndone!\n");
            extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, ";"), extraLines) == -1];
            assignInMyNamespace("rxSymPy.jac2", extraLines);
        } else {
            extraLines <- rxSymPy.jac2
        }
    } else {
        extraLines <- extraLines[regexpr(rex::rex(regJac, any_spaces, "="), extraLines) == -1];
    }
    extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, ";"), extraLines) == -1];
    ret <- sprintf("%s\n%s", rxModelVars(model)$model["normModel"], paste(extraLines, collapse="\n"));
    rxSymPyClean()
    return(ret);
}

##' Remove variables created by RxODE from the sympy environment.
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyClean <- function(){
    for (v in unique(rxSymPy.vars)){
        try({.Jython$exec(sprintf("del %s", v))}, silent = TRUE);
    }
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
}

## Supported Sympy special functions
## besseli -> besseli(nu,z) -> bessel_i(z,nu,1)
## besselj -> besselj(nu,z) -> bessel_j(z,nu)
## besselk -> besselk(nu,z) -> bessel_k(z,nu,1)
## bessely -> bessely(nu,z) -> bessel_y(z,nu)
## beta -> beta
## digamma -> digamma
## erfc -> erfc
## gamma -> gammafn
## polygamma -> polygamma(n, z) returns log(gamma(z)).diff(n + 1) = psigamma(z, n)
## trigamma -> trigamma
