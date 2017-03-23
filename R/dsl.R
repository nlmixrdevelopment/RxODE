regIdentifier1 <- rex::rex(one_of("a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."))
regIdentifier2 <- rex::rex(at_least(".",1), one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."));
regIdentifier <- rex::rex(or(regIdentifier1, regIdentifier2));
regSens <- rex::rex("rx__sens_", capture(regIdentifier), "_BY_",  capture(regIdentifier), "__");
regSensEtaTheta <- rex::rex("rx__sens_", capture(regIdentifier), "_BY_",  capture(regIdentifier),
                            "_BY_",capture(regIdentifier), "__");
regToSens1 <- rex::rex( capture(regIdentifier), or("_", ".", ":"),  capture(regIdentifier));
regToSens2 <- rex::rex( "d/dt(d(", capture(regIdentifier), ")/d(",  capture(regIdentifier), "))");
regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                         group(any_of("0":"9"), ".", some_of("0":"9"))),
                      at_most(group(or("E", "e"), at_most(or("+", "-"), 1), some_of("0":"9")), 1));
regFloat2 <- rex::rex(some_of("0":"9"), or("E", "e"), at_most(or("-", "+"), 1), some_of("0":"9"));
regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
regNum <- rex::rex(at_most("-", 1), or(regDecimalint, regFloat1, regFloat2))
regDDt <- rex::rex(start, "rx__d_dt_", capture(anything), "__", end);
regDfDy <- rex::rex(start, "rx__df_", capture(anything), "_dy_", capture(anything), "__", end);
regThEt <- rex::rex(capture(or("TH", ""), "ETA"), "_",
                    capture("1":"9", any_of("0":"9")), "_")
regDfDyTh <- rex::rex(start, "rx__df_", capture(anything), "_dy_", regThEt, "__", end);
regEta <- rex::rex(start, "ETA[", capture("1":"9", any_of("0":"9")), "]")
regTheta <- rex::rex(start, "THETA[", capture("1":"9", any_of("0":"9")), "]")
regJac <- rex::rex( "df(", capture(regIdentifier), ")/dy(",  capture(or(regIdentifier, group(or("THETA[", "ETA["), "1":"9", any_of("0":"9"), "]"))), ")");
known.print <- c('printf', 'Rprintf', 'print',
                 'jac_printf', 'jac_Rprintf', 'jac_print',
                 'ode_printf', 'ode_Rprintf', 'ode_print',
                 'jac0_printf', 'jac0_Rprintf', 'jac0_print',
                 'ode_printf', 'ode_Rprintf', 'ode_print',
                 'ode0_printf', 'ode0_Rprintf', 'ode0_print',
                 'lhs_printf', 'lhs_Rprintf', 'lhs_print')

regPrint <- rex::rex(start, or(known.print), or(group("(", anything, ")", any_spaces, at_most(";", 1), any_spaces),
                                                group(any_spaces, at_most(";", 1), any_spaces)),
                     end)

regIni0 <- rex::rex(start, "rx_", capture(anything), "_ini_0__", end);


## Start DSL based on http://adv-r.had.co.nz/dsl.html
## These operators are called to create the language and are not called in tests.
## nocov start
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
rxSympyFEnv$c <- function(...){
    eval(parse(text=sprintf("c(%s)",paste(paste0("rxToSymPy(",c(...),")"),collapse=","))))
}
rxSympyFEnv$"/" <- divOp();
rxSympyFEnv$"[" <- function(name, val){
    n <- toupper(name)
    err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(n == c("THETA", "ETA")) && is.numeric(val)){
        if (round(val) == val && val > 0){
            return(sprintf("%s_%s_", n, val));
        } else {
            stop(err);
        }
    } else {
        stop(err)
    }
};

sympyRxFEnv$"/" <- binaryOp(" / ");
rxSympyFEnv$"^" <- binaryOp("**")
rxSympyFEnv$"**" <- binaryOp("**")
sympyRxFEnv$"**" <- binaryOp("^")
sympyRxFEnv$"^" <- binaryOp("^")

functionIgnore <- function(){
    function(...){
        return("")
    }
}

for (p in known.print){
    rxSympyFEnv[[p]] <- functionIgnore();
}

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
sympyRxFEnv$loggamma <- functionOp("lgamma");
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

## equivalent functions
sympy.equiv.f <- c("abs", "acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
                   "cos", "cosh", "digamma", "erf", "erfc", "exp", "factorial",
                   "gamma", "log", "log10", "sin", "sinh", "sqrt", "tan",
                   "tanh", "trigamma")
for (f in sympy.equiv.f){
    rxSympyFEnv[[f]] <- functionOp(f);
    sympyRxFEnv[[f]] <- functionOp(f);
}

rxSympyFEnv$structure <- function(one, ..., .Names){
    eval(parse(text=sprintf("rxToSymPy(%s)", deparse(sprintf("%s", one)))));
}

sympyRxFEnv$structure <- function(one, ..., .Names){
    eval(parse(text=sprintf("rxFromSymPy(%s)", deparse(sprintf("%s", one)))));
}

rxSympyAllowDiff <- FALSE;
rxSympyFEnv$diff <- function(fn, x){
    if (rxSympyAllowDiff){
        sprintf("diff(%s,%s)", fn, x);
    } else {
        stop("diff is not suported in RxODE.");
    }
}

rxSympyFEnv$psigamma <- function(z, n){
    paste0("polygamma(", n, ", ", z, ")");
}
sympyRxFEnv$polygamma <- function(n, z){
    paste0("psigamma(", z, ", ", n, ")");
}

## Add sympy->C mini DSL for omega parsing

rxSympyC <- new.env(parent = emptyenv())
rxSympyC$"**" <- function(a, b){
    sprintf("pow(%s, %s)", a, b);
}
rxSympyC$"^" <- function(a, b){
    sprintf("pow(%s, %s)", a, b);
}
rxSympyC$exp <- functionOp("exp");
rxSympyC$exp <- functionOp("log");
rxSympyC$"(" <- unaryOp("(", ")")
for (op in c("+", "-", "*", "/")){
    rxSympyC[[op]] <- binaryOp(paste0(" ", op, " "));
}

unknownCSympy <- function(op){
    force(op)
    function(...){
        stop(sprintf("RxODE doesn't support '%s' translation for Omega translation.", op));
    }
}

sympyCEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownCSympy), calls)
    callEnv <- list2env(callList);
    rxSympyFEnv <- cloneEnv(rxSympyC, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    n2 <- gsub(rex::rex("t", capture(numbers)), "REAL(theta)[\\1]", n2)
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSympyC);
    return(symbol.env)
}

sympyC <- function(x){
    if (class(substitute(x)) == "character"){
        return(eval(parse(text=sprintf("RxODE:::sympyC(quote(%s))", x))))
    } else if (class(substitute(x)) == "name"){
        return(eval(parse(text=sprintf("RxODE:::sympyC(%s)", deparse(x)))));
    } else {
        return(eval(x, sympyCEnv(x)))
    }
}


## nocov end

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
            identical(x[[1]], quote(paste0)) ||
            identical(x[[1]], quote(sympy)) ||
            identical(x[[1]], quote(rxToSymPy)) ||
            identical(x[[1]], quote(rxFromSymPy))){
            txt <- sprintf("%s", eval(x, envir));
            if (regexpr("=", txt) != -1){
                txt <- deparse(txt);
            } else {
                txt <- paste0("quote(", txt, ")")
            }
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
        if (identical(c(...), c(0))){
            return(sprintf("rx_%s_ini_0__", op));
        } else {
            stop(sprintf("RxODE doesn't know how to translate '%s' to sympy.", op));
        }
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
    ## Replace f with rx_pred_
    n2[n2 == "f"] <- "rx_pred_"
    ## Replace print functions with nothing.
    n2[n2 %in% c('print', 'jac_print', 'ode_print', 'jac0_print', 'ode_print', 'ode0_print', 'lhs_print')] <- "";
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSympyFEnv);
    return(symbol.env)
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
    n2 <- gsub(regIni0, "\\1(0)",
               gsub(regDfDy, "df(\\1)/dy(\\2)",
                    gsub(regDfDyTh, "df(\\1)/dy(\\2[\\3])",
                         gsub(regDDt, "d/dt(\\1)",
                              gsub(rex::rex(start, regThEt, end), "\\1[\\2]", names)))));
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
            names(x) <- NULL;
            txt <- strsplit(gsub(";", "\n", x), "\n+")[[1]];
            txt <- strsplit(txt, "=", txt);
            vars <- c();
            addNames <- TRUE;
            txt <- unlist(lapply(txt, function(x){
                var <- eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", x[1])));
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    eq <- eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", x[2])));
                    return(sprintf("%s = %s", var, eq));
                } else {
                    addNames <<- FALSE
                    return(var);
                }
            }));
            ## txt <- txt[txt != ""];
            if (addNames){
                names(txt) <- vars;
            }
            return(txt);
        } else {
            return(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", deparse(paste(as.vector(x), collapse="\n"))))));
        }
    } else if (class(substitute(x)) == "name"){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (any(cls == c("list", "rxDll", "RxCompilationManager", "RxODE", "solveRxDll"))){
            ret <- strsplit(rxNorm(x),"\n")[[1]];
            ret <- rxRmIni(ret);
            ret <- eval(parse(text=sprintf("RxODE::rxToSymPy(%s,envir=envir)", deparse(paste0(as.vector(ret), collapse="\n")))), envir=envir);
            return(ret);
        } else if (cls == "character" && length(cls) == 1){
            return(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", deparse(as.vector(x))))));
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
                var <- eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", x[1])));
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    eq <- eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", x[2])));
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
            return(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", deparse(paste(x, collapse="\n"))))));
        }
    } else if (class(substitute(x)) == "name"){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (cls == "character" && length(cls) == 1){
            return(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", deparse(x)))));
        } else {
            expr <- evalPrints(substitute(x), envir=envir)
            return(eval(expr, rxEnv(expr)))
        }
    } else {
        expr <- evalPrints(substitute(x), envir=envir)
        return(eval(expr, rxEnv(expr)))
    }
}

## Start error function DSL

rxErrEnvF <- new.env(parent = emptyenv())
for (op in c("+", "-", "*", "/", "^", "**",
             "!=", "==", "&", "&&", "|", "||")){
    op2 <- op
    if (op == "**"){
        op2 <- "^";
    }
    rxErrEnvF[[op]] <- binaryOp(paste0(" ", op2, " "));
}
for (op in c("=", "~", "<-")){
    rxErrEnvF[[op]] <- binaryOp(paste0(" = "));
}
rxErrEnvF$"{" <- function(...){
    return(sprintf("{\n%s;\n}", paste(unlist(list(...)), collapse=";\n")))
}
rxErrEnvF$"(" <- unaryOp("(", ")");
rxErrEnvF$"[" <- function(name, val){
    n <- toupper(name)
    err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(n == c("THETA", "ETA")) && is.numeric(val)){
        if (round(val) == val && val > 0){
            if (n == "THETA" && as.numeric(val) <= length(rxErrEnv.init)){
                return(sprintf("THETA[%s]", val));
            } else {
                return(sprintf("%s[%s]", n, val));
            }
        } else {
            stop(err);
        }
    } else {
        stop(err)
    }
}

rxErrEnvF$"if" <- function(lg, tr, fl){
    if (missing(fl)){
        return(sprintf("if (%s) %s", lg, tr))
    } else {
        return(sprintf("if (%s) %s else %s", lg, tr, fl))
    }
}
rxErrEnv.theta <- 1;
rxErrEnv.diag.xform <- "sqrt";
rxErrEnv.diag.est <- c();
rxErrEnv.ret <- "rx_r_";
rxErrEnv.init <- NULL;

rxErrEnvF$add <- function(est){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The add(.) can only be in an error function.")
    }
    theta <- sprintf("THETA[%s]", rxErrEnv.theta);
    est <- as.numeric(est);
    theta.est <- theta;
    if (rxErrEnv.diag.xform == "sqrt"){
        ret <- (sprintf("(%s)^2", theta.est))
    } else if (rxErrEnv.diag.xform == "log"){
        ret <- (sprintf("exp(%s)", theta.est))
    } else {
        ret <- (sprintf("%s", theta.est));
    }
    tmp <- rxErrEnv.diag.est;
    tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- as.numeric(est);
    assignInMyNamespace("rxErrEnv.diag.est", tmp);
    assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
    return(ret);
}
rxErrEnvF$"for" <- function(...){stop("for not supported (yet)")}
rxErrEnvF$"return" <- function(est){
    if (rxErrEnv.ret == ""){
        stop("The PK function should not return anything.")
    }
    return(sprintf("%s = %s", rxErrEnv.ret, est));
}

## rxErrEnvF$c <- function(...){
##     print(sprintf("c(%s)",paste(paste0("rxParseErr(",c(...),")"),collapse=",")))
##     eval(parse(text=sprintf("c(%s)",paste(paste0("rxParseErr(",c(...),")"),collapse=","))))
## }

rxErrEnvF$prop <- function(est){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The prop(.) can only be in an error function.")
    }
    ret <- ""
    theta <- sprintf("THETA[%s]", rxErrEnv.theta);
    est <- as.numeric(est);
    theta.est <- theta;
    if (rxErrEnv.diag.xform == "sqrt"){
        ret <- (sprintf("rx_pred_^2 * (%s)^2", theta.est))
    } else if (rxErrEnv.diag.xform == "log"){
        ret <- (sprintf("rx_pred_^2 * exp(%s)", theta.est))
    } else {
        ret <- (sprintf("rx_pred_^2 * %s", theta.est))
    }
    tmp <- rxErrEnv.diag.est;
    tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- as.numeric(est);
    assignInMyNamespace("rxErrEnv.diag.est", tmp);
    assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
    return(ret);
}

rxErrEnv <- function(expr){
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, functionOp), calls)
    callEnv <- list2env(callList);

    ## Known functions
    rxErrFEnv <- cloneEnv(rxErrEnvF, callEnv);

    ## Symbols
    names <- allNames(expr)
    n1 <- names;
    n2 <- names;
    n2[n2 == "time"] <- "t";
    n2[n2 == "err"] <- "rx_r_";
    n2[n2 == "error"] <- "rx_r_";
    n2[n2 == "f"] <- "rx_pred_";
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxErrFEnv);
    return(symbol.env)
}

##' Parse PK function for inclusion in RxODE
##'
##' @param x PK function
##' @inheritParams rxParseErr
##' @return RxODE transformed text.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxParsePk <- function(x, init=NULL){
    return(rxParseErr(x, init=init, ret=""));
}
##' Prepare Pred function for inclusion in RxODE
##'
##' @param x pred function
##' @inheritParams rxParseErr
##' @return RxODE transformed text.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxParsePred <- function(x, init=NULL){
    return(rxParseErr(x, ret="rx_pred_", init=init));
}
##' Prepare Error function for inclusion in RxODE
##'
##' @param x error function
##' @param base.theta Base theta to start numbering add(.) and prop(.) from.
##' @param diag.xform Diagonal form of variance parameters
##' @param ret Intenral return type.  Should not be changed by the user...
##' @param init Initilization vector
##' @return RxODE transformed text
##' @keywords internal
##' @author Matthew L. Fidler
##' @export
rxParseErr <- function(x, base.theta, diag.xform=c("sqrt", "log", "identity"),
                       ret="rx_r_", init=NULL){
    if (!missing(diag.xform)){
        diag.xform <- match.arg(diag.xform)
        assignInMyNamespace("rxErrEnv.diag.xform", diag.xform);
    }
    if (!missing(base.theta)){
        assignInMyNamespace("rxErrEnv.theta", base.theta);
    }
    if (!missing(ret)){
        assignInMyNamespace("rxErrEnv.ret", ret);
    }
    if (!missing(init)){
        assignInMyNamespace("rxErrEnv.init", init);
    }
    if (!missing(init)){
        assignInMyNamespace("rxErrEnv.init", init);
    }
    if (class(x) == "function"){
        x <- rxAddReturn(x, ret != "");
    }
    if (class(substitute(x)) == "character"){
        ret <- eval(parse(text=sprintf("RxODE:::rxParseErr(quote({%s}))", x)));
        ret <- substring(ret, 3, nchar(ret) - 2)
        if (regexpr("else if", ret) != -1){
            stop("else if expressions not supported (yet).");
        }
        assignInMyNamespace("rxErrEnv.diag.est", c());
        assignInMyNamespace("rxErrEnv.diag.xform", "sqrt");
        assignInMyNamespace("rxErrEnv.theta", 1)
        assignInMyNamespace("rxErrEnv.ret", "rx_r_");
        assignInMyNamespace("rxErrEnv.init", NULL);
        return(ret)
    } else if (class(substitute(x)) == "name"){
        ret <- eval(parse(text=sprintf("RxODE:::rxParseErr(%s)", deparse(x))))
        if (regexpr("else if", ret) != -1){
            stop("else if expressions not supported (yet).");
        }
        assignInMyNamespace("rxErrEnv.diag.est", c());
        assignInMyNamespace("rxErrEnv.diag.xform", "sqrt");
        assignInMyNamespace("rxErrEnv.theta", 1)
        assignInMyNamespace("rxErrEnv.ret", "rx_r_");
        assignInMyNamespace("rxErrEnv.init", NULL);
        return(ret);
    } else {
        ret <- c();
        if (class(x) == "character"){
            ret <- eval(parse(text=sprintf("RxODE:::rxParseErr(quote({%s}))", paste(x, collapse="\n"))));
            ret <- substring(ret, 3, nchar(ret) - 2);
        } else {
            ret <- eval(x, rxErrEnv(x));
        }
        attr(ret, "ini") = rxErrEnv.diag.est;
        assignInMyNamespace("rxErrEnv.diag.est", c());
        assignInMyNamespace("rxErrEnv.diag.xform", "sqrt");
        assignInMyNamespace("rxErrEnv.theta", 1)
        assignInMyNamespace("rxErrEnv.ret", "rx_r_");
        assignInMyNamespace("rxErrEnv.init", NULL);
        if (regexpr("else if", ret) != -1){
            stop("else if expressions not supported (yet).");
        }
        return(ret);
    }
}

## Stop DSL

