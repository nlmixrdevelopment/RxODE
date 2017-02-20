utils::globalVariables(c(".Jython"))

regIdentifier1 <- rex::rex(one_of("a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."))
regIdentifier2 <- rex::rex(at_least(".",1), one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."));
regIdentifier <- rex::rex(or(regIdentifier1, regIdentifier2));
regSens <- rex::rex("rx__sens_", capture(regIdentifier), "_BY_",  capture(regIdentifier), "__");
regSensEtaTheta <- rex::rex("rx__sens_", capture(regIdentifier), "_BY_",  capture(regIdentifier),
                            "_BY_",capture(regIdentifier), "__");
regToSens1 <- rex::rex( capture(regIdentifier), or("_", ".", ":"),  capture(regIdentifier));
regToSens2 <- rex::rex( "d/dt(d(", capture(regIdentifier), ")/d(",  capture(regIdentifier), "))");
regJac <- rex::rex( "df(", capture(regIdentifier), ")/dy(",  capture(regIdentifier), ")");
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

## Based on normalized grammar output by parser which has
## if (expr) {
## }
## else {
## }
## for the parsing of if/else statements
regIf <- rex::rex(start, any_spaces, "if", any_spaces, "(", capture(anything), ")", any_spaces, "{", any_spaces, end);
regElse <- rex::rex(start, any_spaces, "else", any_spaces, "{", any_spaces, end);
regEnd <- rex::rex(start, any_spaces, "}", any_spaces, end);
regIfOrElse <- rex::rex(or(regIf, regElse))

##' Remove INIs
##' @param x RxODE list of lines to remove
##' @return RxODE lines with inis removed.
##' @author Matthew L. Fidler
##' @keywords internal
rxRmIni <- function(x){
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),any_spaces,"="), x) == -1];
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),"(0)", any_spaces,"="), x) == -1];
    return(x);
}
##' Remove Jacobian
##' @param x RxODE list of lines to remove
##' @return RxODE lines with df/dy removed.
##' @author Matthew L. Fidler
##' @keywords internal
rxRmJac <- function(x){
    return(x[regexpr(rex::rex(regJac, any_spaces, "="), x) == -1])
}
##' Remove sensitivity equations
##' @param x RxODE lines to remove
##' @return Lines with d/dt(rx_sens_...._) removed.
##' @author Matthew L. Fidler
rxRmSens <- function(x){
    return(x[regexpr(rex::rex(start, any_spaces, "d/dt(", any_spaces, regSens, any_spaces, ")", any_spaces, "="), x) == -1]);
}
##' Remove print statements
##' @param x RxODE lines to remove
##' @return RxODE with print lines removed.
##' @author Matthew L. Fidler
rxRmPrint <- function(x){
    return(x[regexpr(regPrint, x) == -1]);
}

##' Expand if/else clauses into mutiple different types of lines.
##'
##'
##' @param model Model can be a character, or a RxODE model.  It needs
##'     to have normalized syntax, that is \code{if (...)\{} has to be
##'     on the same line.  The \code{else} statement must be on its
##'     own line with the closing bracket of the \code{if} statement
##'     on the previous line.  This \code{else} statment must also
##'     contain the opening bracket, like the code \code{else \{}
##' @param removeInis A boolean indicating if parameter
##'     initializations should be removed from the model.
##' @param removePrint A boolean indicating if printing statements
##'     should be removed from the model.
##' @return A named character vector. The names of the vector are the
##'     logical conditions, the values are the lines that satisfy the
##'     logical conditions.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxExpandIfElse <- function(model, removeInis=TRUE, removePrint=TRUE){
    ## expand if/else blocks into a list with lines for conditions that are true
    x <- strsplit(rxNorm(model, FALSE), "\n")[[1]];
    if (removeInis){
        x <- rxRmIni(x);
    }
    if (removePrint){
        x <- rxRmPrint(x);
    }
    model <- x;
    w1 <- which(regexpr(regIfOrElse, model) != -1);
    w2 <- which(regexpr(regEnd, model) != -1);
    if (length(w1) > 0){
        curr.expr <- c("");
        lst <- list();
        last <- "";
        known <- list();
        for (i in 1:length(model)){
            if (any(i == w1)){
                if (regexpr(regElse, model[i]) != -1){
                    curr.expr[length(curr.expr) + 1] <- sprintf("!(%s)", last);
                } else {
                    curr.expr[length(curr.expr) + 1] <- gsub(regIf, "!(\\1)", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"), collapse=" && "), curr.expr[-1]);
                    curr.expr[length(curr.expr)] <- gsub(regIf, "\\1", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"), collapse=" && "), curr.expr[-1]);
                }
                lst[[i]] <- "control";
            } else if (any(i == w2)){
                last <- curr.expr[length(curr.expr)];
                curr.expr <- curr.expr[seq(1, length(curr.expr) - 1)];
                lst[[i]] <- "control";
            } else {
                lst[[i]] <- curr.expr;
            }
        }
        ret <- list();
        rm <- c();
        for (i in 1:length(known)){
            mod <- c();
            for (j in 1:length(model)){
                if (identical(lst[[j]], c(""))){
                    mod[length(mod) + 1] <- model[j];
                } else {
                    i1 <- lst[[j]][-1];
                    i2 <- known[[i]][-1];
                    i3 <- i2[1:min(length(i1), length(i2))];
                    if (!identical(i2, i3)){
                        ## Find the expression to remove
                        for (k in 1:length(known)){
                            i4 <- known[[k]][-1];
                            if (identical(i4, i3)){
                                rm <- c(rm, known[[k]][1]);
                            }
                        }
                    }
                    if (identical(i1, i3)){
                        mod[length(mod) + 1] <- model[j];
                    }
                }
            }
            ret[[known[[i]][1]]] <- paste(mod, collapse="\n");
        }
        ret <- unlist(ret);
        ret <- ret[!(names(ret) %in% rm)];
        return(ret);
    } else {
        return(paste(model, collapse="\n"));
    }
}
rxExpandIfElse.slow <- NULL


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
    n2 <- gsub(regDfDy, "df(\\1)/dy(\\2)",
               gsub(regDfDyTh, "df(\\1)/dy(\\2[\\3])",
                    gsub(regDDt, "d/dt(\\1)",
                         gsub(rex::rex(start, regThEt, end), "\\1[\\2]", names))));
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
            ## txt <- txt[txt != ""];
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
            ret <- strsplit(rxNorm(x),"\n")[[1]];
            ret <- rxRmIni(ret);
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
    if (class(model) == "character" && length(model) > 1){
        vars <- model;
    } else {
        vars <- c(rxParams(model), rxState(model),
                  "podo", "t", "time", "tlast");
    }
    vars <- sapply(vars, function(x){return(rxToSymPy(x))});
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
    setup <- rxToSymPy(model);
    rxSymPyVars(model)
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, names(setup)))
    for (line in setup){
        tmp <- line;
        names(tmp) <- NULL;
        .Jython$exec(tmp);
    }
    return(invisible());
}

##' Setup sympy envirnoment if needed
##'
##' @param model RxODE lines to setup
##' @author Matthew L. Fidler
##' @keywords internal
##' @return model lines
##' @export
rxSymPySetupIf <- function(model){
    if (class(model) != "character"){
        return(rxSymPySetup(model));
    } else {
        lastLine <- sub(rex::rex(start, any_spaces, capture(anything), any_spaces, end),
                        "\\1", strsplit(model[length(model)], "=")[[1]][1])
        if (!rxSymPyExists(rxToSymPy(lastLine))){
            model.setup <- paste(model, collapse="\n");
            rxSymPySetup(model.setup)
        }
        return(model)
    }
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
        return(rxSymPySetupIf(rxSymPyDfDyFull(model, vars, rxCondition(model))));
    } else {
        if (!is.null(model)){
            rxSymPySetup(model);
        }
        var1 <- rxToSymPy(sprintf("d/dt(%s)", df));
        var <- rxToSymPy(sprintf("df(%s)/dy(%s)", df, dy));
        line <- sprintf("diff(%s,%s)", var1, rxToSymPy(dy));
        line <- rSymPy::sympy(line);
        .Jython$exec(sprintf("%s=%s", var, line));
        ret <- sprintf("df(%s)/dy(%s) = %s", df, dy, rxFromSymPy(line));
        return(ret);
    }
}

rxSymPyDfDyFull <- function(model, vars, cond){
    if (class(vars) == "logical"){
        if (vars){
            jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), rxParams(model)),
                               stringsAsFactors=FALSE);
        } else  {
            jac <- expand.grid(s1=rxState(model), s2=rxState(model),
                               stringsAsFactors=FALSE);
        }
    }
    jac <- with(jac, data.frame(jac,
                                rx=sprintf("df(%s)/dy(%s)", s1, s2),
                                sym=sprintf("__d_df_%s_dy_%s__", s1, s2),
                                stringsAsFactors=FALSE))
    rxSymPySetup(model);
    extraLines <- c();
    rxCat("Calculate Jacobian.");
    for (dfdy in jac$rx){
        if (!any(dfdy == rxDfdy(model))){
            extraLines[length(extraLines) + 1] <- with(jac[jac$rx == dfdy, ], rxSymPyDfDy(NULL, s1, s2));
            rxCat(".")
        }
    }
    rxCat("done!\n");
    return(extraLines);
}

rxSymPyDfDyFull.slow <- NULL;

##' Calculate the full jacobain for a model
##'
##' This expand the model to caluclate the jacobian.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return RxODE syntax for model with jacobian specified.
##' @author Matthew L. Fidler
rxSymPyJacobian <- function(model){
    cnd <- rxNorm(model, TRUE); ## Get the conditional statements
    extraLines <- c();
    if (is.null(cnd)){
        extraLines <- rxSymPyDfDy(model, vars=FALSE);
    } else {
        extraLines <- c();
        for (i in cnd){
            rxCat("################################################################################\n");
            rxCat(sprintf("## Calculate for %s\n", i));
            rxCat("################################################################################\n");
            rxCondition(model, i);
            extraLines <- c(extraLines,
                            sprintf("if %s {", i),
                            rxSymPyDfDy(model, vars=FALSE),
                            "}")
        }
        rxCondition(model, FALSE);
    }
    ## extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, or(";","")), extraLines) == -1];
    model <- sprintf("%s\n%s", rxNorm(model), paste(extraLines, collapse="\n"));
    rxSymPyClean();
    return(model);
}

rxSymPyJacobian.slow <- NULL

##' Does the varaible exists in the sympy python environment
##'
##' @param var Variable to test if it exists.
##' @return boolean
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyExists <- function(var){
    return(regexpr(rex::rex("'", var, "'"), rSymPy::sympy("dir()")) != -1)
}
##' Delete variable if exists.
##'
##' @param var Variable to delete.
##' @author Matthew L. Fidler
##' @keywords internal
rxSymPyClear <- function(var){
    if (rxSymPyExists(var)){
        .Jython$exec(sprintf("del %s", var));
    }
}

rxSymPySensitivityFull.text <- "Calculate Sensitivities."

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivityFull <- function(state, calcSens, model, cond){
    rxCat(rxSymPySensitivityFull.text);
    all.sens <- extraLines <- c();
    for (s1 in state){
        for (sns in calcSens){
            tmp <- c()
            vars <- c();
            for (s2 in state){
                vars <- c(vars, sprintf("rx__sens_%s_BY_%s__", s2, rxToSymPy(sns)));
                extra <- sprintf("df(%s)/dy(%s)*rx__sens_%s_BY_%s__", s1, rxToSymPy(s2), s2, rxToSymPy(sns))
                tmp <- c(tmp, extra);
            }
            tmp <- c(tmp, sprintf("df(%s)/dy(%s)", s1, rxToSymPy(sns)));
            rxSymPyVars(vars);
            all.sens <- c(all.sens, vars);
            line <- rxToSymPy(paste(tmp, collapse=" + "))
            line <- rSymPy::sympy(line);
            var.rx <- sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns))
            var <- rxToSymPy(var.rx)
            .Jython$exec(sprintf("%s=%s", var, line));
            assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
            extraLines[length(extraLines) + 1] <- sprintf("%s=%s", var.rx, rxFromSymPy(line));
            rxCat(".")
        }
    }
    rxCat("\ndone!\n");
    return(list(all.sens=all.sens, extraLines=extraLines))
}

rxSymPySensitivityFull.slow <- NULL;


rxSymPySensitivity2Full_ <- function(state, s1, eta, sns, all.sens){
    v1 <- rxToSymPy(sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns)));
    tmp <- c(sprintf("diff(%s,%s)",v1, rxToSymPy(eta)));
    vars <- c();
    for (s2 in state){
        extra <- sprintf("diff(%s,%s)*rx__sens_%s_BY_%s__", v1, rxToSymPy(s2), s2, rxToSymPy(eta))
        v2 <- sprintf("rx__sens_%s_BY_%s_BY_%s__", s2, rxToSymPy(eta), rxToSymPy(sns));
        vars <- c(vars, v2);
        tmp <- rxToSymPy(sprintf("df(%s)/dy(%s)", s1, rxToSymPy(s2)))
        extra2 <- sprintf("%s*%s", tmp, v2)
        tmp <- c(tmp, extra, extra2);
    }
    rxSymPyVars(vars);
    all.sens <- c(all.sens, vars);
    line <- paste(tmp, collapse=" + ");
    line <- rSymPy::sympy(line);
    var.rx <- sprintf("d/dt(rx__sens_%s_BY_%s_BY_%s__)", s1, rxToSymPy(eta), rxToSymPy(sns))
    var <- rxToSymPy(var.rx)
    .Jython$exec(sprintf("%s=%s", var, line));
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
    rxCat(".");
    return(list(all.sens=all.sens, line=sprintf("%s=%s", var.rx, rxFromSymPy(line))));
}

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivity2Full <- function(state, etas, thetas, model, cond){
    all.sens <- extraLines <- c();
    rxCat(rxSymPySensitivityFull.text);
    for (s1 in state){
        for (eta in etas){
            if (identical(etas, thetas)){
                tmp <- rxSymPySensitivity2Full_(state, s1, eta, eta, all.sens);
            } else {
                for (sns in thetas){
                    tmp <- rxSymPySensitivity2Full_(state, s1, eta, sns, all.sens);
                }
            }
            extraLines[length(extraLines) + 1] <- tmp$line;
            all.sens <- tmp$all.sens;
        }
    }
    cat("\ndone!\n");
    return(list(all.sens=all.sens, extraLines=extraLines))
}
rxSymPySensitivity2Full.slow <- NULL;

rxSymPySensitivity.single <- function(model, calcSens, calcJac){
    rxSymPySetupIf(model);
    state <- rxState(model);
    extraLines <- rxSymPyDfDy(model, vars=TRUE);
    if (class(calcSens) == "list" && all(c("eta","theta") %in% names(calcSens))){
        eta <- calcSens$eta;
        theta <- calcSens$theta;
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d(state)/d(eta)) .");
        on.exit({assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate Sensitivites.")})
        ## Calculate dx/dn
        tmp <- rxSymPySensitivityFull(state, eta, model, rxCondition(model));
        all.sens <- tmp$all.sens;
        extraLines <- c(extraLines, rxSymPySetupIf(tmp$extraLines));
        ## Calculate dx/dT
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d(state)/d(theta)) .");
        tmp <- rxSymPySensitivityFull(state, theta, model, rxCondition(model));
        all.sens <- c(all.sens, tmp$all.sens);
        extraLines <- c(extraLines, tmp$extraLines);
        ## Calculate d^2x/dx^2
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d^2(state)/d(eta)^2) .");
        tmp <- rxSymPySensitivity2Full(state, eta, eta, model, rxCondition(model));
        all.sens <- c(all.sens, tmp$all.sens);
        extraLines <- c(extraLines, tmp$extraLines);
        ## Calculate d^2x/(dn dT)
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d^2(state)/d(eta)d(theta)) .");
        tmp <- rxSymPySensitivity2Full(state, eta, theta, model, rxCondition(model));
        all.sens <- c(all.sens, tmp$all.sens);
        extraLines <- c(extraLines, tmp$extraLines);
    } else {
        tmp <- rxSymPySensitivityFull(state, calcSens, model, rxCondition(model))
        extraLines <- c(extraLines, rxSymPySetupIf(tmp$extraLines));
        all.sens <- tmp$all.sens;
    }
    if (calcJac){
        rxCat("Expanding Jacobian for sensitivities.")
        jac2 <- expand.grid(s1=unique(all.sens), s2=unique(c(all.sens, rxState(model))),
                            stringsAsFactors=FALSE)
        for (i in 1:length(jac2$s1)){
            extraLines[length(extraLines) + 1] <- rxSymPyDfDy(NULL, jac2$s1[i], jac2$s2[i]);
            rxCat(".");
            if (i %% 5 == 0){
                rxCat(i);
            }
            if (i %% 50 == 0){
                rxCat("\n");
            }
        }
        rxCat("\ndone!\n");
        ## extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", end), extraLines) == -1];
    } else {
        extraLines <- rxRmJac(extraLines);
    }
    extraLines <- extraLines[regexpr(rex::rex(any_spaces, regJac, any_spaces, "=", any_spaces,
                                              "0", any_spaces, or(";", ""), any_spaces), extraLines) == -1];
    ## cat(paste(extraLines, collapese="\n"), "\n")
    return(extraLines);
}

##' Calculate the sensitivity equations for a model
##'
##' This expands the model to calculate sensitivities.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @param calcSens Either a logical or list of sensitivity parameters
##'     to calculate. When \code{TRUE}, calculate the sensitivities of
##'     all the known parameters.  When \code{FALSE} raise an error.
##' @param calcJac A boolean that determines if the jacobian should be
##'     calculated.
##' @param keepState State parameters to keep the sensitivites for.
##' @return Model syntax that includes the sensitivity parameters.
##' @author Matthew L. Fidler
##' @export
rxSymPySensitivity <- function(model, calcSens, calcJac=FALSE, keepState=NULL){
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
    cnd <- rxNorm(model, TRUE); ## Get the conditional statements
    extraLines <- c();
    if (is.null(cnd)){
        extraLines <- rxSymPySensitivity.single(model, calcSens, calcJac);
    } else {
        extraLines <- c();
        for (i in cnd){
            rxCat("################################################################################\n");
            rxCat(sprintf("## Calculate for %s\n", i));
            rxCat("################################################################################\n");
            rxCondition(model, i);
            extraLines <- c(extraLines,
                            sprintf("if %s {", i),
                            rxSymPySensitivity.single(model, calcSens, calcJac),
                            "}")
        }
        rxCondition(model, FALSE);
    }
    if (!is.null(keepState)){
        w1 <- which(regexpr(rex::rex("d/dt(", regSens), extraLines) != -1);
        if (length(w1) > 0){
            w2 <- which(regexpr(rex::rex("d/dt(rx__sens_", or(keepState)), extraLines[w1]) == -1);
            if (length(w2) > 0){
                extraLines <- extraLines[-w1[w2]];
            }
        }
    }
    ret <- sprintf("%s\n%s", rxNorm(model), paste(extraLines, collapse="\n"));
    return(ret);
}

rxSymPySensitivity.slow <- NULL;

##' Remove variables created by RxODE from the sympy environment.
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyClean <- function(){
    for (v in unique(rxSymPy.vars)){
        rxSymPyClear(v);
    }
    assignInMyNamespace("rxSymPy.vars", c());
}

##' Setup Pred function based on RxODE object
##'
##' @param obj RxODE object
##' @param predfn Prediction function
##' @return RxODE object expanded with predfn and with calculated sensitivities.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPySetupPred <- function(obj, predfn, pkpars, errfn){
    rxSymPyVars(obj);
    on.exit({rxSymPyClean()});
    if (!missing(pkpars)){
        txt <- deparse(body(pkpars));
        if (txt[1] == "{"){
            txt <- txt[-c(1, length(txt))];
        }
        newmod <- rxGetModel(paste0(paste(txt, collapse="\n"), "\n", rxNorm(obj)));
        rxSymPySetupIf(newmod);
        collapseModel <- c()
        for (v in rxState(newmod)){
            tmp <- rSymPy::sympy(rxToSymPy(sprintf("d/dt(%s)", v)));
            tmp <- rxFromSymPy(tmp);
            collapseModel[length(collapseModel) + 1] <- sprintf("d/dt(%s)=%s", v, tmp);
        }
        obj <- rxGetModel(paste(collapseModel, collapse="\n"));
        etas <- rxParams(obj)
        etas <- etas[regexpr(rex::rex(start, "ETA[", "1":"9", any_of("0":"9"), "]"), etas) != -1];
        if (length(etas) > 0){
            calcSens = etas;
        } else {
            calcSens = TRUE;
        }
    } else {
        calcSens = TRUE;
    }
    txt <- deparse(body(predfn));
    if (txt[1] == "{"){
        txt <- txt[-c(1, length(txt))];
    }
    if (regexpr(rex::rex(or(boundary, start), "return(", anything, ")"), txt[length(txt)], perl=TRUE) == -1){
        ## Add return statement
        txt[length(txt)] <- gsub(rex::rex(start, any_spaces, capture(anything), or(";", ""), any_spaces, end), "return(\\1);", txt);
    }
    ## change return(x) to rx_pred = X
    txt <- paste(gsub(rex::rex(or(boundary, start), "return(", capture(anything), ")"), "rx_pred_ = \\1", txt), collapse="\n");
    newmod <- rxGetModel(paste0(rxNorm(obj), "\n", txt), calcSens=calcSens);
    rxSymPySetupIf(newmod);
    ## rxSymPySetup(txt);
    extraLines <- c();
    ## FIXME conditional predfn
    for (state in rxState(newmod)){
        newLine <- rSymPy::sympy(sprintf("diff(rx_pred_,%s)", state));
        newLine <- rxFromSymPy(newLine);
        if (newLine != "0"){
            for (var in calcSens){
                newLine2 <- rSymPy::sympy(sprintf("diff(rx_pred_,%s)", rxToSymPy(var)));
                newLine2 <- rxFromSymPy(newLine2);
                ## (-d(eps)/d(eta)) simialr to Equation 19 in Almquist 2015
                line <- rSymPy::sympy(sprintf("simplify(-(%s))", rxToSymPy(sprintf("(%s)*rx__sens_%s_BY_%s__+(%s)",
                                                                                newLine, state, rxToSymPy(var),
                                                                                newLine2))));
                line <- rxFromSymPy(line);
                ## Chain rule dpred/dstate * dstate/dx = dpred/dx
                extraLines[length(extraLines) + 1] <- sprintf("rx__sens_rx_pred__BY_%s__ = %s", rxToSymPy(var), line);
            }
        }
    }
    newmod <- rxGetModel(paste(c(rxNorm(newmod), extraLines), collapse="\n"))
    if (length(extraLines) == 0){
        stop("Your prediction function does not depend on any of the state variables.")
    }
    if (!missing(errfn)){
        txt <- deparse(body(errfn));
        if (txt[1] == "{"){
            txt <- txt[-c(1, length(txt))];
        }
        if (regexpr(rex::rex(or(boundary, start), "return(", anything, ")"), txt[length(txt)], perl=TRUE) == -1){
            ## Add return statement
            txt[length(txt)] <- gsub(rex::rex(start, any_spaces, capture(anything), or(";", ""), any_spaces, end), "return(\\1);", txt);
        }
        ## change return(x) to rx_r_ = X
        txt <- paste(gsub(rex::rex(or(boundary, start), "return(", capture(anything), ")"), "rx_r_ = (\\1)^2", txt), collapse="\n");
        txt <- rxToSymPy(txt);
        txt <- rxFromSymPy(txt);
        newmod <- rxGetModel(paste0(rxNorm(newmod), "\n", txt));
        rxSymPySetupIf(newmod);
        ## FIXME conditional errfn
        for (state in rxState(newmod)){
            newLine <- rSymPy::sympy(sprintf("diff(rx_r_,%s)", state));
            newLine <- rxFromSymPy(newLine);
            if (newLine != "0"){
                for (var in calcSens){
                    newLine2 <- rSymPy::sympy(sprintf("diff(rx_r_,%s)", rxToSymPy(var)));
                    newLine2 <- rxFromSymPy(newLine2);
                    ## (-d(eps)/d(eta)) simialr to Equation 19 in Almquist 2015
                    line <- rSymPy::sympy(sprintf("simplify(%s)", rxToSymPy(sprintf("(%s)*rx__sens_%s_BY_%s__+(%s)",
                                                                                       newLine, state, rxToSymPy(var),
                                                                                       newLine2))));
                    line <- rxFromSymPy(line);
                    ## Chain rule dpred/dstate * dstate/dx = dpred/dx
                    extraLines[length(extraLines) + 1] <- sprintf("rx__sens_rx_r__BY_%s__ = %s", rxToSymPy(var), line);
                }
            }
        }
        newmod <- rxGetModel(paste(c(rxNorm(newmod), extraLines), collapse="\n"))
        if (length(extraLines) == 0){
            stop("Your error function does not depend on any of the state variables.")
        }
    }
    ## txt <- paste0(rxNorm(obj), "\n", txt);
    return(RxODE(rxNorm(newmod)));
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
