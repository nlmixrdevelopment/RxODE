.regIdentifier1 <- rex::rex(one_of("a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."))
.regIdentifier2 <- rex::rex(at_least(".",1), one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."));
.regIdentifier <- rex::rex(or(.regIdentifier1, .regIdentifier2));
regSens <- rex::rex("rx__sens_", capture(.regIdentifier), "_BY_",  capture(.regIdentifier), "__");
regSensEtaTheta <- rex::rex("rx__sens_", capture(.regIdentifier), "_BY_",  capture(.regIdentifier),
                            "_BY_",capture(.regIdentifier), "__");
regToSens1 <- rex::rex( capture(.regIdentifier), or("_", ".", ":"),  capture(.regIdentifier));
regToSens2 <- rex::rex( "d/dt(d(", capture(.regIdentifier), ")/d(",  capture(.regIdentifier), "))");
regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                         group(any_of("0":"9"), ".", some_of("0":"9"))),
                      maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9"))));
regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"));
regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
regNum <- rex::rex(maybe("-"), or(regDecimalint, regFloat1, regFloat2))
regDDt <- rex::rex(start, "rx__d_dt_", capture(anything), "__", end);
regDfDy <- rex::rex(start, "rx__df_", capture(anything), "_dy_", capture(anything), "__", end);
regThEt <- rex::rex(capture(or("TH", ""), "ETA"), "_",
                    capture("1":"9", any_of("0":"9")), "_")
regDfDyTh <- rex::rex(start, "rx__df_", capture(anything), "_dy_", regThEt, "__", end);
regEta <- rex::rex(start, "ETA[", capture("1":"9", any_of("0":"9")), "]")
regTheta <- rex::rex(start, "THETA[", capture("1":"9", any_of("0":"9")), "]")
regJac <- rex::rex( "df(", capture(.regIdentifier), ")/dy(",  capture(or(.regIdentifier, group(or("THETA[", "ETA["), "1":"9", any_of("0":"9"), "]"))), ")");
.regRate <- rex::rex(start, "rx_rate_",capture(anything),"_");
.regDur <- rex::rex(start, "rx_dur_",capture(anything),"_");
.regLag <- rex::rex(start, "rx_lag_",capture(anything),"_");
.regF <- rex::rex(start, "rx_f_",capture(anything),"_");
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
regIni <- rex::rex(or(group(one_of("_."), "0"), "0", "(0)", "[0]", "{0}"), end);


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
        x <- .rxRmIni(x);
    }
    if (removePrint){
        x <- .rxRmPrint(x);
    }
    model <- x;
    w1 <- which(regexpr(regIfOrElse, model) != -1);
    w2 <- which(regexpr(regEnd, model) != -1);
    if (length(w1) > 0){
        curr.expr <- c("");
        lst <- list();
        last <- "";
        known <- list();
        for (i in seq_along(model)){
            if (any(i == w1)){
                if (regexpr(regElse, model[i]) != -1){
                    curr.expr[length(curr.expr) + 1] <- sprintf("!(%s)", last);
                } else {
                    curr.expr[length(curr.expr) + 1] <- gsub(regIf, "!(\\1)", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"),
                                                          collapse=" && "), curr.expr[-1]);
                    curr.expr[length(curr.expr)] <- gsub(regIf, "\\1", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"),
                                                          collapse=" && "), curr.expr[-1]);
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
        for (i in seq_along(known)){
            mod <- c();
            for (j in seq_along(model)){
                if (identical(lst[[j]], c(""))){
                    mod[length(mod) + 1] <- model[j];
                } else {
                    i1 <- lst[[j]][-1];
                    i2 <- known[[i]][-1];
                    i3 <- i2[seq(1,min(length(i1), length(i2)))];
                    if (!identical(i2, i3)){
                        ## Find the expression to remove
                        for (k in seq_along(known)){
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

##' Remove INIs
##' @param x RxODE list of lines to remove
##' @return RxODE lines with inis removed.
##' @author Matthew L. Fidler
##' @keywords internal
.rxRmIni <- function(x){
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),any_spaces,or("=", "~")), x) == -1];
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),"(0)", any_spaces,or("=", "~")), x) == -1];
    return(x);
}

##' Remove print statements
##' @param x RxODE lines to remove
##' @return RxODE with print lines removed.
##' @author Matthew L. Fidler
.rxRmPrint <- function(x){
    return(x[regexpr(getFromNamespace("regPrint", "RxODE"), x) == -1]);
}

##' Add a return statment to a function.
##'
##' @param fn Function to deparse
##' @param ret boolean stating if a return statement will be added.
##' @return Function with parens removed and add a return statment.
##' @author Matthew L. Fidler
rxAddReturn <- function(fn, ret=TRUE){
    txt <- deparse(body(fn));
    if (txt[1] == "{"){
        txt <- txt[-c(1, length(txt))];
    }
    ## FIXME, naieve assumption about functions.
    if (ret){
        if (regexpr(rex::rex(or(boundary, start), "return(", anything, ")"), txt[length(txt)], perl=TRUE) == -1){
            ## Add return statement
            if (regexpr(rex::rex(or("=", "~", "<-", "}")), txt[length(txt)]) == -1){
                txt[length(txt)] <- gsub(rex::rex(start, any_spaces, capture(anything), or(";", ""), any_spaces, end), "return(\\1);", txt[length(txt)]);
            }
        }
    }
    return(paste(txt, collapse="\n"))
}



## Start DSL based on http://adv-r.had.co.nz/dsl.html
## These operators are called to create the language and are not called in tests.
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
        } else if (is(e1,"numeric")){
            paste0("S(", e1, ")/", e2)
        } else {
            paste0(e1, "/", e2)
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

dsl.strip.paren <- function(x){
    strip.it <- function(x){
        if (is.call(x)){
            if (length(x) == 1){
                return(x)
            } else if (identical(x[[1]], quote(`(`))){
                return(strip.it(x[[2]]));
            } else {
                return(x)
            }
        } else {
            return(x);
        }
    }
    ret <- eval(parse(text=sprintf("strip.it(quote(%s))", x)))
    ret <- paste(deparse(ret), collapse=" ")
    return(ret)
}

dsl.to.pow <- function(a, b){
    a <- dsl.strip.paren(a);
    b <- dsl.strip.paren(b);
    num <- suppressWarnings({as.numeric(b)});
    if (is.na(num)){
        return(sprintf("Rx_pow(%s, %s)", a, b));
    } else if (abs(num - round(num)) < .Machine$double.eps^0.5){
        return(sprintf("Rx_pow_di(%s, %s)", a, b));
    } else if (abs(num - 0.5) < .Machine$double.eps^0.5){
        return(sprintf("sqrt(%s)", a));
    } else {
        return(sprintf("Rx_pow(%s, %s)", a, b));
    }
}

## Add sympy->C mini DSL for omega parsing

symengineC <- new.env(parent = emptyenv())
symengineC$"**" <- dsl.to.pow
symengineC$"^" <- dsl.to.pow

symengineC$S <- function(x){
    sprintf("%s", x);
}

for (f in c("acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
            "cos", "cosh", "digamma", "erf", "erfc", "exp", "factorial",
            "gamma", "sin", "sinh", "sqrt", "tan",
            "tanh", "trigamma", "rxTBS", "rxTBSd")){
    symengineC[[f]] <- functionOp(f);
}
symengineC$"(" <- unaryOp("(", ")")
for (op in c("+", "-", "*")){
    symengineC[[op]] <- binaryOp(paste0(" ", op, " "));
}

symengineC[["/"]] <- function(e1, e2){
    sprintf("%s /( (%s == 0) ? %s : %s)", e1, e2, .Machine$double.eps, e2)
}


unknownCsymengine <- function(op){
    force(op)
    function(...){
        stop(sprintf("RxODE doesn't support '%s' translation for Omega translation.", op));
    }
}

symengineCEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownCsymengine), calls)
    callEnv <- list2env(callList);
    rxSymPyFEnv <- cloneEnv(symengineC, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    n2 <- gsub(rex::rex("t", capture(numbers)), "REAL(theta)[\\1]", n2)
    n2 <- gsub(rex::rex("pi"), "M_PI", n2)
    n2 <- gsub(rex::rex("rx_SymPy_Res_"), "", n2)
    n2 <- gsub("None", "NA_REAL", n2);
    w <- n2[n2 == "t"];
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=symengineC);
    return(symbol.env)
}

seC <- function(x){
    expr <-eval(parse(text=sprintf("quote(%s)", as.character(x))))
    .ret <- eval(expr, symengineCEnv(expr))
}

sympyTransit4 <- function(t, n, mtt, bio, podo="podo", tlast="tlast"){
    ktr <- paste0("((", n, " + 1)/(", mtt, "))");
    lktr <- paste0("(log((", n, ") + 1) - log(", mtt, "))");
    tc <- paste0("((", t, ")-(", tlast, "))");
    paste0("exp(log((", bio,") * (", podo, ")) + ", lktr, " + (",
           n, ") * ", "(", lktr," + log(", t, ")) - ",
           ktr," * (", t,") - log(gamma(1 + (", n, "))))");
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
            if (regexpr("[=~]", txt) != -1){
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

unknownSymPy <- function(op){
    force(op)
    function(...){
        if (identical(c(...), c(0))){
            return(sprintf("rx_%s_ini_0__", op));
        } else {
            stop(sprintf("RxODE doesn't know how to translate '%s' to SymPy.", op));
        }
    }
}

unknownRx <- function(op){
    force(op)
    function(...){
        if (any(op == c("Derivative", "D", "diff"))){
            .lst <- list(...)
            if (length(.lst) == 2){
                .w <- .lst[[1]];
                .a <- .lst[[2]];
                ## diff(rxTBS(a,lambda,yj),a)
                .reg <- rex::rex("rxTBS(", any_spaces, capture(.a), any_spaces, ",",
                                 capture(except_some_of(",)")), ",", capture(except_some_of(",)")), ")")
                if (regexpr(.reg, .w, perl=TRUE) != -1){
                    .ret <- sub(.reg, "rxTBSd(\\1,\\2,\\3)", .w);
                    return(.ret);
                }
                ## diff(rxTBS2d(a,lambda,yj),a)
                .reg <- rex::rex("rxTBSd(", any_spaces, capture(.a), any_spaces, ",",
                                 capture(except_some_of(",)")), ",", capture(except_some_of(",)")), ")")
                if (regexpr(.reg, .w, perl=TRUE) != -1){
                    return(sub(.reg, "rxTBSd2(\\1,\\2,\\3)", .w));
                }
            }
        }
        stop(sprintf("RxODE doesn't know how to translate '%s' to a RxODE compatible function.", op));
    }
}

cloneEnv <- function(env, parent = parent.env(env)) {
    list2env(as.list(env), parent = parent)
}

sympyEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownSymPy), calls)
    callEnv <- list2env(callList);
    rxSymPyFEnv <- cloneEnv(rxSymPyFEnv, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    n2[n2 == "time"] <- "t";
    ## Replace f with rx_pred_
    n2[n2 == "f"] <- "rx_pred_f_"
    n2 <- gsub(rex::rex("."), "__DoT__", n2)
    ## Replace print functions with nothing.
    n2[regexpr(regPrint, n2) != -1] <- "";
    res <- rxSymPyReserved()
    res <- res[res != "pi"];
    w <- which(n2 %in% res);
    n2[w] <- sprintf("rx_SymPy_Res_%s", n2[w]);
    ## n2 <- gsub(rex::rex("rx_underscore_"), "_", n2);
    n2[n2 == "M_E"] <- "E";
    n2[n2 == "M_PI"] <- "pi";
    n2[n2 == "M_PI_2"] <- "pi/2";
    n2[n2 == "M_PI_4"] <- "pi/4";
    n2[n2 == "M_1_PI"] <- "1/pi";
    n2[n2 == "M_2_PI"] <- "2/pi";
    n2[n2 == "M_2PI"] <- "2*pi";
    n2[n2 == "M_SQRT_PI"] <- "sqrt(pi)";
    n2[n2 == "M_2_SQRTPI"] <- "2/sqrt(pi)";
    n2[n2 == "M_1_SQRT_2PI"] <- "1/sqrt(2*pi)";
    n2[n2 == "M_SQRT_2"] <- "sqrt(2)";
    n2[n2 == "M_SQRT_3"] <- "sqrt(3)";
    n2[n2 == "M_SQRT_32"] <- "sqrt(32)";
    n2[n2 == "M_SQRT_2dPI"] <- "sqrt(2/pi)";
    n2[n2 == "M_LN_SQRT_PI"] <- "log(sqrt(pi))";
    n2[n2 == "M_LN_SQRT_2PI"] <- "log(sqrt(2*pi))";
    n2[n2 == "M_LN_SQRT_PId2"] <- "log(sqrt(pi/2))";
    n2[n2 == "M_SQRT2"] <- "sqrt(2)";
    n2[n2 == "M_SQRT3"] <- "sqrt(3)";
    n2[n2 == "M_SQRT32"] <- "sqrt(32)";
    n2[n2 == "M_LOG10_2"] <- "log10(2)";
    n2[n2 == "M_LOG2E"] <- "1/log(2)"
    n2[n2 == "M_LOG10E"] <- "log10(E)"
    n2[n2 == "M_LN2"] <- "log(2)"
    n2[n2 == "M_LN10"] <- "log(10)"
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSymPyFEnv);
    return(symbol.env)
}

rxEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownRx), calls)
    callEnv <- list2env(callList);
    rxSymPyFEnv <- cloneEnv(sympyRxFEnv, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- gsub(.regRate, "rate(\\1)",
               gsub(.regDur,"dur(\\1)",
               gsub(.regLag,"alag(\\1)",
               gsub(.regF, "f(\\1)",
               gsub(regIni0, "\\1(0)",
               gsub(regDfDy, "df(\\1)/dy(\\2)",
               gsub(regDfDyTh, "df(\\1)/dy(\\2[\\3])",
               gsub(regDDt, "d/dt(\\1)",
               gsub(rex::rex(start, regThEt, end), "\\1[\\2]", names)))))))));
    ## n2 <- gsub(regRate, "rate(\\1)", n2);
    n2 <- gsub(rex::rex("rx_SymPy_Res_"), "", n2)
    n2[n2 == "E"] <- "M_E";
    n2[n2 == "pi"] <- "M_PI";
    n2[n2 == "time"] <- "t";
    n2 <- gsub(rex::rex("__DoT__"), ".", n2)
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSymPyFEnv);
}

.exists2 <- function(x, where){
    .nc <- try(nchar(x) < 1000, silent=TRUE);
    if (inherits(.nc, "try-error")) .nc <- FALSE
    if (rxIs(.nc, "logical")) .nc <- FALSE
    if (.nc){
        return(exists(x, where))
    } else {
        return(FALSE);
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
    return(sprintf("{\n%s\n}", paste(unlist(list(...)), collapse="\n")))
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
rxErrEnv.diag.est <- c();
rxErrEnv.ret <- "rx_r_";
rxErrEnv.init <- NULL;
rxErrEnv.lambda <- NULL;
rxErrEnv.yj <- NULL;

rxErrEnvF$lnorm <- function(est){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The lnorm(.) can only be in an error function.")
    }
    if (!is.null(rxErrEnv.lambda)){
        if (rxErrEnv.lambda != "0" && rxErrEnv.yj != "0"){
            stop("The lnorm(.) cannot be used with other data transformations.")
        }
    }
    estN <- suppressWarnings(as.numeric(est));
    if (is.na(estN)){
        ret <- (sprintf("(%s)^2", est))
        assignInMyNamespace("rxErrEnv.lambda", "0");
        assignInMyNamespace("rxErrEnv.yj", "0");
    } else {
        theta <- sprintf("THETA[%s]", rxErrEnv.theta);
        est <- estN;
        theta.est <- theta;
        ret <- (sprintf("(%s)^2", theta.est))
        tmp <- rxErrEnv.diag.est;
        tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- as.numeric(est);
        assignInMyNamespace("rxErrEnv.diag.est", tmp);
        assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
        assignInMyNamespace("rxErrEnv.lambda", "0");
        assignInMyNamespace("rxErrEnv.yj", "0");
    }
    return(ret);
}

rxErrEnvF$dlnorm <- rxErrEnvF$lnorm
rxErrEnvF$logn <- rxErrEnvF$lnorm

rxErrEnvF$tbs <- function(lambda){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The tbs(.) can only be in an error function.")
    }
    if (!is.null(rxErrEnv.lambda)){
        if (rxErrEnv.yj != "0" & rxErrEnv.lambda != "0" & rxErrEnv.lambda != "1"){
            stop("The tbs(.) cannot be used with other data transformations.")
        }
    }
    estN <- suppressWarnings(as.numeric(lambda));
    if (is.na(estN)){
        assignInMyNamespace("rxErrEnv.lambda", lambda);
        assignInMyNamespace("rxErrEnv.yj", "0");
    } else {
        tmp <- rxErrEnv.diag.est;
        tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- estN;
        assignInMyNamespace("rxErrEnv.lambda", sprintf("THETA[%s]", rxErrEnv.theta));
        assignInMyNamespace("rxErrEnv.diag.est", tmp);
        assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
        assignInMyNamespace("rxErrEnv.yj", "0");
    }
    return("0");
}

rxErrEnvF$boxCox <- rxErrEnvF$tbs

rxErrEnvF$tbsYj <- function(lambda){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The tbsYj(.) can only be in an error function.")
    }
    if (!is.null(rxErrEnv.lambda)){
        if (rxErrEnv.yj != "1" & rxErrEnv.lambda != "0" & rxErrEnv.lambda != "1"){
            stop("The tbsYj(.) cannot be used with other data transformations.")
        }
    }
    estN <- suppressWarnings(as.numeric(lambda));
    if (is.na(estN)){
        assignInMyNamespace("rxErrEnv.lambda", lambda);
        assignInMyNamespace("rxErrEnv.yj", "1");
    } else {
        tmp <- rxErrEnv.diag.est;
        tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- estN;
        assignInMyNamespace("rxErrEnv.lambda", sprintf("THETA[%s]", rxErrEnv.theta));
        assignInMyNamespace("rxErrEnv.diag.est", tmp);
        assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
        assignInMyNamespace("rxErrEnv.yj", "1");
    }
    return("0");
}

rxErrEnvF$yeoJohnson <- rxErrEnvF$tbsYj

rxErrEnvF$add <- function(est){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The add(.) can only be in an error function.")
    }
    estN <- suppressWarnings(as.numeric(est));
    if (is.na(estN)){
        ret <- (sprintf("(%s)^2", est))
    } else {
        theta <- sprintf("THETA[%s]", rxErrEnv.theta);
        est <- estN;
        theta.est <- theta;
        ret <- (sprintf("(%s)^2", theta.est))
        tmp <- rxErrEnv.diag.est;
        tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- as.numeric(est);
        assignInMyNamespace("rxErrEnv.diag.est", tmp);
        assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
        assignInMyNamespace("rxErrEnv.lambda", "1");
        assignInMyNamespace("rxErrEnv.yj", "0");
    }
    return(ret);
}

rxErrEnvF$norm <- rxErrEnvF$add
rxErrEnvF$dnorm <- rxErrEnvF$add

rxErrEnvF$"for" <- function(...){stop("'for' is not supported (yet).")}
rxErrEnvF$`return` <- function(est){
    if (rxErrEnv.ret == ""){
        stop("The PK function should not return anything.")
    }
    .extra <- ""
    force(est)
    if (rxErrEnv.ret == "rx_r_"){
        if (is.null(rxErrEnv.lambda)){
            .lambda <- "1"
        } else {
            .lambda <- rxErrEnv.lambda
        }
        if (is.null(rxErrEnv.yj)){
            .yj <- "0"
        } else {
            .yj <- rxErrEnv.yj
        }
        if (.yj == "0" & .lambda == "1"){
            .yj <- "2";
            .lambda <- "1";
        }
        if (.yj == "0" & .lambda == "0"){
            .yj <- "3";
            .lambda <- "0";
        }
        .extra <- sprintf("rx_yj_~%s;\nrx_lambda_~%s;\n", .yj, .lambda);
        assignInMyNamespace("rxErrEnv.yj", NULL);
        assignInMyNamespace("rxErrEnv.lambda", NULL);
    }
    return(sprintf("%s%s=%s;", .extra, rxErrEnv.ret, est));
}

## rxErrEnvF$c <- function(...){
##     print(sprintf("c(%s)",paste(paste0("rxParseErr(",c(...),")"),collapse=",")))
##     eval(parse(text=sprintf("c(%s)",paste(paste0("rxParseErr(",c(...),")"),collapse=","))))
## }

rxErrEnvF$`|`  <- binaryOp(" | ")
rxErrEnvF$`||`  <- binaryOp(" || ")
rxErrEnvF$`&&`  <- binaryOp(" && ")
rxErrEnvF$`<=`  <- binaryOp(" <= ")
rxErrEnvF$`>=`  <- binaryOp(" >= ")
rxErrEnvF$`==`  <- binaryOp(" == ")

rxErrEnvF$prop <- function(est){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The prop(.) can only be in an error function.")
    }
    estN <- suppressWarnings(as.numeric(est));
    if (is.na(estN)){
        ret <- (sprintf("(rx_pred_f_)^2 * (%s)^2", est))
    } else {
        est <- estN
        ret <- ""
        theta <- sprintf("THETA[%s]", rxErrEnv.theta);
        theta.est <- theta;
        ret <- (sprintf("(rx_pred_f_)^2*(%s)^2", theta.est))
        tmp <- rxErrEnv.diag.est;
        tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- as.numeric(est);
        assignInMyNamespace("rxErrEnv.diag.est", tmp);
        assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 1);
    }
    return(ret);
}

rxErrEnvF$pow <- function(est, pow){
    if (rxErrEnv.ret != "rx_r_"){
        stop("The pow(.) can only be in an error function.")
    }
    estN <- suppressWarnings(as.numeric(est));
    if (is.na(estN)){
        ret <- (sprintf("(rx_pred_f_)^(2*%s) * (%s)^2", pow, est))
    } else {
        est <- estN
        ret <- ""
        theta <- sprintf("THETA[%s]", rxErrEnv.theta);
        theta2 <- sprintf("THETA[%s]", rxErrEnv.theta + 1);
        theta.est <- theta;
        ret <- (sprintf("(rx_pred_f_)^(2*%s) * (%s)^2", theta2, theta.est))
        tmp <- rxErrEnv.diag.est;
        tmp[sprintf("THETA[%s]", rxErrEnv.theta)] <- as.numeric(est);
        tmp[sprintf("THETA[%s]", rxErrEnv.theta + 1)] <- as.numeric(pow);
        assignInMyNamespace("rxErrEnv.diag.est", tmp);
        assignInMyNamespace("rxErrEnv.theta", rxErrEnv.theta + 2);
    }
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
    if (any(n2 == "err")){stop("Use return() for errors.")}
    if (any(n2 == "error")){stop("Use return() for errors.")}
    if (any(n2 == "rx_r")){stop("Use return() for errors.")}
    ## n2[n2 == "err"] <- "rx_r_";
    ## n2[n2 == "error"] <- "rx_r_";
    n2[n2 == "f"] <- "rx_pred_f_";
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
rxParsePred <- function(x, init=NULL, err=NULL){
    if (is.null(err)){
        return(rxParseErr(x, ret="rx_pred_", init=init));
    } else {
        .ini <- attr(err, "ini");
        .errs <- rxExpandIfElse(rxGetModel(err));
        .prd <- rxParseErr(x, ret="rx_pred_", init=init)
        .prd <- rxExpandIfElse(rxGetModel(.prd));
        if (length(.prd) > 1){
            if (length(.prd) == length(.errs)){
                .prd <- .prd[names(.errs)];
                if (any(is.na(.prd))){
                    stop("The errors and predictions need to have the same conditions (if/then statements).")
                }
            } else if (length(.errs) != 1){
                stop("Do not know how to handle this error/pred combination")
            }
        }
        .ret <- sapply(seq(1, max(length(.errs), length(.prd))), function(en){
            .e <- .errs[min(length(.errs), en)];
            .p <- .prd[min(length(.prd), en)];
            .reg <- rex::rex("rx_pred_", any_spaces, "=",any_spaces, capture(except_any_of(";\n")), any_of(";\n"))
            if (regexpr(rex::rex("rx_yj_~2;\nrx_lambda_~1;\n"), .e) != -1){
                .p <- gsub(.reg, "rx_pred_f_~\\1;\nrx_pred_ = \\1", .p)
            } else if (regexpr(rex::rex("rx_yj_~3;\nrx_lambda_~0;\n"), .e) != -1){
                .p <- gsub(.reg, "rx_pred_f_~\\1;\nrx_pred_ = log(\\1)", .p);
            } else {
                .p <- gsub(.reg, "rx_pred_f_~\\1;\nrx_pred_ = rxTBS(\\1, rx_lambda_, rx_yj_)", .p);
            }
            return(gsub("rx_r_", sprintf("%s\nrx_r_", .p), .e));
        })
        if (length(.ret) == 1L){
            .ret <- setNames(.ret, NULL)
            attr(.ret, "ini") <- .ini
            return(.ret);
        } else {
            if (length(.errs) >= length(.prd)){
                .n <- names(.errs);
            } else {
                .n <- names(.prd);
            }
            .ord <- order(sapply(.n, nchar));
            .n <- .n[.ord];
            .ret <- .ret[.ord];
            .ret <- paste(sapply(seq_along(.n), function(x){
                if (x > 1){
                    if (.n[x] == sprintf("(!%s)", .n[x - 1])){
                        return(sprintf("else {\n%s\n}", .ret[x]));
                    }
                }
                return(sprintf("if %s {\n%s\n}", .n[x], .ret[x]));
            }), collapse="\n");
            attr(.ret, "ini") <- .ini
            return(.ret);
        }
    }
}
##' Prepare Error function for inclusion in RxODE
##'
##' @param x error function
##' @param base.theta Base theta to start numbering add(.) and prop(.) from.
##' @param ret Intenral return type.  Should not be changed by the user...
##' @param init Initilization vector
##' @return RxODE transformed text
##' @keywords internal
##' @author Matthew L. Fidler
##' @export
rxParseErr <- function(x, base.theta, ret="rx_r_", init=NULL){
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
    if (is(x,"function")){
        x <- rxAddReturn(x, ret != "");
    }
    if (is(substitute(x),"character")){
        ret <- eval(parse(text=sprintf("RxODE:::rxParseErr(quote({%s}))", x)));
        ret <- substring(ret, 3, nchar(ret) - 2)
        if (regexpr("else if", ret) != -1){
            stop("else if expressions not supported (yet).");
        }
        assignInMyNamespace("rxErrEnv.diag.est", c());
        assignInMyNamespace("rxErrEnv.theta", 1)
        assignInMyNamespace("rxErrEnv.ret", "rx_r_");
        assignInMyNamespace("rxErrEnv.init", NULL);
        return(ret)
    } else if (is(substitute(x),"name")){
        ret <- eval(parse(text=sprintf("RxODE:::rxParseErr(%s)", deparse(x))))
        assignInMyNamespace("rxErrEnv.diag.est", c());
        assignInMyNamespace("rxErrEnv.theta", 1)
        assignInMyNamespace("rxErrEnv.ret", "rx_r_");
        assignInMyNamespace("rxErrEnv.init", NULL);
        return(ret);
    } else {
        ret <- c();
        if (is(x,"character")){
            ret <- eval(parse(text=sprintf("RxODE:::rxParseErr(quote({%s}))", paste(x, collapse="\n"))));
            ret <- substring(ret, 3, nchar(ret) - 2);
        } else {
            ret <- eval(x, rxErrEnv(x));
        }
        attr(ret, "ini") = rxErrEnv.diag.est;
        assignInMyNamespace("rxErrEnv.diag.est", c());
        assignInMyNamespace("rxErrEnv.theta", 1)
        assignInMyNamespace("rxErrEnv.ret", "rx_r_");
        assignInMyNamespace("rxErrEnv.init", NULL);
        if (regexpr("else if", ret) != -1){
            stop("else if expressions not supported (yet).");
        }
        return(ret);
    }
}

sumProdEnv <- new.env(parent = emptyenv())

rxSumProdSum <- FALSE
rxSumProdProd <- FALSE


sumProdEnv$"^" <- binaryOp("^")
sumProdEnv$"**" <- binaryOp("^")

sumProdEnv[["*"]] <- function(a, b){
    if (rxSumProdProd){
        a <- dsl.strip.paren(a)
        b <- dsl.strip.paren(b)
        ## log transformed
        sprintf("prod(%s, %s)", sub(rex::rex(start, "prod(", capture(anything), ")", end), "\\1", a), b)
    } else {
        sprintf("%s * %s", a, b)
    }
}

sumProdEnv[["/"]] <- function(a, b){
    if (rxSumProdProd){
        a <- dsl.strip.paren(a)
        b <- dsl.strip.paren(b)
        sprintf("prod(%s, 1/%s)", sub(rex::rex(start, "prod(", capture(anything), ")", end), "\\1", a), b)
    } else {
        sprintf("%s / %s", a, b)
    }
}

sumProdEnv[["+"]] <- function(a, b){
    if (rxSumProdSum){
        a <- dsl.strip.paren(a)
        if (!missing(b)){
            b <- dsl.strip.paren(b)
            return(sprintf("sum(%s, %s)", sub(rex::rex(start, "sum(", capture(anything), ")", end), "\\1", a), b))
        } else {
            return(a)
        }
    } else {
        if (!missing(b)){
            b <- dsl.strip.paren(b)
            return(sprintf("%s + %s", a, b))
        } else {
            return(a)
        }
    }
}

sumProdEnv[["-"]] <- function(a, b){
    if (rxSumProdSum){
        a <- dsl.strip.paren(a)
        if (!missing(b)){
            b <- dsl.strip.paren(b)
            sprintf("sum(%s, -%s)", sub(rex::rex(start, "sum(", capture(anything), ")", end), "\\1", a), b)
        } else {
            paste0("-", a)
        }
    } else {
        if (!missing(b)){
            b <- dsl.strip.paren(b)
            return(sprintf("%s - %s", a, b))
        } else {
            return(paste0("-", a));
        }
    }
}

sumProdEnv[["("]] <- function(a){
    return(sprintf("%s", a));
}

sumProdEnv$"[" <- function(name, val){
    n <- toupper(name)
    err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(n == c("THETA", "ETA")) && is.numeric(val)){
        if (round(val) == val && val > 0){
            return(sprintf("%s[%s]", n, val));
        } else {
            stop(err);
        }
    } else {
        stop(err)
    }
}

sumProdRxEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, functionOp), calls)
    callEnv <- list2env(callList);
    currEnv <- cloneEnv(sumProdEnv, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=currEnv);
    return(symbol.env)
}

rxSumProd <- function(x){
    return(eval(x, sumProdRxEnv(x)))
}
##' Recast model in terms of sum/prod
##'
##' @param model RxODE model
##' @param expand Boolean indicating if the expression is expanded.
##' @param sum Use sum(...)
##' @param prod Use prod(...)
##' @return model string with prod(.) and sum(.) for all these
##'     operations.
##' @author Matthew L. Fidler
##' @export
rxSumProdModel <- function(model, expand=FALSE, sum=TRUE, prod=TRUE){
    ## Sum for pairwise is equivalent to regular sum under 8 elements.
    assignInMyNamespace("rxSumProdSum", sum)
    assignInMyNamespace("rxSumProdProd", prod)
    lines <- strsplit(rxNorm(model), "\n")[[1]];
    for (i in seq_along(lines)){
        if (regexpr("[=~]", lines[i]) != -1){
            type <- sub(".*([=~]).*", "\\1", lines[i]);
            l0 <- strsplit(lines[i], "[=~]")[[1]];
            l2 <- substr(l0[2], 1, nchar(l0[2]) - 1);
            if (expand){
                l2 <- rxToSE(l2)
                l2 <- symengine::S(l2);
                l2 <- rxFromSE(l2)
            }
            l0[2] <- eval(parse(text=sprintf("rxSumProd(quote(%s))", l2)))
            lines[i] <- paste0(paste0(l0[1], type, l0[2]));
        }
    }
    mod <- paste(lines, collapse="\n")
    return(mod);
}

## rxSplitPlusQ(quote(2*THETA[3]^2*centr*rx__sens_centr_BY_ETA_2__BY_THETA_3___*exp(-2*ETA[2]-2*THETA[2])-
##                    4*THETA[3]^2*centr*rx__sens_centr_BY_THETA_3___*exp(-2*ETA[2]-2*THETA[2])+
##                    2*THETA[3]^2*rx__sens_centr_BY_ETA_2___*rx__sens_centr_BY_THETA_3___*exp(-2*ETA[2]-2*THETA[2])-
##                    4*THETA[3]*centr^2*exp(-2*ETA[2]-2*THETA[2])+
##                    4*THETA[3]*centr*rx__sens_centr_BY_ETA_2___*exp(-2*ETA[2]-2*THETA[2])))

## Stop DSL

rm(f);
rm(op);
