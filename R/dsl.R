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
        } else if (class(e1) == "numeric"){
            paste0("S(", e1, ") / ", e2)
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
rxSymPyFEnv <- new.env(parent = emptyenv())
for (op in c("+", "-", "*")){
    rxSymPyFEnv[[op]] <- binaryOp(paste0(" ", op, " "));
    sympyRxFEnv[[op]] <- binaryOp(paste0(" ", op, " "));
}

rxSymPyFEnv[["~"]] <- binaryOp(" = ");

rxSymPyFEnv$c <- function(...){
    eval(parse(text=sprintf("c(%s)",paste(paste0("rxToSymPy(",c(...),")"),collapse=","))))
}
rxSymPyFEnv$"/" <- divOp();
rxSymPyFEnv$"[" <- function(name, val){
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
rxSymPyFEnv$"^" <- binaryOp("**")
rxSymPyFEnv$"**" <- binaryOp("**")
sympyRxFEnv$"**" <- binaryOp("^")
sympyRxFEnv$"^" <- binaryOp("^")

functionIgnore <- function(){
    function(...){
        return("")
    }
}

for (p in known.print){
    rxSymPyFEnv[[p]] <- functionIgnore();
}

rxSymPyFEnv$dt <- function(e1){
    paste0("__dt__", e1, "__");
}

rxSymPyFEnv$df <- function(e1){
    paste0("__df_", e1, "_")
}

rxSymPyFEnv$dy <- function(e1){
    paste0("_dy_", e1, "__")
}

## rx -> sympy
rxSymPyFEnv$bessel_i <- besselOp("i");
rxSymPyFEnv$bessel_j <- besselOp("j");
rxSymPyFEnv$bessel_k <- besselOp("k");
rxSymPyFEnv$bessel_y <- besselOp("y");

## sympy -> rx
sympyRxFEnv$besseli <- besselOp2("i");
sympyRxFEnv$besselj <- besselOp2("j");
sympyRxFEnv$besselk <- besselOp2("k");
sympyRxFEnv$bessely <- besselOp2("y");


## rxSymPyFEnv$"[" <- binaryOp("_")
##
rxSymPyFEnv$"(" <- unaryOp("(", ")")
sympyRxFEnv$"(" <- unaryOp("(", ")")

## pow -> **
rxSymPyFEnv$pow <- binaryOp2("**");
sympyRxFEnv$pow <- binaryOp2("^");

## gammafn -> gamma
rxSymPyFEnv$gammafn <- functionOp("gamma")
rxSymPyFEnv$lgammafn <- functionOp2("log(gamma(", "))")
rxSymPyFEnv$lgamma <- functionOp2("log(gamma(", "))")
rxSymPyFEnv$tetragamma <- functionOp2("psigamma(", ", 2)")
rxSymPyFEnv$pentagamma <- functionOp2("psigamma(", ", 3)")
rxSymPyFEnv$lbeta <- functionOp2("log(beta(", "))")
rxSymPyFEnv$lgamma1p <- functionOp2("log(gamma((", ")+1))")
rxSymPyFEnv$cospi <- functionOp2("cos(pi * (", "))")
rxSymPyFEnv$sinpi <- functionOp2("sin(pi * (", "))")
rxSymPyFEnv$tanpi <- functionOp2("tan(pi * (", "))")
rxSymPyFEnv$logspace_add <- binaryOp2(" + ");
rxSymPyFEnv$logspace_sub <- binaryOp2(" - ");
sympyRxFEnv$loggamma <- functionOp("lgamma");
## Following R functions are not translated
## ftrunc, fround, fprec, fsign, sign, fmin2, fmax2, imin2, imax2, logspace_sum
## choose lchoose

rxSymPyFEnv$R_pow <- binaryOp2("**");
rxSymPyFEnv$R_pow_di <- binaryOp2("**");
rxSymPyFEnv$log1p <- functionOp2("log(1 + (", "))");
rxSymPyFEnv$log1pmx <- functionBrewx("(log(1 + (<%=x%>))-(<%=x%>))");
rxSymPyFEnv$expm1 <- functionOp2("(exp(", ")-1)");

rxSymPyFEnv$choose <- functionBrewxy("(factorial(<%=x%>)/(factorial(<%=y%>)*factorial((<%=x%>)-(<%=y%>))))");
rxSymPyFEnv$lchoose <- functionBrewxy("(log(gamma((<%=x%>)+1))-log(gamma((<%=y%>)+1))-log(gamma((<%=x%>)-(<%=y%>)+1)))");

rxPrintOp <- function(op){
    force(op)
    function(...){
        txt <- do.call(op, list(...));
        eval(bquote(rxToSymPy(.(txt))))
    }
}
## rxSymPyFEnv$sprintf <- rxPrintOp("sprintf")
## rxSymPyFEnv$paste <- rxPrintOp("paste")
## rxSymPyFEnv$paste0 <- rxPrintOp("paste0")

## equivalent functions
sympy.equiv.f <- c("abs", "acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
                   "cos", "cosh", "digamma", "erf", "erfc", "exp", "factorial",
                   "gamma", "log", "log10", "sin", "sinh", "sqrt", "tan",
                   "tanh", "trigamma")
for (f in sympy.equiv.f){
    rxSymPyFEnv[[f]] <- functionOp(f);
    sympyRxFEnv[[f]] <- functionOp(f);
}

rxSymPyAbsLog <- FALSE
rxSymPyLogSign <- c()

rxSymPyExpThetas <- c()
rxSymPyExpEtas <- c()

sympyRxFEnv$exp <- function(arg){
    ## Check for log-scale parameters.
    theta.reg <- rex::rex(boundary, or("THETA", "theta"), or("[", ")"));
    eta.reg <- rex::rex(boundary, or("ETA", "eta"), or("[", ")"));
    tmp <- suppressWarnings({as.numeric(gsub(rex::rex(start, capture(any_numbers), or(")", "]"), anything),"\\1",
                                             strsplit(paste0(arg), theta.reg, perl=TRUE)[[1]]))});
    tmp <- sort(unique(c(rxSymPyExpThetas, tmp[!is.na(tmp)])));
    assignInMyNamespace("rxSymPyExpThetas", tmp);
    ##
    tmp <- gsub(theta.reg, "", paste0(arg))
    tmp <- suppressWarnings({as.numeric(gsub(rex::rex(start, capture(any_numbers), or(")", "]"), anything),"\\1",
                                             strsplit((tmp), eta.reg, perl=TRUE)[[1]]))});
    tmp <- sort(unique(c(rxSymPyExpEtas, tmp[!is.na(tmp)])));
    assignInMyNamespace("rxSymPyExpEtas", tmp);
    return(paste0("exp(", arg, ")"))
}

sympyRxFEnv$log <- function(arg){
    if (rxSymPyAbsLog){
        tmp <- rxSymPyLogSign
        argn <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", arg)))
        if (length(argn) > 1){
            tmp[length(tmp) + 1] <- paste0("(", arg, ")");
        } else {
            tmp[length(tmp) + 1] <- arg;
        }
        assignInMyNamespace("rxSymPyLogSign", tmp);
        return(paste0("abs_log(", arg, ")"))
    } else {
        return(paste0("log(", arg, ")"))
    }
}

rxSymPyFEnv$structure <- function(one, ..., .Names){
    eval(parse(text=sprintf("rxToSymPy(%s)", deparse(sprintf("%s", one)))));
}

sympyRxFEnv$structure <- function(one, ..., .Names){
    eval(parse(text=sprintf("rxFromSymPy(%s)", deparse(sprintf("%s", one)))));
}

sympyRxFEnv$Subs <- function(expr, what, with){
    what <- strsplit(substring(what, 2, nchar(what) - 1), ",")[[1]];
    with <- strsplit(substring(with, 2, nchar(with) - 1), ",")[[1]];
    for (i in 1:length(what)){
        expr <- gsub(rex::rex(boundary, what[i], boundary), with[i], expr, perl=TRUE);
    }
    return(expr);
}

sympyRxFEnv$subs <- sympyRxFEnv$Subs;

rxSymPyAllowDiff <- FALSE;

rxSymPyDiff <- function(name){
    force(name)
    function(fn, x){
        stop(sprintf("'%s' is not suported in RxODE.", name));
    }
}

rxSymPyFEnv$diff <- rxSymPyDiff("diff");
rxSymPyFEnv$D <- rxSymPyDiff("D");
rxSymPyFEnv$Derivative <- rxSymPyDiff("Derivative");

rxSymPyFEnv$psigamma <- function(z, n){
    paste0("polygamma(", n, ", ", z, ")");
}
sympyRxFEnv$polygamma <- function(n, z){
    paste0("psigamma(", z, ", ", n, ")");
}

## Add sympy->C mini DSL for omega parsing

rxSymPyC <- new.env(parent = emptyenv())
rxSymPyC$"**" <- function(a, b){
    sprintf("pow(%s, %s)", a, b);
}
rxSymPyC$"^" <- function(a, b){
    sprintf("pow(%s, %s)", a, b);
}

rxSymPyC$S <- function(x){
    sprintf("%s", x);
}

for (f in sympy.equiv.f){
    rxSymPyC[[f]] <- functionOp(f);
}
rxSymPyC$"(" <- unaryOp("(", ")")
for (op in c("+", "-", "*", "/")){
    rxSymPyC[[op]] <- binaryOp(paste0(" ", op, " "));
}


unknownCSymPy <- function(op){
    force(op)
    function(...){
        stop(sprintf("RxODE doesn't support '%s' translation for Omega translation.", op));
    }
}

sympyCEnv <- function(expr){
    ## Known functions
    calls <- allCalls(expr)
    callList <- setNames(lapply(calls, unknownCSymPy), calls)
    callEnv <- list2env(callList);
    rxSymPyFEnv <- cloneEnv(rxSymPyC, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    n2 <- gsub(rex::rex("t", capture(numbers)), "REAL(theta)[\\1]", n2)
    n2 <- gsub(rex::rex("pi"), "M_PI", n2)
    n2 <- gsub(rex::rex(start, "rx_SymPy_Res_"), "", n2)
    n2 <- gsub("None", "NA_REAL", n2);
    w <- n2[n2 == "t"];
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSymPyC);
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

rxSymPyFEnv$transit <- function(n, mtt, bio){
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

rxDefinedDerivatives <- new.env(parent = emptyenv())

rxDefinedDerivatives$solveLinB <- function(fn, var){
    fn <- fn[-1];
    diff1 <- fn[3];
    diff2 <- fn[4];
    if (diff2 != "0"){
        stop("Cannot currently take the third order derivitave of solveLinB.")
    } else{
        vals <- paste(fn[-seq(1:4)]);
        w <- which(vals == var)
        if (diff1 == "0"){
            ret <- sprintf("solveLinB(%s,%s,0,%s)", paste(fn[1:2], collapse=","), w, paste(vals, collapse=","))
        } else {
            ret <- sprintf("solveLinB(%s,%s,%s)", paste(fn[1:3], collapse=","), w, paste(vals, collapse=","));
        }
        return(ret);
    }
}

changeDerivs <- function(fn, var, var2=NULL){
    ## Fn is a vector fn[1] == function name,the rest are the argments
    if (length(var) > 1){
        env <- rxEnv(var)
        fnl <- as.list(var[-1])
        if (any(names(sympyRxFEnv) == var[1])){
            fne <- sympyRxFEnv[[var[1]]]
            var <- do.call(fne, fnl)
        } else {
            stop("Cannot figure out how to deparse the deriavative");
        }
    }
    if (!is.null(var2)){
        if (length(var2) > 1){
            env <- rxEnv(var2)
            fnl <- as.list(var2[-1])
            if (any(names(sympyRxFEnv) == var2[1])){
                fne <- sympyRxFEnv[[var2[1]]]
                var2 <- do.call(fne, fnl)
            } else {
                stop("Cannot figure out how to deparse the deriavative");
            }
        }
    }
    if (any(names(rxDefinedDerivatives) == fn[1])){
        fne <- rxDefinedDerivatives[[fn[1]]];
        ret <- do.call(fne, list(fn, var));
        if (!is.null(var2)){
            ## Send through parser recursively...
            ret <- sprintf("D(%s,%s)", ret, var2);
            return(rxFromSymPy(ret));
        }
        return(ret)
    } else {
        stop(sprintf("RxODE does not know how to take a deriavite of '%s'", fn[1]));
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
        } else if ((identical(x[[1]], quote(Derivative)) ||
                   identical(x[[1]], quote(D)) ||
                   identical(x[[1]], quote(diff))) &&
                   length(x) == 3){
            return(changeDerivs(as.character(x[[2]]), as.character(x[[3]])));
        } else if ((identical(x[[1]], quote(Derivative)) ||
                    identical(x[[1]], quote(D)) ||
                    identical(x[[1]], quote(diff))) &&
                   length(x) == 4){
            changeDerivs(as.character(x[[2]]), as.character(x[[3]]), as.character(x[[4]]));
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

for (f in names(rxDefinedDerivatives)){
    rxSymPyC[[f]] <- functionOp(f);
    sympyRxFEnv[[f]] <- functionOp(f);
    rxSymPyFEnv[[f]] <- functionOp(f);
}

unknownSymPy <- function(op){
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
    callList <- setNames(lapply(calls, unknownSymPy), calls)
    callEnv <- list2env(callList);
    rxSymPyFEnv <- cloneEnv(rxSymPyFEnv, callEnv);
    names <- allNames(expr)
    ## Replace time with t.
    n1 <- names;
    n2 <- names;
    n2[n2 == "time"] <- "t";
    ## Replace f with rx_pred_
    n2[n2 == "f"] <- "rx_pred_"
    n2 <- gsub(rex::rex("."), "__DoT__", n2)
    ## Replace print functions with nothing.
    n2[regexpr(regPrint, n2) != -1] <- "";
    res <- rxSymPyReserved()
    res <- res[res != "pi"];
    w <- which(n2 %in% res);
    n2[w] <- sprintf("rx_SymPy_Res_%s", n2[w]);
    n2 <- gsub(rex::rex("rx_underscore_"), "_", n2);
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
    n2 <- gsub(rex::rex(start, "rx_SymPy_Res_"), "",
               gsub(regIni0, "\\1(0)",
                    gsub(regDfDy, "df(\\1)/dy(\\2)",
                         gsub(regDfDyTh, "df(\\1)/dy(\\2[\\3])",
                              gsub(regDDt, "d/dt(\\1)",
                                   gsub(rex::rex(start, regThEt, end), "\\1[\\2]", names))))));
    n2[n2 == "time"] <- "t";
    n2 <- gsub(rex::rex("__DoT__"), ".", n2)
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSymPyFEnv);
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
            txt <- strsplit(txt, rex::rex(or("=", "~", "<-")));
            vars <- c();
            addNames <- TRUE;
            txt <- unlist(lapply(txt, function(x){
                tmp <- x[1]
                var <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", x[1]))));
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    eq <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", x[2]))));
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
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", paste(deparse(paste(as.vector(x), collapse="\n"), collapse=""))))))
            return(txt);
        }
    } else if (class(substitute(x)) == "name"){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (any(cls == c("list", "rxDll", "RxCompilationManager", "RxODE", "solveRxDll"))){
            ret <- strsplit(rxNorm(x),"\n")[[1]];
            ret <- rxRmIni(ret);
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s,envir=envir)", paste(deparse(paste0(as.vector(ret), collapse="\n")), collapse=""))), envir=envir));
            return(txt);
        } else if (cls == "character" && length(cls) == 1){
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", paste(deparse(as.vector(x)), collapse="")))));
            return(txt);
        } else {
            expr <- evalPrints(substitute(x), envir=envir)
            txt <- eval(expr, sympyEnv(expr))
            return(paste0(txt))
        }
    } else {
        expr <- evalPrints(substitute(x), envir=envir)
        txt <- eval(expr, sympyEnv(expr));
        return(paste0(txt))
    }
}

##' @rdname rxToSymPy
##' @export
rxFromSymPy <- function(x, envir=parent.frame(1)) {
    if (class(substitute(x)) == "character"){
        if (length(x) == 1){
            txt <- strsplit(x, "\n+")[[1]];
            txt <- strsplit(txt, "[=~]", txt);
            vars <- c();
            addNames <- TRUE;
            txt <- unlist(lapply(txt, function(x){
                var <- paste0(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", x[1]))));
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    eq <- paste0(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", x[2]))));
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
            txt <- sprintf(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", deparse(paste(x, collapse="\n"))))));
            return(txt);
        }
    } else if (class(substitute(x)) == "name"){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (cls == "character" && length(cls) == 1){
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", deparse(x)))));
            return(txt);
        } else {
            expr <- evalPrints(substitute(x), envir=envir)
            txt <- eval(expr, rxEnv(expr));
            return(paste0(txt))
        }
    } else {
        expr <- evalPrints(substitute(x), envir=envir)
        txt <- eval(expr, rxEnv(expr))
        return(paste0(txt))
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
    estN <- suppressWarnings(as.numeric(est));
    if (is.na(estN)){
        if (rxErrEnv.diag.xform == "sqrt"){
            ret <- (sprintf("(%s)^2", est))
        } else if (rxErrEnv.diag.xform == "log"){
            ret <- (sprintf("exp(%s)", est))
        } else {
            ret <- (sprintf("%s", est));
        }
    } else {
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
    }
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
    estN <- suppressWarnings(as.numeric(est));
    if (is.na(estN)){
        if (rxErrEnv.diag.xform == "sqrt"){
            ret <- (sprintf("rx_pred_^2 * (%s)^2", est))
        } else if (rxErrEnv.diag.xform == "log"){
            ret <- (sprintf("rx_pred_^2 * exp(%s)", est))
        } else {
            ret <- (sprintf("rx_pred_^2 * %s", est))
        }
    } else {
        est <- estN
        ret <- ""
        theta <- sprintf("THETA[%s]", rxErrEnv.theta);
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

rxSimpleExprP <- function(x){
    if (is.name(x) || is.atomic(x)){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

##' This function splits a function based on + or - terms
##'
##' It uses the parser and does not disturb terms within other
##' functions.  For example:
##'
##' a*exp(b+c)+d*log(e-f)-g*f
##'
##' would return
##'
##' c("a * exp(b + c)", "d * log(e - f)", "- g * f")
##'
##' @param x Quoted R expression for splitting
##' @param level Internal level of parsing
##' @return character vector of the split expressions
##' @author Matthew L. Fidler
rxSplitPlusQ <- function(x, level=0){
    if (class(x) == "character" && level == 0){
        return(eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", x))))
    }
    if (is.name(x) || is.atomic(x)){
        if (level == 0){
            return(paste(deparse(x), collapse=""));
        } else {
            return(character())
        }
    } else if (is.call(x)) { # Call recurse_call recursively
        if ((identical(x[[1]], quote(`+`)) ||
             identical(x[[1]], quote(`-`))) && level == 0){
            if (length(x) == 3){
                if (identical(x[[1]], quote(`+`))){
                    one <- paste(deparse(x[[3]]), collapse="");
                } else {
                    one <- paste("-", paste(deparse(x[[3]]), collapse=""));
                }
                tmp <- rxSplitPlusQ(x[[2]], level=0);
                if (length(tmp) > 0){
                    return(c(tmp, one))
                } else {
                    tmp <- paste(deparse(x[[2]]), collapse="");
                    return(c(tmp, one))
                }
            } else {
                ## Unary + or -
                if (identical(x[[1]], quote(`+`))){
                    one <- paste(deparse(x[[2]]), collapse="");
                } else {
                    one <- paste("-", paste(deparse(x[[2]]), collapse=""));
                }
                return(one);
            }
        } else {
            tmp <- unlist(lapply(x, rxSplitPlusQ, level=1));
            if (level == 0){
                if (length(tmp) == 0){
                    tmp <- paste(deparse(x), collapse="")
                }
            }
            return(tmp)
        }
    } else if (is.pairlist(x)) {
        ## Call recurse_call recursively
        tmp <- unlist(lapply(x, rxSplitPlusQ, level=level));
        if (level == 0){
            if (length(tmp) == 0){
                tmp <- paste(deparse(x), collapse="");
            }
        }
        return(tmp)
    } else { # User supplied incorrect input
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
}

## rxSplitPlusQ(quote(2*THETA[3]^2*centr*rx__sens_centr_BY_ETA_2__BY_THETA_3___*exp(-2*ETA[2]-2*THETA[2])-
##                    4*THETA[3]^2*centr*rx__sens_centr_BY_THETA_3___*exp(-2*ETA[2]-2*THETA[2])+
##                    2*THETA[3]^2*rx__sens_centr_BY_ETA_2___*rx__sens_centr_BY_THETA_3___*exp(-2*ETA[2]-2*THETA[2])-
##                    4*THETA[3]*centr^2*exp(-2*ETA[2]-2*THETA[2])+
##                    4*THETA[3]*centr*rx__sens_centr_BY_ETA_2___*exp(-2*ETA[2]-2*THETA[2])))

## Stop DSL

rm(f);
rm(op);
rm(p);
