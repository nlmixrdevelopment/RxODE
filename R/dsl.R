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

sympyRxFEnv <- new.env(parent = emptyenv())
rxSymPyFEnv <- new.env(parent = emptyenv())

sympyRxFEnv$"(" <- unaryOp("(", ")")
rxSymPyFEnv$"(" <- unaryOp("(", ")")

rxSymPyFEnv$solveLinB <- function(...){
    return("rx1c")
}
rxSymPyFEnv$diff <- rxSymPyFEnv$Derivative
rxSymPyFEnv$D <- rxSymPyFEnv$Derivative;

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

sympyRxFEnv$"[" <- function(name, val){
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

dsl.strip.paren <- function(x){
    strip.it <- function(x){
        if (is.call(x)){
            if (length(x) == 1){
                return(x)
            } else if (as.character(x[[1]]) == "("){
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



sympyRxFEnv$dt <- function(w) {return(sprintf("dt(%s)", w))}
rxSymPyFEnv$"^" <- binaryOp("**")
rxSymPyFEnv$"**" <- binaryOp("**")
sympyRxFEnv$"**" <- dsl.to.pow
sympyRxFEnv$"^" <- dsl.to.pow

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

rxSymPyFEnv$rate <- function(e1){
    paste0("rx_rate_",e1,"_");
}
rxSymPyFEnv$r <- rxSymPyFEnv$rate

rxSymPyFEnv$dur <- function(e1){
    paste0("rx_dur_",e1,"_");
}
rxSymPyFEnv$d <- rxSymPyFEnv$dur

rxSymPyFEnv$f <- function(e1){
    paste0("rx_f_",e1,"_");
}
rxSymPyFEnv$F <- rxSymPyFEnv$f

rxSymPyFEnv$lag <- function(e1){
    paste0("rx_lag_",e1,"_");
}

rxSymPyFEnv$alag <- rxSymPyFEnv$lag

rxSymPyFEnv$cmt  <- function(e1){
    return("");
}

rxSymPyFEnv$dvid  <- function(...){
    return("");
}

.rxMtimes <- c();

## rxSymPyFEnv$mtime <- function(e2){
##     assignInMyNamespace(".rxMtimes", unique(c(paste0(e2), .rxMtimes)));
##     return(e2);
## }

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
sympyRxFEnv$Abs <- functionOp("abs")

## Following R functions are not translated
## ftrunc, fround, fprec, fsign, sign, fmin2, fmax2, imin2, imax2, logspace_sum
## choose lchoose

rxSymPyFEnv$R_pow <- binaryOp2("**");
rxSymPyFEnv$R_pow_di <- binaryOp2("**");
rxSymPyFEnv$Rx_pow <- binaryOp2("**");
rxSymPyFEnv$Rx_pow_di <- binaryOp2("**");
rxSymPyFEnv$log1p <- functionOp2("log(1 + (", "))");
rxSymPyFEnv$log1pmx <- functionBrewx("(log(1 + (<%=x%>))-(<%=x%>))");
rxSymPyFEnv$expm1 <- functionOp2("(exp(", ")-1)");
rxSymPyFEnv$log10 <- function(x){
    if (x == "E" || x == "exp(1)" || x == "e"){
        return("1/log(10)")
    } else {
        return(paste0("(log(", x, ")/log(10))"))
    }
}
rxSymPyFEnv$abs <- function(e1) {
    .e1 <- paste(e1);
    if (.e1 == "0") return("rx_eff_abs_0__");
    return(paste0("sqrt((", e1, ")*(", e1, "))"))
}


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
sympy.equiv.f <- c("acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
                   "cos", "cosh", "digamma", "erf", "erfc", "exp", "factorial",
                   "gamma", "sin", "sinh", "sqrt", "tan",
                   "tanh", "trigamma", "rxTBS", "rxTBSd")
for (f in sympy.equiv.f){
    rxSymPyFEnv[[f]] <- functionOp(f);
    sympyRxFEnv[[f]] <- functionOp(f);
}

rxSymPyFEnv[["log"]] <- functionOp("log");
sympyRxFEnv[["log"]] <- functionOp("log");

dsl.factor.pi.1 <- function(x){
    mult.split <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s),mult=TRUE)", x)));
    w <- which(mult.split == "M_PI");
    if (length(w) == 1){
        r <- gsub(rex::rex(" * 1/"), " * ", paste(mult.split[-w], collapse=" * "));
        if (r == "")
            r <- "1"
        return(r);
    }
    w <- which(mult.split == "M_2PI");
    if (length(w) >= 1){
        mult.split[w[1]] <- "2"
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    w <- which(mult.split == "M_PI_2");
    if (length(w) >=  1){
        mult.split[w[1]] <- "1/2"
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    w <- which(mult.split == "M_PI_4");
    if (length(w) >= 1){
        mult.split[w[1]] <- "1/4"
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    w <- which(mult.split == "M_1_PI")
    if (length(w) >= 1){
        mult.split[w[1]] <- "Rx_pow_di(M_PI,-2)";
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    w <- which(mult.split == "M_2_PI")
    if (length(w) >= 1){
        mult.split[w[1]] <- "2 * Rx_pow_di(M_PI,-2)";
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    w <- which(mult.split == "M_2_SQRTPI")
    if (length(w) >= 1){
        mult.split[w[1]] <- "2 * Rx_pow(M_PI,-1.5)";
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    w <- which(mult.split == "M_SQRT_PI")
    if (length(w) >= 1){
        mult.split[w[1]] <- "2 * Rx_pow(M_PI,-1.5)";
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")));
    }
    reg <- rex::rex(start, any_spaces, capture("Rx_pow", or("_di", ""), "("), any_spaces,
                    capture(except_some_of(", "), "PI", except_some_of(", ")), any_spaces, ",", any_spaces,
                    capture(except_some_of(", )")), any_spaces, ")")
    w <- which(regexpr(reg, mult.split) != -1)
    if (length(w) >= 1){
        num <- suppressWarnings(as.numeric(gsub(reg, "\\3", mult.split[w])));
        if (is.na(num)){
            pow <- paste0(gsub(reg, "\\3", mult.split[w]), " - 1");
            mult.split[w] <- gsub(reg, sprintf("\\1\\2,%s)", pow), mult.split[w])
        } else {
            pow <- num - 1;
            if (pow == 1){
                mult.split[w] <- gsub(reg, "\\2", mult.split[w])
            } else {
                mult.split[w] <- gsub(reg, sprintf("\\1\\2,%s)", pow), mult.split[w])
            }
        }
        return(gsub(rex::rex(" * 1/"), " * ", paste(mult.split, collapse=" * ")))
    }
    return(NA)
}


dsl.factor.pi <- function(x, op="cos"){
    factored <- sapply(eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", x))), dsl.factor.pi.1);
    if (any(is.na(factored))){
        return(sprintf("%s(%s)", op, x))
    } else {
        return(sprintf("%spi(%s)", op, gsub(rex::rex("+ -"), "-", paste(factored, collapse=" + "))))
    }

}

sympyRxFEnv$cos <- function(x){
    ## Convert to cospi
    return(dsl.factor.pi(x))
}
sympyRxFEnv$sin <- function(x){
    ## Convert to cospi
    dsl.factor.pi(x, "sin")
}

sympyRxFEnv$tan <- function(x){
    ## Convert to cospi
    dsl.factor.pi(x, "tan")
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


dsl.handle.log <- function(x, abs=FALSE){
    ## Put more accurate functions in...
    ## lgamma1p
    x <- dsl.strip.paren(x);
    num <- suppressWarnings(as.numeric(x))
    if (!is.na(num)){
        if (num == 2){
            return("M_LN2");
        } else if (num == 10){
            return("M_LN10");
        }
    }
    if (any(x == c("M_SQRT_PI", "sqrt(M_PI)"))){
        return("M_LN_SQRT_PI");
    } else if (any(x == c("sqrt(2 * M_PI)", "sqrt(M_PI * 2)", "sqrt(M_2PI)"))){
        return("M_LN_SQRT_2PI");
    } else if (any(x == c("sqrt(M_PI/2)", "sqrt(M_PI_2)"))){
        return("M_LN_SQRT_PId2");
    }
    pls <- try(eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", x))), silent=TRUE);
    if (inherits(pls, "try-error")){
        return(sprintf("%slog(%s)", ifelse(abs, "abs_", ""), x))
    } else if (length(pls) == 1) {
        reg <- rex::rex(start, any_spaces, "gamma(", capture(anything), ")", any_spaces, end)
        new <- sub(reg, "\\1", pls);
        if (nchar(new) != nchar(pls)){
            new.expand <- rxSymPyExpand(new);
            pls <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", new.expand)));
            pls.n <- suppressWarnings({as.numeric(pls)})
            w <- which(!is.na(pls.n));
            if (length(w) == 1){
                n <- pls.n[w] - 1;
                if (n == 0){
                    pls <- pls[-w];
                } else {
                    pls[w] <- paste(n);
                }
                ret <- gsub(rex::rex("+-"), "-", paste(sapply(pls, rxFromSymPy), collapse="+"))
                return(sprintf("lgamma1p(%s))", ret))
                return(sprintf("lgamma1p(%s)", ret))
            }
        }
    }
    ## log1pexp
    w <- which(pls == "1");
    if (length(pls) == 2 && length(w) == 1){
        pls <- pls[-w];
        reg <- rex::rex(start, any_spaces, "exp(", capture(anything), ")", any_spaces, end)
        new <- sub(reg, "log1pexp(\\1)", pls);
        if (nchar(new) != nchar(pls)){
            return(new)
        } else {
            return(sprintf("%slog1p(%s)", ifelse(abs, "abs_", ""), pls))
        }

    }
    ## log1p
    ## log(1+x) or log(x + 1) or log(y + 1 + z)
    ## First expand
    x.expand <- rxSymPyExpand(x);
    ## Then split on +
    pls <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", x.expand)));
    ## Let see if any of these are numeric...
    pls.n <- suppressWarnings({as.numeric(pls)})
    w <- which(!is.na(pls.n));
    if (length(w) == 1){
        n <- pls.n[w] - 1;
        if (n == 0){
            pls <- pls[-w];
        } else {
            pls[w] <- paste(n);
        }
        for (i in seq_along(pls)){
            tmp <- pls[i]
            pls[i] <- rxFromSymPy(tmp);
        }
        ret <- gsub(rex::rex("+-"), "-", paste(pls, collapse="+"))
        return(sprintf("log1p(%s)", ret))
    }
    p1 <- rex::rex(start, any_spaces, "1", any_spaces, "+", any_spaces)
    p2 <- rex::rex(any_spaces, "+", any_spaces, "1", any_spaces, end)
    if (regexpr(p2, x) != -1){
        x2 <- gsub(p2, "", x);
        return(sprintf("log1p(%s)", x2))
    } else if (regexpr(p1, x) != -1){
        x2 <- gsub(p1, "", x);
        return(sprintf("log1p(%s)", x2))
    } else {
        return(sprintf("%slog(%s)", ifelse(abs, "abs_", ""), x));
    }
}
rxSymPyFEnv$log2 <- function(e1){
    if (e1 == "E" || e1 == "exp(1)"){
        return("1/log(2)");
    } else {
        return(paste0("log(", e1, ")/log(2)"));
    }
}

sympyRxFEnv$log10 <- function(e1){
    if (e1 == "M_E" || e1 == "exp(1)" || e1 == "e"){
        return("M_LOG10E");
    } else {
        return(paste0("log(", e1, ")/log(10)"));
    }
}
sympyRxFEnv[["/"]] <- function(e1, e2, sep="/"){
    if (e2 == "M_LN2"){
        if (e1 == 1){
            ## log(e)/log(2)
            return("M_LOG2E");
        }
        reg <- rex::rex(start, any_spaces, "log(");
        if (regexpr(reg, e1) != -1){
            return(sub(reg, "log2(", e1));
        }
        return(paste0(e1, " * M_LOG2E"))
    }
    if (e1 == "M_PI"){
        if (e2 == 2){
            return("M_PI_2")
        } else if (e2 == 4){
            return("M_PI_4")
        }
    }
    if (e1 == 1){
        if (e2 == "M_PI"){
            return("M_1_PI")
        } else if (e2 == "sqrt(2)"){
            return("M_SQRT1_2");
        } else if (e2 == "sqrt(M_2PI)"){
            return("M_1_SQRT_2PI");
        }
    }
    if (e1 == 2){
        if (e2 == "M_PI"){
            return("M_2_PI")
        } else if (e2 == "M_SQRT_PI"){
            return("M_2_SQRTPI")
        }
    }
    if (e2 == "M_LN10"){
        if (e1 == 1){
            ## log(e)/log(10)
            return("M_LOG10E");
        }
        reg <- rex::rex(start, any_spaces, "log(");
        if (regexpr(reg, e1) != -1){
            return(sub(reg, "log10(", e1));
        }
        return(paste0(e1, " * M_LOG10E"))
    }
    return(paste0(e1, sep, e2))
}
sympyRxFEnv[["*"]] <- function(e1, e2, sep=" * "){
    if ((e1 == 2 && e2 == "M_PI") ||
        (e2 == 2 && e1 == "M_PI")){
        return("M_2PI")
    } else {
        ##  You shouldn't need to parse additive expressions first, since R handles multiplication before +/-
        if (e2 == 2){
            mult.split <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s),mult=TRUE)", e1)))
            w.pi <- which(mult.split == "M_PI");
            if (length(w.pi) == 1){
                return(gsub(rex::rex(sep, "1/"), sep, paste0(paste(mult.split[-w.pi], collapse=sep), sep, "M_2PI")));
            }
        }
        if (e2 == "M_PI"){
            mult.split <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s),mult=TRUE)", e1)))
            w.2 <- which(mult.split == "2");
            if (length(w.2) == 1){
                return(gsub(rex::rex(sep, "1/"), sep, paste0(paste(mult.split[-w.2], collapse=sep), sep, "M_2PI")));
            }
        }
        return(paste0(e1, sep, e2))
    }
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
        return(dsl.handle.log(arg, abs=TRUE));
    } else {
        return(dsl.handle.log(arg));
    }
}

sympyRxFEnv$sqrt <- function(arg){
    arg <- dsl.strip.paren(arg);
    if (any(arg == c("2/M_PI", "M_2_PI"))){
        return("M_SQRT_2dPI");
    } else if (arg == "M_PI"){
        return("M_SQRT_PI");
    } else {
        num <- suppressWarnings(as.numeric(arg));
        if (!is.na(num)){
            if (num == 32){
                return("M_SQRT_32")
            } else if (num == 3){
                return("M_SQRT_3");
            } else if (num == 2){
                return("M_SQRT_2");
            }
        }
        return(sprintf("sqrt(%s)", arg))
    }
}

rxSymPyFEnv$structure <- function(one, ..., .Names){
    eval(parse(text=sprintf("rxToSymPy(%s)", deparse(sprintf("%s", one)))));
}

sympyRxFEnv$structure <- function(one, ..., .Names){
    eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", deparse(sprintf("%s", one)))));
}

sympyRxFEnv$Subs <- function(...){
    .lst <- list(...);
    if (length(.lst) == 3){
        .what <- dsl.strip.paren(.lst[[2]]);
        .with <- dsl.strip.paren(.lst[[3]]);
        .expr <- dsl.strip.paren(.lst[[1]]);
        if (regexpr("[,]", .what) != -1){
            stop("Can't handle Subs with multiple replacements (yet)")
        }
        .subs <- function(x){
            if (all(as.character(x) == .what)){
                return(eval(parse(text=sprintf("quote(%s)", .with))));
            } else if (is.call(x)) {
                as.call(lapply(x, .subs))
            } else if (is.pairlist(x)) {
                as.pairlist(lapply(x, .subs))
            } else {
                return(x);
            }
        }
        .expr <- eval(parse(text=sprintf(".subs(quote(%s))", .expr)))
        .expr <- paste(deparse(.expr), collapse=" ");
        return(.expr);
    }
}

sympyRxFEnv$subs <- sympyRxFEnv$Subs;

rxSymPyFEnv$psigamma <- function(z, n){
    paste0("polygamma(", n, ", ", z, ")");
}
sympyRxFEnv$polygamma <- function(n, z){
    paste0("psigamma(", z, ", ", n, ")");
}

## Add sympy->C mini DSL for omega parsing

rxSymPyC <- new.env(parent = emptyenv())
rxSymPyC$"**" <- function(a, b){
    a <- dsl.strip.paren(a);
    b <- dsl.strip.paren(b);
    num <- suppressWarnings({as.numeric(b)});
    if (is.na(num)){
        return(sprinf("R_pow(%s, %s)", a, b))
    } else if (num == round(num)) {
        return(sprintf("R_pow_di(%s, %s)", a, b));
    } else if (num == 0.5){
        return(sprintf("sqrt(%s)", a));
    } else {
        return(sprintf("Rx_pow(%s, %s)", a, b));
    }
}
rxSymPyC$"^" <- rxSymPyC$"**"

rxSymPyC$S <- function(x){
    sprintf("%s", x);
}

for (f in sympy.equiv.f){
    rxSymPyC[[f]] <- functionOp(f);
}
rxSymPyC$"(" <- unaryOp("(", ")")
for (op in c("+", "-", "*")){
    rxSymPyC[[op]] <- binaryOp(paste0(" ", op, " "));
}

rxSymPyC[["/"]] <- function(e1, e2){
    sprintf("%s /( (%s == 0) ? %s : %s)", e1, e2, .Machine$double.eps, e2)
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
    n2 <- gsub(rex::rex("rx_SymPy_Res_"), "", n2)
    n2 <- gsub("None", "NA_REAL", n2);
    w <- n2[n2 == "t"];
    symbol.list <- setNames(as.list(n2), n1);
    symbol.env <- list2env(symbol.list, parent=rxSymPyC);
    return(symbol.env)
}

sympyC <- function(x){
    if (is(substitute(x),"character")){
        return(eval(parse(text=sprintf("RxODE:::sympyC(quote(%s))", x))))
    } else if (is(substitute(x),"name")){
        return(eval(parse(text=sprintf("RxODE:::sympyC(%s)", deparse(x)))));
    } else {
        return(eval(x, sympyCEnv(x)))
    }
}


## nocov end

sympyTransit4 <- function(t, n, mtt, bio, podo="podo", tlast="tlast"){
    ktr <- paste0("((", n, " + 1)/(", mtt, "))");
    lktr <- paste0("(log((", n, ") + 1) - log(", mtt, "))");
    tc <- paste0("((", t, ")-(", tlast, "))");
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
    if (is(substitute(x),"character")){
        force(x);
        if (length(x) == 1){
            names(x) <- NULL;
            txt  <- gsub(rex::rex(boundary,or("cmt","dvid"),"(",except_any_of(")"),")"), "", x,perl=TRUE)
            txt <- strsplit(gsub(";", "\n", txt), "\n+")[[1]];
            txt <- strsplit(txt, rex::rex(or("=", "~", "<-")));
            tmp <- unlist(lapply(txt, function(x){length(x)}));
            if (length(tmp) > 1){
                if (all(tmp == 1)){
                    txt <- paste(txt, collapse=" ");
                } else if (any(tmp == 2) && any(tmp == 1)){
                    txt2 <- list()
                    for (i in seq_along(txt)){
                        if (length(txt[[i]]) == 2){
                            txt2[[length(txt2) + 1]] <- txt[[i]]
                        } else {
                            tmp <- txt2[[length(txt2)]];
                            tmp[2] <- gsub(" +", "", paste(tmp[2], txt[[i]]));
                            txt2[[length(txt2)]] <- tmp
                        }
                    }
                    txt <- txt2
                }
            }
            vars <- c();
            addNames <- TRUE;
            txt <- unlist(lapply(txt, function(x){
                tmp <- sub(rex::rex(any_spaces, end), "", sub(rex::rex(start, any_spaces), "", x[1]))
                if (.exists2(tmp, envir)){
                    res <- rxSymPyReserved()
                    if (any(tmp == res)){
                        var <- paste0("rx_SymPy_Res_", tmp)
                    } else {
                        var <- tmp;
                    }
                } else {
                    var <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", tmp)), envir=envir));
                }
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    tmp <- sub(rex::rex(any_spaces, end), "", sub(rex::rex(start, any_spaces), "", x[2]))
                    if (.exists2(tmp, envir)){
                        res <- rxSymPyReserved()
                        if (any(tmp == res)){
                            eq <- paste0("rx_SymPy_Res_", tmp)
                        } else {
                            eq <- tmp;
                        }
                    } else {
                        eq <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", tmp)), envir=envir));
                    }
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
            force(x);
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", paste(deparse(paste(as.vector(x), collapse="\n")))))), collapse="")
            return(txt);
        }
    } else if (is(substitute(x),"name")){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (any(cls == c("list", "rxDll", "RxCompilationManager", "RxODE", "solveRxDll", "rxModelVars"))){
            force(x);
            ret <- strsplit(rxNorm(x),"\n")[[1]];
            ret <- .rxRmIni(ret);
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", paste(deparse(paste0(as.vector(ret), collapse="\n")), collapse=""))), envir=envir));
            return(txt);
        } else if (cls == "character" && length(cls) == 1){
            force(x)
            txt <- paste0(eval(parse(text=sprintf("RxODE::rxToSymPy(%s)", paste(deparse(as.vector(x)), collapse="")))));
            return(txt);
        } else {
            expr <- evalPrints(substitute(x), envir=envir)
            txt  <- eval(expr, sympyEnv(expr))
            txt  <- gsub(rex::rex(" * 1/"), "/", txt);
            return(paste0(txt))
        }
    } else {
        expr <- evalPrints(substitute(x), envir=envir)
        txt <- eval(expr, sympyEnv(expr));
        txt <- gsub(rex::rex(" * 1/"), "/", txt);
        return(paste0(txt))
    }
}

##' @rdname rxToSymPy
##' @export
rxFromSymPy <- function(x, envir=parent.frame(1)) {
    if (is(substitute(x),"character")){
        force(x);
        if (length(x) == 1){
            names(x) <- NULL;
            txt <- strsplit(x, "\n+")[[1]];
            txt <- strsplit(txt, "[=~]", txt);
            vars <- c();
            addNames <- TRUE;
            tmp <- unlist(lapply(txt, function(x){length(x)}));
            if (length(tmp) > 1){
                if (all(tmp == 1)){
                    txt <- paste(txt, collapse=" ");
                } else if (any(tmp == 2) && any(tmp == 1)){
                    txt2 <- list()
                    for (i in seq_along(txt)){
                        if (length(txt[[i]]) == 2){
                            txt2[[length(txt2) + 1]] <- txt[[i]]
                        } else {
                            tmp <- txt2[[length(txt2)]];
                            tmp[2] <- gsub(" +", "", paste(tmp[2], txt[[i]]));
                            txt2[[length(txt2)]] <- tmp
                        }
                    }
                    txt <- txt2
                }
            }
            txt <- unlist(lapply(txt, function(x){
                tmp <- sub(rex::rex(any_spaces, end), "", sub(rex::rex(start, any_spaces), "", x[1]))
                if (.exists2(tmp, envir)){
                    var <- sub(rex::rex(start, "rx_SymPy_Res_"), "", tmp)
                } else {
                    var <- paste0(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", tmp)), envir=envir));
                }
                if (length(x) == 2){
                    vars <<- c(vars, var);
                    tmp <- sub(rex::rex(any_spaces, end), "", sub(rex::rex(start, any_spaces), "", x[2]))
                    if (.exists2(tmp, envir)){
                        e1 <- sub(rex::rex(start, "rx_SymPy_Res_"), "", tmp)
                    } else {
                        eq <- paste0(eval(parse(text=sprintf("RxODE::rxFromSymPy(%s)", x[2])), envir=envir));
                    }
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
    } else if (is(substitute(x),"name")){
        cls <- tryCatch({class(x)}, error=function(e){return("error")});
        if (cls == "character" && length(cls) == 1){
            force(x);
            names(x) <- NULL
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
##' @param ret Internal return type.  Should not be changed by the user...
##' @param init Initialization vector
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
##' @param mult boolean to split based on * and / expressions instead.
##'     By default this is turned off.
##' @return character vector of the split expressions
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxSplitPlusQ <- function(x, level=0, mult=FALSE){
    if (is(x,"character") && level == 0){
        return(eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", x))))
    }
    if (is.name(x) || is.atomic(x)){
        if (level == 0){
            return(paste(deparse(x), collapse=""));
        } else {
            return(character())
        }
    } else if (is.call(x)) { # Call recurse_call recursively
        if ((mult && ((identical(x[[1]], quote(`*`)) ||
                       identical(x[[1]], quote(`/`))) && level == 0)) ||
            (!mult && ((identical(x[[1]], quote(`+`)) ||
                        identical(x[[1]], quote(`-`))) && level == 0))){
            if (length(x) == 3){
                if (identical(x[[1]], quote(`+`))){
                    one <- paste(deparse(x[[3]]), collapse="");
                } else if (!mult){
                    one <- paste("-", paste(deparse(x[[3]]), collapse=""));
                } else if (identical(x[[1]], quote(`*`))){
                    one <- paste(deparse(x[[3]]), collapse="");
                } else if (mult){
                    one <- paste("1/", paste(deparse(x[[3]]), collapse=""));
                }
                tmp <- rxSplitPlusQ(x[[2]], level=0, mult=mult);
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
            tmp <- unlist(lapply(x, rxSplitPlusQ, level=1, mult=mult));
            if (level == 0){
                if (length(tmp) == 0){
                    tmp <- paste(deparse(x), collapse="")
                }
            }
            return(tmp)
        }
    } else if (is.pairlist(x)) {
        ## Call recurse_call recursively
        tmp <- unlist(lapply(x, rxSplitPlusQ, level=level, mult=mult));
        if (level == 0){
            if (length(tmp) == 0){
                tmp <- paste(deparse(x), collapse="");
            }
        }
        return(tmp)
    } else { # User supplied incorrect input
        stop("Don't know how to handle type '", typeof(x), "'.",
             call. = FALSE)
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
    rxSymPySetup(model);
    cnd <- rxNorm(model, TRUE);
    lines <- strsplit(rxNorm(model), "\n")[[1]];
    for (i in seq_along(lines)){
        if (regexpr("[=~]", lines[i])){
            type <- sub(".*([=~]).*", "\\1", lines[i]);
            l0 <- strsplit(lines[i], "[=~]")[[1]];
            l2 <- substr(l0[2], 1, nchar(l0[2]) - 1);
            if (expand){
                l2 <- rxSymPy(sprintf("expand(%s)", rxToSymPy(l2)))
                l2 <- rxFromSymPy(l2)
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
rm(p);
