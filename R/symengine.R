## first.arg second.arg and type of function
.rxSEsingle <- list("gammafn"=c("gamma(", ")", "gamma"),
                    "lgammafn"=c("loggamma(", ")", "lgamma"),
                    "lgamma"=c("loggamma(", ")", "lgamma"),
                    "digamma"=c("polygamma(0,", ")", "psigamma"),
                    "trigamma"=c("polygamma(1,", ")", "psigamma"),
                    "tetragamma"=c("polygamma(2,", ")", "psigamma"),
                    "pentagamma"=c("polygamma(3,", ")", "psigamma"),
                    "cospi"=c("cos(pi*(", "))", "cos"),
                    "sinpi"=c("sin(pi*(", "))", "sin"),
                    "tanpi"=c("tan(pi*(", "))", "tan"),
                    "log1p"=c("log(1+", ")", "log"),
                    "expm1"=c("(exp(", ")-1)", "exp"),
                    "factorial"=c("gamma(", "+1)", "gamma"),
                    "lgamma1p"=c("loggamma(", "+1)", "lgamma"),
                    "expm1"=c("(exp(", ")-1)", "exp"),
                    "log10"=c("log(", ")/log(10)", "log"),
                    "log2"=c("log(", ")/log(2)", "log"),
                    "log1pexp"=c("log(1+exp(", "))", "log1pexp"),
                    "!"=c("rxNot(", ")", "")
                    )

.SEsingle <- list("loggamma"=c("lgamma(", ")"),
                  "rxNot"=c("(!(", "))"))

.rxSEdouble <- list("pow"=c("(", ")^(", ")"),
                    "R_pow"=c("(", ")^(", ")"),
                    "R_pow_di"=c("(", ")^(", ")"),
                    "Rx_pow_di"=c("(", ")^(", ")"),
                    "Rx_pow"=c("(", ")^(", ")"),
                    "lbeta"=c("log(beta(", ",", "))"),
                    "=="=c("rxEq(", ",", ")"),
                    "!="=c("rxNeq(", ",", ")"),
                    ">="=c("rxGeq(", ",", ")"),
                    "<="=c("rxLeq(", ",", ")"),
                    "<"=c("rxLt(", ",", ")"),
                    ">"=c("rxGt(", ",", ")"),
                    "&&"=c("rxAnd(", ",", ")"),
                    "||"=c("rxOr(", ",", ")"),
                    "&"=c("rxAnd(", ",", ")"),
                    "|"=c("rxOr(", ",", ")")
                    )

.SEdouble <- list("lbeta"=c("lbeta(", ",", ")"),
                  "rxEq"=c("(", "==", ")"),
                  "rxNeq"=c("(", "!=", ")"),
                  "rxGeq"=c("(", ">=", ")"),
                  "rxLeq"=c("(", "<=", ")"),
                  "rxGt"=c("(", ">", ")"),
                  "rxLt"=c("(", "<", ")"),
                  "rxAnd"=c("(", "&&", ")"),
                  "rxOr"=c("(", "||", ")"))

## atan2
.rxSEeq <- c("acos"=1, "acosh"=1, "asin"=1, "atan"=1,
             "atanh"=1, "beta"=2,
             "cos"=1, "cosh"=1, "erf"=1, "erfc"=1,
             "exp"=1, "gamma"=1, "sin"=1, "sinh"=1,
             "sqrt"=1, "tan"=1, "tanh"=1, "log"=1, "abs"=1, "asinh"=1,
             "rxTBS"=3, "rxTBSd"=3, "rxTBSd2"=3,
             "solveLinB"=19)

.rxSEeqUsr <- c()

.rxCcode <- c()
.symengineFs <- new.env(parent=emptyenv())

.extraCnow <- "";
.extraC <- function(extraC=NULL){
    if (!is.null(extraC)){
        if (file.exists(extraC)){
            .ret <- sprintf("#include \"%s\"\n", extraC);
        } else {
            .ret <- paste(extraC, collapse="\n");
        }
    } else {
        .ret <- ""
    }
    if (length(.rxCcode) > 0L){
        .ret <- sprintf("%s\n%s\n", .ret, paste(.rxCcode, collapse="\n"));
    }
    assignInMyNamespace(".extraCnow", .ret);
    return(invisible());
}

##' Add user function to RxODE
##'
##' This adds a user function to RxODE that can be called.  If needed,
##' these functions can be differentiated by numerical differences or
##' by adding the derivatives to RxODE's internal derivative table
##' with \code{\link{rxD}}
##'
##' @param name This gives the name of the user function
##' @param args This gives the arguments of the user function
##' @param cCode This is the C-code for the new function
##' @return nothing if successful
##' @author Matthew L. Fidler
##' @examples
##'
##'
##' \dontrun{
##' ## Right now RxODE is not aware of the function f
##' ## Therefore it cannot translate it to symengine or
##' ## Compile a model with it.
##'
##' RxODE("a=fun(a,b,c)")
##'
##' ## Note for this approach to work, it cannot interfere with C
##' ## function names or reserved RxODE specical terms.  Therefore
##' ## f(x) would not work since f is an alias for bioaviability.
##' }
##'
##' fun <- "
##' double fun(double a, double b, double c) {
##'   return a*a+b*a+c;
##' }
##' " ## C-code for function
##'
##' rxFun("fun",c("a","b","c"), fun) ## Added function
##'
##' ## Now RxODE knows how to translate this function to symengine
##'
##' rxToSE("f(a,b,c)")
##'
##' ## And will take a central difference when calculating derivatives
##'
##' rxFromSE("Derivative(f(a,b,c),a)")
##'
##' ## Of course, you could specify the derivative table manually
##' rxD("fun", list(function(a,b,c){
##'   paste0("2*",a,"+",b);
##' },
##'     function(a,b,c){
##'      return(a)
##'     },
##'     function(a,b,c){
##'       return("0.0")
##' } ))
##'
##'
##' @export
rxFun <- function(name, args, cCode){
    if (!is.character(name) || length(name) != 1L)
        stop("name argument must be a length-one character vector")
    if (missing(cCode)) stop("A new function requires a C function so it can be used in RxODE")
    if (any(name == names(.rxSEeqUsr))){
        stop(sprintf("Already defined user function '%s', remove it fist (rxRmFun).", name));
    }
    assignInMyNamespace(".rxSEeqUsr", c(.rxSEeqUsr, setNames(length(args), name)))
    assignInMyNamespace(".rxCcode", c(.rxCcode, setNames(cCode, name)))
    .symengineFs[[name]] <- symengine::FunctionSymbol(name, args)
    return(invisible())
}

##' @rdname
##' @export
rxRmFun <- function(name){
    if (!is.character(name) || length(name) != 1L)
        stop("name argument must be a length-one character vector")
    if (!any(name == names(.rxSEeqUsr))){
        stop(sprintf("No user function '%s' to remove", name))
    }
    .w <- which(name == names(.rxSEeqUsr))
    if (length(.w) == 1L) assignInMyNamespace(".rxSEeqUsr", .rxSEeqUsr[-.w])
    .w <- which(name == names(.rxCcode))
    if (length(.w) == 1L) assignInMyNamespace(".rxCcode", .rxCcode[-.w])
    if (exists(name, envir=.rxD)) rm(list=name, envir=.rxD);
    return(invisible())
}

.SE1p <-c("loggamma"="lgamma1p",
          "log"="log1p")

.SE1m <-c("cos"="cospi",
          "sin"="sinpi",
          "tan"="tanpi");

## "rxTBS", "rxTBSd"
.rxSEcnt <- c("M_E" = "E",
              "M_PI" = "pi",
              "M_PI_2" = "pi/2",
              "M_PI_4" = "pi/4",
              "M_1_PI" = "1/pi",
              "M_2_PI" = "2/pi",
              "M_2PI" = "2*pi",
              "M_SQRT_PI" = "sqrt(pi)",
              "M_2_SQRTPI" = "2/sqrt(pi)",
              "M_1_SQRT_2PI" = "1/sqrt(2*pi)",
              "M_SQRT_2" = "sqrt(2)",
              "M_SQRT_3" = "sqrt(3)",
              "M_SQRT_32" = "sqrt(32)",
              "M_SQRT_2dPI" = "sqrt(2/pi)",
              "M_LN_SQRT_PI" = "log(sqrt(pi))",
              "M_LN_SQRT_2PI" = "log(sqrt(2*pi))",
              "M_LN_SQRT_PId2" = "log(sqrt(pi/2))",
              "M_SQRT2" = "sqrt(2)",
              "M_SQRT3" = "sqrt(3)",
              "M_SQRT32" = "sqrt(32)",
              "M_LOG10_2" = "log(2)/log(10)",
              "M_LOG2E" = "1/log(2)",
              "M_LOG10E" = "1/log(10)",
              "M_LN2" = "log(2)",
              "M_LN10" = "log(10)");

## "rxTBS", "rxTBSd"

.rxSEreserved <- list("e"="M_E",
                      "E"="M_E",
                      "EulerGamma"=0.57721566490153286060651209008240243104215933593992,
                      "Catalan"=0.915965594177219015054603514932384110774,
                      "GoldenRatio"="(1+sqrt(5)/2)",
                      "I"=1i);

## diff(rxTBS(a,lambda,yj),a)
.rxD <- new.env(parent=emptyenv());
## This environment is a derivative table;
## For example:
## Derivative(f(a,b,c), a) = fa()
## Derivative(f(a,b,c), b) = fb()
## Derivative(f(a,b,c), c) = fc()
## Then
##
## .rxD$f <- list(fa(a,b,c), fb(a,b,c), fc(a,b,c))
##
##  fa translates the arguments to the derivative with respect to a
##  fb translates the arguments to the derivative with respect to b
##
## If any of the list is NULL then RxODE won't know how to take a
## derivative with respect to the argument.
##
## If the list is shorter than the length of the arguments then the
## argument then the derivative of arguments that are not specified
## cannot be taken.
.rxD$rxTBS <- list(function(a, lambda, yj){
    paste0("rxTBSd(", a, ",", lambda, ",", yj, ")")
})

.rxD$rxTBSd <- list(function(a, lambda, yj){
    paste0("rxTBSd2(", a, ",", lambda, ",", yj, ")")
})

.rxD$..k <- 10
.rxD$..tol <- 1e-4
## Approx a==b by
##(1-tanh(k*(a-b))^2)
.rxD$rxEq <- list(
    function(a, b){
    .ab <- paste0("(", a, "-", b, ")");
    return(paste0("(", -2 * .rxD$..k, "*tanh(", .rxD$..k, "*", .ab, ")+",
                  2*.rxD$..k, "*tanh(", .rxD$..k, "*", .ab, ")^3)"))
}, function(a, b){
    .ab <- paste0("(", a, "-", b, ")");
    return(paste0("(", 2 * .rxD$..k, "*tanh(", .rxD$..k, "*", .ab, ")-",
                  2*.rxD$..k, "*tanh(", .rxD$..k, "*", .ab, ")^3)"))
})

.rxD$rxGeq <- list(
    function(a, b){
    .delta <- atanh(2*.rxD$..tol-1);
    ## approx is (1/2+1/2*tanh(k*(a-b)-delta))
    .ab <- paste0("(", a, "-", b, ")");
    ## (1/2)*k + (-1/2)*k*tanh(-delta + k*(a - b))^2
    return(paste0("(", .rxD$..k / 2, "-", .rxD$..k / 2, "*tanh(", -.delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
}, function(a, b){
    .delta <- atanh(2*.rxD$..tol-1);
    ## approx is (1/2+1/2*tanh(k*(a-b)-delta))
    .ab <- paste0("(", a, "-", b, ")");
    ## (1/2)*k + (-1/2)*k*tanh(-delta + k*(a - b))^2
    return(paste0("(", -.rxD$..k / 2, "+", .rxD$..k / 2, "*tanh(", -.delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
})

.rxD$rxLeq <- list(
    function(a, b){
    .delta <- atanh(2*.rxD$..tol-1);
    ## approx is (1/2-1/2*tanh(k*(a-b)+delta))
    .ab <- paste0("(", a, "-", b, ")");
    return(paste0("(", -.rxD$..k / 2, "+", .rxD$..k / 2, "*tanh(", .delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
}, function(a, b){
    .delta <- atanh(2*.rxD$..tol-1);
    ## approx is (1/2-1/2*tanh(k*(a-b)+delta))
    .ab <- paste0("(", a, "-", b, ")");
    return(paste0("(", .rxD$..k / 2, "-", .rxD$..k / 2, "*tanh(", .delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
})


.rxD$rxLt <- list(
    function(a, b){
    ## Approx is 1/2-1/2*tanh(k*(a-b)-delta)
    .delta <- atanh(2*.rxD$..tol-1);
    .ab <- paste0("(", a, "-", b, ")");
    ## (-1/2)*k + (1/2)*k*tanh(-delta + k*(a - b))^2
    return(paste0("(", -.rxD$..k / 2, "+", .rxD$..k / 2, "*tanh(", -.delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
},
function(a, b){
    ## Approx is 1/2-1/2*tanh(k*(a-b)-delta)
    .delta <- atanh(2*.rxD$..tol-1);
    .ab <- paste0("(", a, "-", b, ")");
    ## (-1/2)*k + (1/2)*k*tanh(-delta + k*(a - b))^2
    return(paste0("(", .rxD$..k / 2, "-", .rxD$..k / 2, "*tanh(", -.delta,
                  "+", .rxD$..k, "*", .ab, ")^2)"))
})


.rxD$rxGt <- list(
    function(a, b){
    ## delta <- atanh(2*tol-1);
    ## 1/2+1/2*tanh(k*(a-b)+delta)
    .delta <- atanh(2*.rxD$..tol-1);
    .ab <- paste0("(", a, "-", b, ")");
    ## (1/2)*k + (-1/2)*k*tanh(delta + k*(a - b))^2
    return(paste0("(", .rxD$..k / 2, "-", .rxD$..k / 2, "*tanh(", .delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
},
function(a, b){
    ## delta <- atanh(2*tol-1);
    ## 1/2+1/2*tanh(k*(a-b)+delta)
    .delta <- atanh(2*.rxD$..tol-1);
    .ab <- paste0("(", a, "-", b, ")");
    ## (-1/2)*k + (1/2)*k*tanh(delta + k*(a - b))^2
    return(paste0("(", -.rxD$..k / 2, "+", .rxD$..k / 2, "*tanh(", .delta,
           "+", .rxD$..k, "*", .ab, ")^2)"))
})

.rxD$rxAnd <- list(
    function(a, b){
    ## a*b
    return(b)
}, function(a, b){
    ## a*b
    return(a)
})

.rxD$rxOr <- list(
    function(a, b){
    ## Using DeMorgan's Theorem
    ## a+b = 1-(1-a)*(1-b)
    return(paste0("(1-(", b, "))"))
}, function(a, b){
    return(paste0("(1-(", a, "))"))
})


.rxD$rxNot <- list(
    function(a){
    ## 1 - a
    return("(-1)")
})

## Approx a>=b by
## 1/2-1/2*tanh(k*x+delta)=1-tol
## 1/2-1+tol=1/2*tanh(k*x+delta)
## atanh(2*tol-1)= delta
## 1/2-1/2*tanh(k*(a-b)+delta)


##' Add to RxODE's derivative tables
##'
##' @param name Function Name
##' @param derivatives A list of functions. Each function takes the
##'     same number of arguments as the original function.  The first
##'     function will construct the derivative with respect to the
##'     first argument; The second function will construct the
##'     derivitive with respect to the second argument, and so on.
##' @return If successful, nothing
##' @author Matthew Fidler
##' @export
##' @examples
##' ## Add an arbitrary list of derivative functions
##' ## In this case the fun(x,y) is assumed to be 0.5*x^2+0.5*y^2
##'
##' rxD("fun", list(function(x,y){return(x)},
##'                 function(x,y){return(y)}
##'                 ))
##'
rxD <- function(name, derivatives){
    if (!inherits(derivatives, "list") || length(derivatives) == 0L){
        stop("Derivatives must be a list of functions with at least 1 element");
    }
    if (!all(sapply(derivatives, function(x) (inherits(x, "function") || is.null(x))))){
        stop("Derivatives must be a list of functions with at least 1 element")
    }
    if (exists(name, envir=.rxD)){
        warning(sprintf("Replacing defined derivatives for `%s`", name))
    }
    assign(name, derivatives, envir=.rxD)
    return(invisible())
}


##' RxODE to symengine.R
##'
##' @param x expression
##' @param envir \code{NULL}
##' @return
##' @author Matthew L. Fidler
##' @export
rxToSE <- function(x, envir=NULL){
    if (is(substitute(x),"character")){
        force(x);
    } else if (is(substitute(x), "{")){
        x <- deparse(substitute(x));
        if (x[1] == "{"){
            x <- x[-1];
            x <- x[-length(x)];
        }
        x <- paste(x, collapse="\n");
    } else {
        .xc <- as.character(substitute(x));
        x <- substitute(x);
        if (length(.xc == 1)){
            .found <- FALSE
            .frames <- seq(1, sys.nframe());
            .frames <- .frames[.frames != 0];
            for (.f in .frames){
                .env <- parent.frame(.f);
                if (exists(.xc, envir=.env)){
                    .val2 <- try(get(.xc, envir=.env), silent=TRUE);
                    if (inherits(.val2, "character")){
                        .val2 <- eval(parse(text=paste0("quote({", .val2, "})")))
                        return(.rxToSE(.val2, envir))
                    } else if (inherits(.val2, "numeric") || inherits(.val2, "integer")){
                        return(sprintf("%s", .val2))
                    }
                }
            }
        }
        return(.rxToSE(x, envir))
    }
    return(.rxToSE(eval(parse(text=paste0("quote({", x, "})"))), envir))
}
##'@rdname
##'@export
.rxToSE <- function(x, envir=NULL){
    .cnst <- names(.rxSEreserved)
    .isEnv <- inherits(envir, "rxS") || inherits(envir, "environment")
    if (is.name(x) || is.atomic(x)){
        .ret <- as.character(x)
        if (any(.ret == .cnst)){
            .ret <- paste0("rx_SymPy_Res_", .ret);
            if (.isEnv && is.name(x)){
                if (substr(x, 1, 1) != "."){
                    if (!exists(.ret, envir=envir))
                        assign(.ret, symengine::Symbol(.ret), envir=envir)
                }
            }
            return(.ret);
        } else {
            .ret0 <- .rxSEcnt[.ret];
            if (is.na(.ret0)){
                if (.isEnv && is.name(x)){
                    ## message(.ret)
                    if (substr(x, 1, 1) != "."){
                        if (!exists(.ret, envir=envir))
                            assign(.ret, symengine::Symbol(.ret), envir=envir)
                    }
                }
                return(.ret)
            }
            return(setNames(.ret0, NULL));
        }
    } else if (is.call(x)){
        if (identical(x[[1]], quote(`(`))){
            return(paste0("(", .rxToSE(x[[2]], envir=envir), ")"))
        } else if (identical(x[[1]], quote(`{`))){
            .x2 <- x[-1];
            return(paste(lapply(.x2, .rxToSE, envir=envir),
                         collapse="\n"));
        } else if (identical(x[[1]], quote(`*`)) ||
                   identical(x[[1]], quote(`^`)) ||
                   identical(x[[1]], quote(`+`)) ||
                   identical(x[[1]], quote(`-`)) ||
                   identical(x[[1]], quote(`/`))){
            if (length(x) == 3){
                if (identical(x[[1]], quote(`/`))){
                    .x2 <- x[[2]];
                    .x3 <- x[[3]];
                    ##df(%s)/dy(%s)
                    if (identical(.x2, quote(`d`)) &&
                        identical(.x3[[1]], quote(`dt`))){
                        .state <- .rxToSE(.x3[[2]], envir=envir);
                        return(paste0("rx__d_dt_", .state, "__"));
                    } else {
                        if (length(.x2) == 2 && length(.x3) == 2){
                            if (identical(.x2[[1]], quote(`df`)) &&
                                identical(.x3[[1]], quote(`dy`))){
                                .state <- .rxToSE(.x2[[2]], envir=envir);
                                .var <- .rxToSE(.x3[[2]], envir=envir);
                                return(paste0("rx__df_", .state,
                                              "_dy_", .var, "__"));
                            }
                        }
                        .ret <- paste0(.rxToSE(.x2, envir=envir),
                                      as.character(x[[1]]),
                                      .rxToSE(.x3, envir=envir))
                    }
                } else {
                    .ret <- paste0(.rxToSE(x[[2]], envir=envir),
                                  as.character(x[[1]]),
                                  .rxToSE(x[[3]], envir=envir))
                }
                return(.ret)
            } else {
                ## Unary Operators
                return(paste(as.character(x[[1]]),
                             .rxToSE(x[[2]], envir=envir)))
            }
        } else if (identical(x[[1]], quote(`=`)) ||
                   identical(x[[1]], quote(`<-`)) ||
                   identical(x[[1]], quote(`~`))){
            .var <- .rxToSE(x[[2]], envir=envir);
            if (inherits(x[[3]], "numeric")){
                if (.isEnv){
                    assign(.var, x[[3]], envir=envir)
                }
            } else {
                .expr <- paste0("with(envir,",
                                     .rxToSE(x[[3]],
                                             envir=envir), ")")
                .expr <- eval(parse(text=.expr));
                if (.isEnv){
                    assign(.var, .expr, envir=envir)
                    if (regexpr(rex::rex(or(.regRate,
                                            .regDur,
                                            .regLag,
                                            .regF,
                                            regIni0,
                                            regDDt)), .var) != -1){
                        .rx <- paste0(rxFromSE(.var), "=",
                                      rxFromSE(.expr));
                        if (regexpr(rex::rex(or(regSens, regSensEtaTheta)),
                                    .var) != -1){
                            assign("..sens0", c(envir$..sens0, .rx),
                               envir=envir)
                        } else {
                            assign("..ddt", c(envir$..ddt, .rx),
                               envir=envir)
                        }
                    } else if (regexpr(rex::rex(or(regDfDy,
                                                   regDfDyTh)), .var) != -1){
                        .rx <- paste0(rxFromSE(.var), "=",
                                      rxFromSE(.expr))
                        assign("..jac0", c(envir$..jac0, .rx),
                               envir=envir)
                    } else if (!identical(x[[1]], quote(`~`))) {
                        .rx <- paste0(rxFromSE(.var), "=",
                                      rxFromSE(.expr))
                        if (!any(.var == c("rx_pred_", "rx_r_"))){
                            assign("..lhs", c(envir$..lhs, .rx),
                                   envir=envir)
                        }

                    }
                }
            }
        } else if (identical(x[[1]], quote(`[`))){
            .type <- toupper(as.character(x[[2]]))
            if (any(.type == c("THETA", "ETA"))){
                if (is.numeric(x[[3]])){
                    .num <- x[[3]]
                    if (round(.num) == .num){
                        if (.num > 0){
                            if (.isEnv){
                                if (.type == "THETA"){
                                    if (exists("..maxTheta", envir=envir)){
                                        .m <- get("..maxTheta", envir=envir)
                                    } else {
                                        .m <- 0;
                                    }
                                    .m <- max(.m, .num);
                                    .funs <- envir$..curCall
                                    .funs <- .funs[.funs != ""]
                                    if (length(.funs) > 0){
                                        .funs <- paste(.funs, collapse=".")
                                        envir$..extraTheta[[.funs]] <-
                                            unique(c(envir$..extraTheta[[.funs]],
                                                     .num))
                                    }
                                    assign("..maxTheta",.m, envir=envir)
                                } else {
                                    if (exists("..maxEta", envir=envir)){
                                        .m <- get("..maxEta", envir=envir)
                                    } else {
                                        .m <- 0;
                                    }
                                    .m <- max(.m, .num);
                                    .funs <- envir$..curCall
                                    .funs <- .funs[.funs != ""]
                                    if (length(.funs) > 0){
                                        .funs <- paste(.funs, collapse=".")
                                        envir$..extraEta[[.funs]] <-
                                        unique(c(envir$..extraEta[[.funs]],
                                                 .num))
                                    }
                                    assign("..maxEta", .m, envir=envir)
                                }
                            }
                            return(paste0(.type, "_", .num, "_"))
                        } else {
                            stop("Only THETA[#] or ETA[#] are supported")
                        }
                    } else {
                        stop("Only THETA[#] or ETA[#] are supported")
                    }
                } else {
                    stop("Only THETA[#] or ETA[#] are supported")
                }
            } else {
                stop("Only THETA[#] or ETA[#] are supported")
            }
        } else if (identical(x[[1]], quote(`psigamma`))){
            if (length(x == 3)){
                if (.isEnv){
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "psigamma")
                }
                .a <- .rxToSE(x[[2]], envir=envir);
                .b <- .rxToSE(x[[3]], envir=envir);
                if (.isEnv) envir$..curCall <- .lastCall
                return(paste0("polygamma(", .b, ",", .a, ")"))
            } else {
                stop("psigamma() takes 2 arguments");
            }
        } else if (identical(x[[1]], quote(`log1pmx`))){
            if (length(x == 2)){
                if (.isEnv){
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "log")
                }
                .a <- .rxToSE(x[[2]], envir=envir);
                if (.isEnv) envir$..curCall <- .lastCall
                return(paste0("(log(1+", .a, ")-(", .a, "))"))
            } else {
                stop("log1pmx() only takes 1 argument");
            }
        } else if (identical(x[[1]], quote(`choose`))){
            if (length(x) == 3){
                if (.isEnv){
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "gamma")
                }
                .n <- .rxToSE(x[[2]], envir=envir)
                .k <- .rxToSE(x[[3]], envir=envir)
                if (.isEnv) envir$..curCall <- .lastCall
                return(paste0("gamma(", .n, "+1)/(gamma(",
                              .k, "+1)*gamma(", .n, "-(", .k, ")+1))"));
            } else {
                stop("choose() takes 2 arguments")
            }
        } else if (identical(x[[1]], quote(`lchoose`))){
            if (length(x) == 3){
                if (.isEnv){
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "lgamma")
                }
                .n <- .rxToSE(x[[2]], envir=envir)
                .k <- .rxToSE(x[[3]], envir=envir)
                if (.isEnv) envir$..curCall <- .lastCall
                return(paste0("(loggamma(", .n, "+1)-loggamma(", .k, "+1)-loggamma(", .n, "-(", .k, ")+1))"))
            } else {
                stop("lchoose() takes 2 arguments")
            }
        } else if (identical(x[[1]], quote(`transit`))){
            if (length(x) == 4){
                ##transit(n, mtt, bio)
                if (.isEnv){
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "lgamma")
                }
                .n <- .rxToSE(x[[2]], envir=envir);
                if (.isEnv) {
                    envir$..curCall <- .lastCall
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "log")
                }
                .mtt <- .rxToSE(x[[3]], envir=envir);
                .bio <- .rxToSE(x[[4]], envir=envir);
                if (.isEnv) envir$..curCall <- .lastCall
                return(paste0("exp(log((", .bio, ")*(podo))+log(",
                              .n, " + 1)-log(", .mtt, ")+(", .n,
                              ")*((log(", .n, "+1)-log(", .mtt,
                              "))+log(t))-((", .n, "+1)/(", .mtt,
                              "))*(t)-loggamma(1+", .n, "))"))
            } else if (length(x) == 3){
                if (.isEnv){
                    .lastCall <- envir$..curCall
                    envir$..curCall <- c(envir$..curCall, "lgamma")
                }
                .n <- .rxToSE(x[[2]], envir=envir);
                .mtt <- .rxToSE(x[[3]], envir=envir);
                if (.isEnv) envir$..curCall <- .lastCall
                return(paste0("exp(log(podo)+(log(", .n, "+1)-log(", .mtt, "))+(", .n, ")*((log(", .n, "+1)-log(", .mtt, "))+ log(t))-((", .n, " + 1)/(", .mtt, "))*(t)-loggamma(1+", .n, "))"))
            } else {
                stop("'transit' can only take 2-3 arguments");
            }
        } else {
            if (length(x[[1]]) == 1){
                .x1 <- as.character(x[[1]])
                .xc <- .rxSEsingle[[.x1]];
                if (!is.null(.xc)){
                    if (length(x) == 2){
                        if (.isEnv){
                            .lastCall <- envir$..curCall
                            envir$..curCall <- c(envir$..curCall, .xc[3])
                        }
                        .ret <- paste0(.xc[1], .rxToSE(x[[2]], envir=envir),
                                       .xc[2])
                        if (.isEnv) envir$..curCall <- .lastCall
                        return(.ret)
                    } else {
                        stop(sprintf("%s() only acceps 1 argument", .x1));
                    }
                }
                .xc <- .rxSEdouble[[.x1]];
                if (!is.null(.xc)){
                    .x1 <- as.character(x[[1]])
                    if (length(x) == 3){
                        .ret <- paste0(.xc[1], .rxToSE(x[[2]], envir=envir),
                                      .xc[2],
                                      .rxToSE(x[[3]], envir=envir),
                                      .xc[3]);
                        return(.ret)
                    } else {
                        stop(sprintf("%s() only acceps 2 arguments", .x1));
                    }
                }
            }
            if (.isEnv){
                .lastCall <- envir$..curCall
                envir$..curCall <- c(envir$..curCall,
                                     as.character(x[[1]]))
            }
            .ret0 <- c(list(as.character(x[[1]])), lapply(x[-1], .rxToSE, envir=envir));
            if (.isEnv) envir$..curCall <- .lastCall
            .SEeq <- c(.rxSEeq, .rxSEeqUsr)
            .nargs <- .SEeq[paste(.ret0[[1]])];
            if (!is.na(.nargs)){
                if (.nargs == length(.ret0) - 1){

                    .ret <- paste0(.ret0[[1]], "(")
                    .ret0 <- .ret0[-1];
                    .ret <- paste0(.ret, paste(unlist(.ret0), collapse=","), ")");
                    if (.ret == "exp(1)") return("E");
                    return(.ret)
                } else {
                    stop(sprintf("%s() takes %s arguments (has %s)",
                                 paste(.ret0[[1]]),
                                 .nargs, length(.ret0) - 1))
                }
            } else {
                .fun <- paste(.ret0[[1]])
                .ret0 <- .ret0[-1];
                .ret <- paste0("(", paste(unlist(.ret0), collapse=","), ")");
                if (.ret == "(0)"){
                    return(paste0("rx_", .fun, "_ini_0__"))
                } else if (any(.fun == c("cmt", "dvid"))){
                    return("")
                } else {
                    stop(sprintf("%s() not supported in RxODE", .fun));
                }
            }
        }
    } else {
        stop("Unsupported expression.");
    }
}

.rxUnXi <- function(x){
    gsub("_xi_([[:digit:]]*)", "rx_xi_\\1", x);
}


## 0L = error
## 1L = forward
## 2L = central
.rxFromNumDer <- 0L
.rxDelta <- (.Machine$double.eps)^(1/3)


##'@rdname rxToSE
##'@export
rxFromSE <- function(x, unknownDerivatives=c("forward", "central", "error")){
    .unknown <- c("central"=2L, "forward"=1L, "error"=0L)
    assignInMyNamespace(".rxFromNumDer", .unknown[match.arg(unknownDerivatives)])
    if (is(substitute(x),"character")){
        return(.rxFromSE(eval(parse(text=paste0("quote({", .rxUnXi(x), "})")))))
    } else if (is(substitute(x), "{")){
        x <- deparse(substitute(x));
        if (x[1] == "{"){
            x <- x[-1];
            x <- x[-length(x)];
        }
        x <- .rxUnXi(paste(x, collapse="\n"));
        .ret <- .rxFromSE(eval(parse(text=paste0("quote({", x, "})"))));
        return(.ret)
    } else {
        .xc <- as.character(substitute(x));
        x <- substitute(x);
        if (length(.xc == 1)){
            .found <- FALSE
            .frames <- seq(1, sys.nframe());
            .frames <- .frames[.frames != 0];
            for (.f in .frames){
                .env <- parent.frame(.f);
                if (exists(.xc, envir=.env)){
                    .val2 <- try(get(.xc, envir=.env), silent=TRUE);
                    if (inherits(.val2, "character")){
                        .val2 <- eval(parse(text=paste0("quote({", .rxUnXi(.val2), "})")))
                        .ret <- .rxFromSE(.val2)
                        return(.ret)
                    } else if (inherits(.val2, "numeric") || inherits(.val2, "integer")){
                        return(sprintf("%s", .val2))
                    } else {
                        if (!is.null(attr(class(.val2), "package"))){
                            if (attr(class(.val2), "package") == "symengine"){
                                .val2 <- eval(parse(text=paste0("quote({", .rxUnXi(as.character(.val2)), "})")))
                                .ret <- .rxFromSE(.val2)
                                return(.ret)
                            }
                        }
                    }
                }
            }
        }
        x <- eval(parse(text=paste("quote(", .rxUnXi(paste(deparse(x), collapse=" ")), ")")));
        .ret <- .rxFromSE(x)
        return(.ret)
    }
    x <- .rxUnXi(x);
    .ret <- .rxFromSE(eval(parse(text=paste0("quote({", x, "})"))))
    return(.ret)
}

.stripP <- function(x){
    if (is.call(x)){
        if (length(x) == 1){
            return(x)
        } else if (identical(x[[1]], quote(`(`))){
            return(.stripP(x[[2]]));
        } else {
            return(x)
        }
    } else {
        return(x);
    }
}

.rxM1rmF <- function(x){
    .env <- new.env(parent=emptyenv())
    .env$found <- FALSE
    .f <- function(x, envir){
        if (is.call(x)){
            if (identical(x[[1]], quote(`*`))){
                if (length(x) == 3){
                    .x2 <- as.character(x[[2]])
                    .x3 <- as.character(x[[3]])
                    if (length(.x3) == 1){
                        if (.x3 == "pi"){
                            envir$found <- TRUE
                            return(.rxFromSE(.stripP(x[[2]])));
                        }
                    }
                    if (length(.x2) == 1){
                        if (.x2 == "pi"){
                            envir$found <- TRUE
                            return(.rxFromSE(.stripP(x[[3]])));
                        }
                    }
                    return(paste0(.f(x[[2]], envir), "*",
                                  .f(x[[3]], envir)))
                } else {
                    return(.rxFromSE(x));
                }
            } else if (identical(x[[1]], quote(`/`))){
                .x2 <- as.character(x[[2]])
                .x3 <- as.character(x[[3]])
                if (length(.x2) == 1){
                    if (.x2 == "pi"){
                        envir$found <- TRUE
                        x[[2]] <- 1;
                        return(.rxFromSE(x));
                    }
                    return(paste0(.f(x[[2]], envir), "/",
                                  .f(x[[3]], envir)))
                }
            } else {
                return(.rxFromSE(x))
            }
        } else {
            return(.rxFromSE(x))
        }
    }
    .ret <- .f(x, .env)
    return(list(.ret, .env$found))
}

.rxP1rmF <- function(x){
    .env <- new.env(parent=emptyenv())
    .env$found <- FALSE
    .f <- function(x, envir){
        if (is.call(x)){
            if (identical(x[[1]], quote(`+`))){
                if (length(x) == 3){
                    if (inherits(x[[3]], "numeric")) {
                        if (x[[3]] == 1){
                            envir$found <- TRUE
                            return(.rxFromSE(.stripP(x[[2]])));
                        } else {
                            return(paste0(.f(x[[2]], envir),
                                          "+", as.character(x[[3]])));
                        }
                    }
                    if (inherits(x[[2]], "numeric")) {
                        if (x[[2]] == 1){
                            envir$found <- TRUE
                            return(.rxFromSE(.stripP(x[[3]])));
                        } else {
                            return(paste0(as.character(x[[2]]), "+",
                                          .f(x[[3]], envir)));
                        }
                    }
                    return(paste0(.f(x[[2]], envir), "+",
                                  .f(x[[3]], envir)))
                } else {
                    return(.rxFromSE(x));
                }
            } else if (identical(x[[1]], quote(`-`))){
                if (length(x) == 3){
                    if (inherits(x[[2]], "numeric")){
                        if (x[[2]] == 1){
                            x <- x[-2];
                            envir$found <- TRUE
                            return(.rxFromSE(x));
                        } else {
                            return(.rxFromSE(x));
                        }
                    } else {
                        return(.rxFromSE(x));
                    }
                } else {
                    return(.rxFromSE(x));
                }
            } else {
                return(.rxFromSE(x))
            }
        } else {
            return(.rxFromSE(x))
        }
    }
    .ret <- .f(x, .env)
    return(list(.ret, .env$found))
}

##'@export
##'@rdname rxToSE
.rxFromSE <- function(x){
    .cnst <- setNames(names(.rxSEreserved),
                      paste0("rx_SymPy_Res_", names(.rxSEreserved)))
    if (is.name(x) || is.atomic(x)){
        .ret <- as.character(x)
        .ret0 <- .cnst[.ret];
        if (!is.na(.ret0)){
            return(.ret0)
        }
        .ret0 <- .rxSEreserved[[.ret]];
        if (!is.null(.ret0)){
            if (is.character(.ret0)){
                return(.ret0);
            } else if (is.numeric(.ret0)){
                return(sprintf("%.16f", .ret0));
            }
        }
        if (.ret == "E") return("M_E");
        return(sub(.regRate, "rate(\\1)",
               sub(.regDur,"dur(\\1)",
               sub(.regLag,"alag(\\1)",
               sub(.regF, "f(\\1)",
               sub(regIni0, "\\1(0)",
               sub(regDfDy, "df(\\1)/dy(\\2)",
               sub(regDfDyTh, "df(\\1)/dy(\\2[\\3])",
               sub(regDDt, "d/dt(\\1)",
               sub(rex::rex(start, regThEt, end),
                       "\\1[\\2]", .ret))))))))))
    } else if (is.call(x)){
        if (identical(x[[1]], quote(`(`))){
            return(paste0("(", .rxFromSE(x[[2]]), ")"))
        } else if (identical(x[[1]], quote(`{`))){
            .x2 <- x[-1];
            return(paste(lapply(.x2, .rxFromSE),
                         collapse="\n"));
        } else if (identical(x[[1]], quote(`*`)) ||
                   identical(x[[1]], quote(`^`)) ||
                   identical(x[[1]], quote(`+`)) ||
                   identical(x[[1]], quote(`-`)) ||
                   identical(x[[1]], quote(`/`))){
            ## Unary Operators
            if (length(x) == 3){
                .x2 <- x[[2]];
                .x3 <- x[[3]];
                .ret <- paste0(.rxFromSE(.x2),
                              as.character(x[[1]]),
                              .rxFromSE(.x3))
                ## FIXME parsing to figure out if *2 or *0.5 *0.4 is in
                ## expression
                if (any(.ret == c("pi*2", "2*pi")))
                    return("M_2PI");
                if (any(.ret == c("pi/2", "pi*0.5", "0.5*pi")))
                    return("M_PI_2");
                if (any(.ret == c("pi/4", "pi*0.25", "0.25*pi")))
                    return("M_PI_4");
                if (.ret == "1/pi") return("M_1_PI")
                if (.ret == "2/pi") return("M_2_PI")
                if (any(.ret == c("(M_2_PI)^0.5", "(M_2_PI)^(1/2)",
                                  "M_2_PI^0.5", "M_2_PI^(1/2)")))
                    return("M_SQRT_2dPI")

                if (any(.ret == c("(pi)^0.5", "(pi)^(1/2)",
                                  "pi^0.5", "pi^(1/2)")))
                    return("M_SQRT_PI")
                if (.ret == "log(2)/log(10)") return("M_LOG10_2")
                if (.ret == "1/log(10)") return("M_LOG10E")
                if (.ret == "1/log(2)") return("M_LOG2E")
                if (any(.ret == c("2/M_SQRT_PI", "2/(M_SQRT_PI)")))
                    return("M_2_SQRTPI")
                if (any(.ret == c("1/sqrt(M_2PI)",
                                  "1/(M_2PI^0.5)", "1/(M_2PI^(1/2))",
                                  "1/((M_2PI)^0.5)", "1/((M_2PI)^(1/2))")))
                    return("M_1_SQRT_2PI")
                return(.ret)
            } else {
                .ret <- paste0(as.character(x[[1]]),
                         .rxFromSE(x[[2]]))

                return(.ret)
            }
        } else if (identical(x[[1]], quote(`=`)) ||
                   identical(x[[1]], quote(`<-`)) ||
                   identical(x[[1]], quote(`~`))){
            .var <- .rxFromSE(x[[2]]);
            .val <- .rxFromSE(x[[3]]);
        } else if (identical(x[[1]], quote(`[`))){
            stop("[...] expressions not supported")
        } else if (identical(x[[1]], quote(`polygamma`))){
            if (length(x == 3)){
                .a <- .rxFromSE(x[[2]]);
                .b <- .rxFromSE(x[[3]]);
                if (.a == "0"){
                    return(paste0("digamma(", .b, ")"))
                } else if (.a == "1"){
                    return(paste0("trigamma(", .b, ")"))
                } else if (.a == "2"){
                    return(paste0("tetragamma(", .b, ")"))
                } else if (.a == "3"){
                    return(paste0("pentagamma(", .b, ")"))
                } else {
                    return(paste0("psigamma(", .b, ",",
                                  .a, ")"))
                }
            } else {
                stop("polygamma() takes 2 arguments");
            }
        }  else {
            if (length(x[[1]]) == 1){
                .x1 <- as.character(x[[1]])
                .xc <- .SEsingle[[.x1]];
                if (!is.null(.xc)){
                    if (length(x) == 2){
                        .x2 <- x[[2]];
                        if (length(.x2) != 1){
                            if (identical(.x2[[1]], quote(`+`))){
                                .tmp0 <- .SE1p[.x1];
                                if (!is.na(.tmp0)){
                                    .ret <- .rxP1rmF(.x2)
                                    if (.ret[[2]]){
                                        .r1 <- .ret[[1]];
                                        return(paste0(.tmp0, "(",
                                                      .r1,
                                                      ")"))
                                    }
                                }
                                return(paste0(.xc[1], .ret[[1]], .xc[2]))
                            }
                        }

                        return(paste0(.xc[1], .rxFromSE(x[[2]]), .xc[2]))
                    } else {
                        stop(sprintf("%s() only acceps 1 argument", .x1));
                    }
                }
                .xc <- .SEdouble[[.x1]];
                if (!is.null(.xc)){
                    if (length(x) == 3){
                        .x1 <- .rxFromSE(x[[1]])
                        return(paste0(.xc[1], .rxFromSE(x[[2]]), .xc[2],
                                      .rxFromSE(x[[3]]),
                                      .xc[3]))
                    } else {
                        stop(sprintf("%s() only acceps 2 arguments", .x1));
                    }
                }
            }
            if (length(x) == 2){
                if (identical(x[[1]], quote(`log`))){
                    if (length(x[[2]]) == 3){
                        if (identical(x[[2]][[1]], quote(`beta`))){
                            .tmp <- x[[2]];
                            .tmp[[1]] <- quote(`lbeta`);
                            return(.rxFromSE(.tmp))
                        }
                    }
                }
            }
            .ret0 <- lapply(lapply(x, .stripP), .rxFromSE)
            .SEeq <- c(.rxSEeq, .rxSEeqUsr)
            .nargs <- .SEeq[paste(.ret0[[1]])];
            if (!is.na(.nargs)){
                if (.nargs == length(.ret0) - 1){
                    .x1 <- as.character(.ret0[[1]]);
                    .tmp0 <- .x1
                    if (.nargs == 1){
                        .tmp0 <- .SE1p[.x1];
                        .x2 <- x[[2]];
                        if (!is.na(.tmp0)){
                            .ret <- .rxP1rmF(.x2)
                            if (.ret[[2]]){
                                if (.tmp0 == "log1p"){
                                    .tmp <- eval(parse(text=paste0("quote(", .ret[[1]], ")")))
                                    if (length(.tmp) > 1){
                                        if (identical(.tmp[[1]], quote(`exp`))){
                                            .tmp <- .tmp[[-1]];
                                            .tmp0 <- "log1pexp"
                                            .ret[[1]] <- .rxFromSE(.tmp);
                                        }
                                    }
                                }
                                return(paste0(.tmp0, "(",
                                              .ret[[1]],
                                              ")"))
                            } else {
                                .ret <- paste0(.x1, "(",
                                              .ret[[1]],
                                              ")")
                                if (.ret == "log(2)") return("M_LN2");
                                if (.ret == "log(10)") return("M_LN10");
                                if (.ret == "log(M_SQRT_PI)")
                                    return("M_LN_SQRT_PI")
                                if (any(.ret == c("log(sqrt(M_PI_2))",
                                              "log((M_PI_2)^(1/2))",
                                              "log((M_PI_2)^0.5)",
                                              "log(M_PI_2^(1/2))",
                                              "log(M_PI_2^0.5)")))
                                    return("M_LN_SQRT_PId2")
                                if (any(.ret == c("log(sqrt(M_2PI))",
                                                  "log((M_2PI)^0.5)",
                                                  "log((M_2PI)^(1/2))",
                                                  "log(M_2PI^0.5)",
                                                  "log(M_2PI^(1/2))")))
                                    return("M_LN_SQRT_2PI")
                                return(.ret)
                            }
                        }
                        .tmp0 <- .SE1m[.x1];
                        if (!is.na(.tmp0)){
                            .ret <- .rxM1rmF(.x2)
                            if (.ret[[2]]){
                                return(paste0(.tmp0, "(",
                                              .ret[[1]],
                                              ")"))
                            } else {
                                return(paste0(.x1, "(",
                                              .ret[[1]],
                                              ")"))
                            }
                        }
                        .tmp0 <- .x1
                    }
                    .ret <- paste0(.tmp0, "(")
                    .ret0 <- .ret0[-1];
                    .ret <- paste0(.ret, paste(unlist(.ret0), collapse=","),
                                   ")");
                    if (.ret == "exp(1)") return("M_E");
                    if (.ret == "sin(pi)") return("0");
                    if (.ret == "cos(pi)") return("1");
                    if (.ret == "tan(pi)") return("0");
                    if (.ret == "sqrt(3)") return("M_SQRT_3");
                    if (.ret == "sqrt(2)") return("M_SQRT_2");
                    if (.ret == "sqrt(32)") return("M_SQRT_32");
                    if (.ret == "sqrt(pi)") return("M_SQRT_PI");
                    if (.ret == "sqrt(M_2_PI)")
                        return("M_SQRT_2dPI");
                    return(.ret)
                } else {
                    stop(sprintf("%s() takes %s arguments",
                                 paste(.ret0[[1]]),
                                 .nargs))
                }
            } else if (identical(x[[1]], quote(`Derivative`))){
                if (length(x) == 3){
                    .fun <- as.character(x[[2]])
                    .var <- .rxFromSE(x[[3]])
                    .args <- .fun[-1];
                    .args <- lapply(.args, .rxFromSE)
                    .with <- which(.var == .args)
                    .errD <- function(force=FALSE){
                        if (!force && .rxFromNumDer != 0L){
                            ## Can calculate forward or central
                            ## difference instead.
                            ## Warn
                            if (.rxFromNumDer == 1L){
                                ## Forward
                                .a1 <- .args
                                .fn <- .fun[1];
                                .a1[.with] <- paste0("(", .a1[.with], ")+", .rxDelta)
                                .a2 <- .args
                                return(paste0("(", .fn, "(", paste0(.a1, collapse=",") ,")-",
                                       .fn, "(", paste0(.a2, collapse=","), "))/", .rxDelta))
                            } else if (.rxFromNumDer == 2L) {
                                ## Central
                                .a1 <- .args
                                .fn <- .fun[1];
                                .a1[.with] <- paste0(.a1[.with], "-", (0.5 * .rxDelta))
                                .a2 <- .args
                                .a2[.with] <- paste0(.a2[.with], "+", (0.5 * .rxDelta))
                                return(paste0("(", .fn, "(", paste0(.a1, collapse=",") ,")-",
                                       .fn, "(", paste0(.a2, collapse=","), "))/", .rxDelta))
                            } else {
                                stop("Only forward and central differences are supported")
                            }
                        } else {
                            stop(sprintf("Cannot figure out the `%s` derivative with respect to `%s`", .fun[1], .var[1]))
                        }
                    }
                    if (length(.with) != 1){
                        .errD(force=TRUE)
                    }
                    if (exists(.fun[1], envir=.rxD)){
                        .funLst <- get(.fun[1], envir=.rxD);
                        if (length(.funLst) < .with){
                            return(.errD())
                        }
                        .derFun <- .funLst[[.with]];
                        if (is.null(.derFun)){
                            return(.errD())
                        }
                        return(do.call(.derFun, as.list(.args)));
                    } else {
                        if (.rxFromNumDer == 0L){
                            stop(sprintf("RxODE/symengine does not know how to take a derivative of `%s`", .fun[1]))
                        } else {
                            return(.errD())
                        }
                    }
                } else {
                    stop("Derivative() conversion only takes one function and one argument")
                }
            } else if (identical(x[[1]], quote(`Subs`))){
                .fun <- eval(parse(text=paste0("quote(", .rxFromSE(x[[2]]), ")")))
                .what <- .stripP(x[[3]])
                .with <- .stripP(x[[4]])
                .subs <- function(x){
                    if (identical(x, .what)){
                        return(.with);
                    } else if (is.call(x)) {
                        as.call(lapply(x, .subs))
                    } else if (is.pairlist(x)) {
                        as.pairlist(lapply(x, .subs))
                    } else {
                        return(x);
                    }
                }
                .ret <- .subs(.fun)
                return(.rxFromSE(.ret))
            } else {
                stop(sprintf("%s() not supported in symengine->RxODE", paste(.ret0[[1]])));
            }
        }
    } else {
        stop("Unsupported expression.");
    }
}

.rxFunction <- function(name){
    .f <- function(...){1};
    body(.f) <- bquote(return(symengine::FunctionSymbol(.(name), unlist(list(...)))))
    return(.f)
}

##' Load a model into a symengine environment
##'
##' @param x RxODE object
##' @return RxODE/symengine environment
##' @author Matthew Fidler
##' @export
rxS <- function(x){
    .cnst <- names(.rxSEreserved)
    .env <- new.env(parent = loadNamespace("symengine"))
    .env$..mv <- rxModelVars(x);
    .env$..jac0 <- c();
    .env$..ddt <- c();
    .env$..sens0 <- c();
    .env$..lhs <- c();
    .env$rxTBS <- .rxFunction("rxTBS")
    .env$rxTBSd <- .rxFunction("rxTBSd")
    .env$rxTBSd2 <- .rxFunction("rxTBSd2")
    .env$solveLinB <- .rxFunction("solveLinB");
    for (.f in c("rxEq", "rxNeq", "rxGeq", "rxLeq", "rxLt", "rxGt", "rxAnd", "rxOr", "rxNot"))
        assign(.f, .rxFunction(.f), envir=.env)
    .env$..polygamma <- symengine::S("polygamma(_rx_a, _rx_b)");
    .env$..a <- symengine::Symbol("_rx_a");
    .env$..b <- symengine::Symbol("_rx_b");
    .env$..s0 <- symengine::S("0")
    .env$..extraTheta <- list()
    .env$..extraEta <- list()
    .env$..curCall <- character(0)
    .env$polygamma <- function(a, b){
        symengine::subs(symengine::subs(..polygamma, ..a, a), ..b,  b)
    }
    .env$loggamma <- function(a){
        lgamma(a)
    }
    .pars <- c(rxParams(x), rxState(x),
               "podo", "t", "time", "tlast", "rx1c", "rx__PTR__");
    ## default lambda/yj values
    .env$rx_lambda_ <- symengine::S("1")
    .env$rx_yj_ <- symengine::S("2")

    sapply(names(.rxSEeqUsr), function(x){
        assign(.rxFunction(x), x, envir=.env);
    })
    ## EulerGamma=0.57721566490153286060651209008240243104215933593992
    ## S("I")
    ## S("pi")
    ## S("E")
    ## S("") # EulerGamma
    ## S("Catalan") = 0.915965594177219015054603514932384110774
    ## S("GoldenRatio") = 1+sqrt(5)/2
    ## S("inf")
    ## S("nan")
    sapply(.pars, function(x){
        if (any(.cnst == x)){
            .tmp <- paste0("rx_SymPy_Res_", x);
            assign(.tmp, symengine::Symbol(.tmp), envir=.env)
        } else {
            .tmp <- rxToSE(x);
            assign(.tmp, symengine::Symbol(.tmp), envir=.env)
            assign(x, symengine::Symbol(x), envir=.env)
        }
    })
    .expr <- eval(parse(text=paste0("quote({",rxNorm(x),"})")));
    .ret <- .rxToSE(.expr, .env)
    class(.env) <- "rxS";
    return(.env)
}
