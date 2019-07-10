.rxSEsingle <- list("gammafn"=c("gamma(", ")"),
                    "lgammafn"=c("loggamma(", ")"),
                    "lgamma"=c("loggamma(", ")"),
                    "digamma"=c("polygamma(0,", ")"),
                    "trigamma"=c("polygamma(1,", ")"),
                    "tetragamma"=c("polygamma(2, ", ")"),
                    "pentagamma"=c("polygamma(3, ", ")"),
                    "cospi"=c("cos(pi * (", "))"),
                    "sinpi"=c("sin(pi * (", "))"),
                    "tanpi"=c("tan(pi * (", "))"),
                    "log1p"=c("log(1+", ")"),
                    "expm1"=c("(exp(", ")-1)"),
                    "factorial"=c("gamma(", "+1)"),
                    "lgamma1p"=c("loggamma(", "+1)"),
                    "expm1"=c("(exp(", ")-1)"),
                    "log10"=c("log(", ")/log(10)"),
                    "log2"=c("log(", ")/log(2)")
                    )

.rxSEdouble <- list("pow"=c("(", ")^(", ")"),
                    "R_pow"=c("(", ")^(", ")"),
                    "R_pow_di"=c("(", ")^(", ")"),
                    "Rx_pow_di"=c("(", ")^(", ")"),
                    "Rx_pow"=c("(", ")^(", ")"),
                    "lbeta"=c("log(beta(", ",", "))")
                    )
## atan2
.rxSEeq <- c("acos"=1, "acosh"=1, "asin"=1, "atan"=1,
             "atanh"=1, "beta"=2,
             "cos"=1, "cosh"=1, "erf"=1, "erfc"=1,
             "exp"=1, "gamma"=1, "sin"=1, "sinh"=1,
             "sqrt"=1, "tan"=1, "tanh"=1, "log"=1, "abs"=1, "asinh"=1,
             "rxTBS"=3, "rxTBSd"=3, "rxTBSd2"=3)

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
              "M_LOG10_2" = "log10(2)",
              "M_LOG2E" = "1/log(2)",
              "M_LOG10E" = "log10(E)",
              "M_LN2" = "log(2)",
              "M_LN10" = "log(10)");
## "rxTBS", "rxTBSd"

##' RxODE to symengine.R
##'
##' @param x expression
##' @param envir \code{NULL}
##' @return
##' @export
##' @author Matthew L. Fidler
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
    .cnst <- c("e", "E");
    if (is.name(x) || is.atomic(x)){
        .ret <- as.character(x)
        if (any(.ret == .cnst)){
            return(paste0("rx_SymPy_Res_", .ret));
        } else {
            .ret0 <- .rxSEcnt[.ret];
            if (is.na(.ret0)) return(.ret)
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
                        return(paste0(.rxToSE(.x2, envir=envir),
                                      as.character(x[[1]]),
                                      .rxToSE(.x3, envir=envir)))
                    }
                } else {
                    return(paste0(.rxToSE(x[[2]], envir=envir),
                                  as.character(x[[1]]),
                                  .rxToSE(x[[3]], envir=envir)));
                }
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
            } else {
                .expr <- eval(parse(text=paste0("with(.env, ",
                                                .rxToSE(x[[3]],
                                                       envir=envir), ")")));
                if (inherits(envir, "environment")){
                    assign(.var, .expr, envir=envir)
                }
            }
        } else if (identical(x[[1]], quote(`[`))){
            .type <- toupper(as.character(x[[2]]))
            if (any(.type == c("THETA", "ETA"))){
                if (is.numeric(x[[3]])){
                    .num <- x[[3]]
                    if (round(.num) == .num){
                        if (.num > 0){
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
                .a <- .rxToSE(x[[2]], envir=envir);
                .b <- .rxToSE(x[[3]], envir=envir);
                return(paste0("polygamma(", .b, ",", .a, "))"))
            } else {
                stop("psigamma() takes 2 arguments");
            }
        } else if (identical(x[[1]], quote(`log1pmx`))){
            if (length(x == 2)){
                .a <- .rxToSE(x[[2]], envir=envir);
                return(paste0("(log(1+", .a, ")-(", .a, "))"))
            } else {
                stop("log1pmx() only takes 1 argument");
            }
        } else if (identical(x[[1]], quote(`choose`))){
            if (length(x) == 3){
                .n <- .rxToSE(x[[2]], envir=envir)
                .k <- .rxToSE(x[[3]], envir=envir)
                return(paste0("gamma(", .n, "+1)/(gamma(",
                              .k, "+1)*gamma(", .n, "-(", .k, ")+1))"));
            } else {
                stop("choose() takes 2 arguments")
            }
        } else if (identical(x[[1]], quote(`lchoose`))){
            if (length(x) == 3){
                .n <- .rxToSE(x[[2]], envir=envir)
                .k <- .rxToSE(x[[3]], envir=envir)
                return(paste0("(loggamma(", .n, "+1)-loggamma(", .k, "+1)-loggamma(", .n, "-(", .k, ")+1))"))
            } else {
                stop("lchoose() takes 2 arguments")
            }
        } else if (identical(x[[1]], quote(`transit`))){
            if (length(x) == 4){
                ##transit(n, mtt, bio)
                .n <- as.character(x[[2]]);
                .mtt <- as.character(x[[3]]);
                .bio <- as.character(x[[4]]);
                return(paste0("exp(log((", .bio, ")*(podo))+log(",
                              .n, " + 1)-log(", .mtt, ")+(", .n,
                              ")*((log(", .n, "+1)-log(", .mtt,
                              "))+log(t))-((", .n, "+1)/(", .mtt,
                              "))*(t)-loggamma(1+", .n, "))"))
            } else if (length(x) == 3){
                .n <- as.character(x[[2]]);
                .mtt <- as.character(x[[3]]);
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
                        return(paste0(.xc[1], as.character(x[[2]]), .xc[2]))
                    } else {
                        stop(sprintf("%s() only acceps 1 argument", .x1));
                    }
                }
                .xc <- .rxSEdouble[[.x1]];
                if (!is.null(.xc)){
                    if (length(x) == 3){
                        .x1 <- as.character(x[[1]])
                        return(paste0(.xc[1], as.character(x[[2]]), .xc[2],
                                      as.character(x[[3]]),
                                      .xc[3]))
                    } else {
                        stop(sprintf("%s() only acceps 2 arguments", .x1));
                    }
                }
            }
            .ret0 <- lapply(x, .rxToSE, envir=envir);
            .nargs <- .rxSEeq[paste(.ret0[[1]])];
            if (!is.na(.nargs)){
                if (.nargs == length(.ret0) - 1){
                    .ret <- paste0(.ret0[[1]], "(")
                    .ret0 <- .ret0[-1];
                    .ret <- paste0(.ret, paste(unlist(.ret0), collapse=","), ")");
                    return(.ret)
                } else {
                    stop(sprintf("%s() takes %s arguments",
                                 paste(.ret0[[1]]),
                                 .nargs))
                }
            } else {
                stop(sprintf("%s() not supported in RxODE", paste(.ret0[[1]])));
            }
        }
    } else {
        stop("Unsupported expression.");
    }
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param x
##' @return
##' @author Matthew Fidler
##' @export
rxS <- function(x){
    .symengine <- (attr(x, "package") == "symengine")
    .cnst <- c("e", "E");
    if (.symengine){
    } else {
        .env <- new.env(parent = loadNamespace("symengine"))
        .env$..mv <- rxModelVars(x);
        .pars <- c(rxParams(x), rxState(x), "podo", "t", "time", "tlast", "rx_lambda_", "rx_yj_", "rx1c");
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
                assign(.tmp, S(.tmp), envir=.env)
            } else {
                assign(x, S(x), envir=.env)
            }
        })
        .expr <- eval(parse(text=paste0("quote({",rxNorm(mod),"})")));
        .ret <- rxToSE(.expr, .env)
        class(.env) <- "rxS";
        return(.env)
    }
}
