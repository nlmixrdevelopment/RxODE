.rxSEsingle <- list("gammafn"=c("gamma(", ")"),
                    "lgammafn"=c("loggamma(", ")"),
                    "lgamma"=c("loggamma(", ")"),
                    "digamma"=c("polygamma(0,", ")"),
                    "trigamma"=c("polygamma(1,", ")"),
                    "tetragamma"=c("polygamma(2, ", ")"),
                    "pentagamma"=c("polygamma(3, ", ")"),
                    "lbeta"=c("log(beta(", "))"),
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
                    "Rx_pow"=c("(", ")^(", ")")
                    )

.rxSEeq <- c("acos", "acosh", "asin", "atan", "atan2", "atanh", "beta",
             "cos", "cosh", "erf", "erfc", "exp", "gamma", "sin", "sinh",
             "sqrt", "tan", "tanh", "log", "abs", "asinh")
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
            return(.ret)
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
                        .state <- as.character(.x3[[2]]);
                        return(paste0("rx__d_dt_", .state, "__"));
                    } else {
                        if (length(.x2) == 2 && length(.x3) == 2){
                            if (identical(.x2[[1]], quote(`df`)) &&
                                identical(.x3[[1]], quote(`dy`))){
                                .state <- as.character(.x2[[2]]);
                                .var <- as.character(.x3[[2]]);
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
        } else if (identical(x[[1]], quote(`log1pmx`))){
            .a <- as.character(x[[2]]);
            return(paste0("(log(1+", .a, ")-(", .a, "))"))
        } else if (identical(x[[1]], quote(`choose`))){
            .n <- as.character(x[[2]])
            .k <- as.character(x[[3]])
            return(paste0("gamma(", .n, "+1)/(gamma(",
                          .k, "+1)*gamma(", .n, "-(", .k, ")+1))"));
        } else if (identical(x[[1]], quote(`lchoose`))){
            .n <- as.character(x[[2]])
            .k <- as.character(x[[3]])
            return(paste0("(loggamma(", .n, "+1)-loggamma(", .k, "+1)-loggamma(", .n, "-(", .k, ")+1))"))
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
                if (length(x) == 2){
                    .xc <- .rxSEsingle[[.x1]];
                    if (!is.null(.xc)){
                        return(paste0(.xc[1], as.character(x[[2]]), .xc[2]))
                    }
                } else if (length(x) == 3){
                    .x1 <- as.character(x[[1]])
                    .xc <- .rxSEdouble[[.x1]];
                    if (!is.null(.xc)){
                        return(paste0(.xc[1], as.character(x[[2]]), .xc[2],
                                      as.character(x[[3]]),
                                      .xc[3]))
                    }
                }
            }
            .ret0 <- lapply(x, .rxToSE, envir=envir);
            if (any(paste(.ret0[[1]]) == .rxSEeq)){
                .ret <- paste0(.ret0[[1]], "(")
                .ret0 <- .ret0[-1];
                .ret <- paste0(.ret, paste(unlist(.ret0), collapse=", "), ")");
                return(.ret)
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
