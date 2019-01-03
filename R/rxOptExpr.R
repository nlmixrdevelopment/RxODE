.rxOptFn <- function(fn){
    force(fn);
    function(...){
        .ret <- paste0(fn, "(", paste(unlist(list(...)), collapse=", "), ")")
        .new <- .rxOptEnv$.rep[[.ret]];
        if (!is.null(.new)){
            if (length(.rxOptEnv$.exclude) != 1) .rxOptEnv$.exclude <- ""
            if (.new == .rxOptEnv$.exclude){
                return(.ret)
            } else {
                ## message(sprintf("%s->%s", .ret, .new))
                ## print(.rxOptEnv$.rep)
                return(.new)
            }
        }
        if(is.null(.rxOptEnv$.list[[.ret]])){
            .rxOptEnv$.list[[.ret]] <- 1L;
        } else {
            .rxOptEnv$.list[[.ret]] <- .rxOptEnv$.list[[.ret]] + 1L;
        }
        return(.ret)
    }
}
.rxOptBin <- function(sep) {
    force(sep)
    function(e1, e2) {
        if (missing(e2)){
            if (sep == "+"){
                .ret <- paste0(e1)
            } else {
                .ret <- paste0(gsub(" ", "", sep), e1)
            }
        } else {
            .ret <- paste0(e1, sep, e2)
        }
        .new <- .rxOptEnv$.rep[[.ret]];
        if (!is.null(.new)){
            if (length(.rxOptEnv$.exclude) != 1) .rxOptEnv$.exclude <- ""
            if (.new == .rxOptEnv$.exclude){
                return(.ret)
            } else {
                ## message(sprintf("%s->%s", .ret, .new))
                ## print(.rxOptEnv$.rep)
                return(.new)
            }
        }
        if(is.null(.rxOptEnv$.list[[.ret]])){
            .rxOptEnv$.list[[.ret]] <- 1L;
        } else {
            .rxOptEnv$.list[[.ret]] <- .rxOptEnv$.list[[.ret]] + 1L;
        }
        return(.ret)
    }
}

.rxOptEnv <- new.env(parent = emptyenv())
.rxOptEnv$"^" <- .rxOptBin("^")
.rxOptEnv$"**" <- .rxOptBin("^")

.rxOptEnv[["*"]] <- .rxOptBin("*");

.rxOptEnv[["/"]] <- .rxOptBin("/");

.rxOptEnv[["+"]] <- .rxOptBin("+");

.rxOptEnv[["-"]] <- .rxOptBin("-");

.rxOptEnv$"[" <- function(name, val){
    .n <- toupper(name)
    .err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(.n == c("THETA", "ETA")) && is.numeric(val)){
        if (round(val) == val && val > 0){
            return(sprintf("%s[%s]", .n, val));
        } else {
            stop(.err);
        }
    } else {
        stop(.err)
    }
}
.rxOptEnv$"{" <- function(...){
    return(sprintf("{\n%s\n}", paste(unlist(list(...)), collapse="\n")))
}
.rxOptEnv$"[" <- function(name, val){
    .n <- toupper(name)
    .err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(.n == c("THETA", "ETA")) && is.numeric(val)){
        if (round(val) == val && val > 0){
            return(sprintf("%s[%s]", .n, val));
        } else {
            stop(err);
        }
    } else {
        stop(err)
    }
}

.rxOptEnv$"(" <- unaryOp("(", ")")

.rxOptEnv$.list <- list();
.rxOptEnv$.rep <- list();
.rxOptEnv$.exclude <- "";

.rxOptGetEnv <- function(expr){
    ## Known functions
    .calls <- allCalls(expr)
    .callList <- setNames(lapply(.calls, .rxOptFn), .calls)
    .callEnv <- list2env(.callList);
    .currEnv <- cloneEnv(.rxOptEnv, .callEnv);
    .names <- allNames(expr)
    .n1 <- .names;
    .n2 <-.names;
    .symbolList <- setNames(as.list(.n2), .n1);
    .symbolEnv <- list2env(.symbolList, parent=.currEnv);
    return(.symbolEnv)
}
.rxOptExpr <- function(x){
    return(eval(x, .rxOptGetEnv(x)))
}
##' Optimize RxODE for computer evaluation
##'
##' This optimizes RxODE code for computer evaluation by only
##' calculating redundant expressions once.
##'
##' @param x RxODE model that can be access by rxNorm
##' @return Optimized RxODE model text.  The order and type lhs and
##'     state variables is maintained while the evaluation is sped up.
##'     While parameters names are maintained, their order may be
##'     modified.
##' @author Matthew L. Fidler
##' @export
rxOptExpr <- function(x){
    .rxOptEnv$.list <- list();
    .rxOptEnv$.rep <- list();
    .rxOptEnv$.exclude <- "";
    .lines <- strsplit(rxNorm(x), "\n")[[1]];
    .f <- function(line, onlyRet=FALSE){
        .silent <- (regexpr("[~]", line) != -1)
        .l2 <- strsplit(line, "[=~]")[[1]]
        if (length(.l2) == 2){
            if (regexpr(rex::rex("if", any_spaces, "("), .l2[1]) != -1) return(line);
            .l1 <- gsub(" +", "", .l2[1])
            .rxOptEnv$.exclude <- .l1;
            .ret <- eval(parse(text=sprintf(".rxOptExpr(quote(%s))", gsub(";$", "",.l2[2]))));
            if (.silent){
                return(paste0(.l1, " ~ ", .ret))
            } else {
                if (onlyRet) return(.ret)
                return(paste0(.l1, " = ", .ret))
            }
        } else {
            return(line)
        }
    }
    .ret <- sapply(.lines, .f)
    .rxOptEnv$.list <- .rxOptEnv$.list[which(unlist(.rxOptEnv$.list) > 1L)];
    .exprs <- names(.rxOptEnv$.list)[order(nchar(names(.rxOptEnv$.list)))];
    .exprs <- .exprs[regexpr(rex::rex(start, regNum, end), .exprs, perl=TRUE) == -1]
    .exprs <- .exprs[regexpr(rex::rex(start, or("THETA[", "ETA["), any_numbers, "]", end), .exprs, perl=TRUE) == -1]
    .exprsI <- sprintf("rx_expr_%03d", seq_along(.exprs))
    if (length(.exprs) > 0){
        .rxOptEnv$.rep <- setNames(as.list(.exprsI), .exprs)
        .rxOptEnv$.exclude <- ""
        .lines0 <- sapply(.lines, .f);
        while(TRUE){
            .exprs0 <- sapply(paste("rx_dummy = ", .exprs), .f, onlyRet=TRUE)
            .w <- which(.exprs0 != .exprsI);
            .exprs[.w] <- .exprs0[.w]
            .rxOptEnv$.rep <- setNames(as.list(.exprsI), .exprs)
            .rxOptEnv$.exclude <- ""
            .linesf <- sapply(.lines0, .f);
            if (all(.linesf == .lines0)){
                break;
            } else {
                .lines0 <- .linesf
            }
        }
        .ret <- c(sprintf("rx_expr_%03d ~ %s", seq_along(.exprs), .exprs), .linesf);
        return(paste(.ret, collapse="\n"))
    } else {
        return(x)
    }
}

################################################################################

.rxPowEnv0 <- new.env(parent = emptyenv())
.rxPowEnv0$.theta <- c()
.rxPowEnv0$"[" <- function(name, val){
    .n <- toupper(name)
    .err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(.n == c("THETA", "ETA")) && is.numeric(val)){
        if (.n == "THETA") .rxPowEnv0$.theta <- sort(unique(c(val, .rxPowEnv0$.theta)))
        if (round(val) == val && val > 0){
            return(sprintf("%s[%s]", .n, val));
        } else {
            stop(.err);
        }
    } else {
        stop(.err)
    }
}

.rxPowGetEnv0 <- function(expr){
    ## Known functions
    .calls <- allCalls(expr)
    .callList <- setNames(lapply(.calls, functionOp), .calls)
    .callEnv <- list2env(.callList);
    .currEnv <- cloneEnv(.rxPowEnv0, .callEnv);
    .names <- allNames(expr)
    .n1 <- .names;
    .n2 <- .names;
    .symbolList <- setNames(as.list(.n2), .n1);
    .symbolEnv <- list2env(.symbolList, parent=.currEnv);
    return(.symbolEnv)
}

.rxPowExpr0 <- function(x){
    return(eval(x, .rxPowGetEnv0(x)))
}

.rxGetTheta <- function(x){
    .rxPowEnv0$.theta <- c()
    eval(parse(text=sprintf(".rxPowExpr0(quote(%s))", x)))
    return(.rxPowEnv0$.theta)
}

.rxPowBin <- function(sep) {
    force(sep)
    function(e1, e2) {
        if (is.null(.rxPowEnv$.powTheta)){
            .rxPowEnv$.powTheta <- .rxGetTheta(e2);
        } else {
            .rxPowEnv$.powTheta <- sort(unique(c(.rxGetTheta(e2), .rxGetTheta(e2))));
        }
        .ret <- paste0(e1, sep, e2)
        return(.ret)
    }
}

.rxPowEnv <- new.env(parent = emptyenv())
.rxPowEnv$.powTheta <- c()
.rxPowEnv$"^" <- .rxPowBin("^")
.rxPowEnv$"**" <- .rxPowBin("^")
.rxPowEnv$Rx_pow <- .rxPowBin("^")
.rxPowEnv$Rx_pow_di <- .rxPowBin("^")
.rxPowEnv$pow <- .rxPowBin("^")

.rxPowEnv$factorial <- function(x){
    .rxPowEnv$.factorial <- sort(unique(c(.rxGetTheta(x), .rxPowEnv$.factorial)))
    return(paste0("factorial(", x, ")"))
}

.rxPowEnv$gamma <- function(x){
    .rxPowEnv$.gamma <- sort(unique(c(.rxGetTheta(x), .rxPowEnv$.gamma)))
    return(paste0("gamma(", x, ")"))
}

.rxPowEnv$log <- function(x){
    .rxPowEnv$.log <- sort(unique(c(.rxGetTheta(x), .rxPowEnv$.log)))
    return(paste0("log(", x, ")"))
}

.rxPowEnv$sin <- function(x){
    .rxPowEnv$.sin <- sort(unique(c(.rxGetTheta(x), .rxPowEnv$.sin)))
    return(paste0("sin(", x, ")"))
}

.rxPowEnv$cos <- function(x){
    .rxPowEnv$.cos <- sort(unique(c(.rxGetTheta(x), .rxPowEnv$.cos)))
    return(paste0("cos(", x, ")"))
}

.rxPowEnv$tan <- function(x){
    .rxPowEnv$.tan <- sort(unique(c(.rxGetTheta(x), .rxPowEnv$.tan)))
    return(paste0("tan(", x, ")"))
}

.rxPowEnv[["*"]] <- binaryOp("*");

.rxPowEnv[["/"]] <- binaryOp("/");

.rxPowEnv[["+"]] <- binaryOp("+");

.rxPowEnv[["-"]] <- binaryOp("-");

.rxPowEnv$"[" <- function(name, val){
    .n <- toupper(name)
    .err <- "RxODE only supports THETA[#] and ETA[#] numbers."
    if (any(.n == c("THETA", "ETA")) && is.numeric(val)){
        if (round(val) == val && val > 0){
            return(sprintf("%s[%s]", .n, val));
        } else {
            stop(.err);
        }
    } else {
        stop(.err)
    }
}

.rxPowEnv$"{" <- function(...){
    return(sprintf("{\n%s\n}", paste(unlist(list(...)), collapse="\n")))
}

.rxPowEnv$"(" <- unaryOp("(", ")")

.rxPowGetEnv <- function(expr){
    ## Known functions
    .calls <- allCalls(expr)
    .callList <- setNames(lapply(.calls, functionOp), .calls)
    .callEnv <- list2env(.callList);
    .currEnv <- cloneEnv(.rxPowEnv, .callEnv);
    .names <- allNames(expr)
    .n1 <- .names;
    .n2 <- .names;
    .symbolList <- setNames(as.list(.n2), .n1);
    .symbolEnv <- list2env(.symbolList, parent=.currEnv);
    return(.symbolEnv)
}

.rxPowExpr <- function(x){
    return(eval(x, .rxPowGetEnv(x)))
}

##' Find power THETAs for appropriate scaling
##'
##' @param x RxODE model that can be access by rxNorm
##' @return THETA numbers of x^theta
##' @author Matthew L. Fidler
.rxFindPow <- function(x){
    .rxPowEnv$.powTheta <- c();
    .rxPowEnv$.factorial <- c();
    .rxPowEnv$.gamma <- c();
    .rxPowEnv$.log <- c();
    .rxPowEnv$.sin <- c();
    .rxPowEnv$.cos <- c();
    .rxPowEnv$.tan <- c();
    .lines <- strsplit(rxNorm(x), "\n")[[1]];
    .f <- function(line){
        .l2 <- strsplit(line, "[=~]")[[1]]
        if (length(.l2) == 2){
            if (regexpr(rex::rex("if", any_spaces, "("), .l2[1]) != -1) return(line);
            .l1 <- gsub(" +", "", .l2[1])
            .ret <- eval(parse(text=sprintf(".rxPowExpr(quote(%s))", gsub(";$", "",.l2[2]))));
            return(paste0(.l1, " = ", .ret))
        } else {
            return(line)
        }
    }
    sapply(.lines, .f)
    return(list(powTheta=.rxPowEnv$.powTheta,
                factorial=.rxPowEnv$.factorial,
                gamma=.rxPowEnv$.gamma,
                log=.rxPowEnv$.log,
                sin=.rxPowEnv$.sin,
                cos=.rxPowEnv$.cos,
                tan=.rxPowEnv$.tan));
}
