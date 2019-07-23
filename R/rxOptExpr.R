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

.rxOptEnv[["&&"]] <- .rxOptBin("&&");
.rxOptEnv[["||"]] <- .rxOptBin("||");
.rxOptEnv[["|"]] <- .rxOptBin("|");
.rxOptEnv[["&"]] <- .rxOptBin("&");
.rxOptEnv[["=="]] <- .rxOptBin("==");
.rxOptEnv[["<="]] <- .rxOptBin("<=");
.rxOptEnv[[">="]] <- .rxOptBin(">=");
.rxOptEnv[["<"]] <- .rxOptBin("<");
.rxOptEnv[[">"]] <- .rxOptBin(">");
.rxOptEnv[["!="]] <- .rxOptBin("!=");

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

..rxOpt <- function(x, onlyRet=FALSE, progress=FALSE){
    if (is.name(x) || is.atomic(x)){
        return(as.character(x))
    } else if (is.call(x)){
        .x2 <- x[-1];
        if (identical(x[[1]], quote(`{`))){
            ## Probably can implement a progress bar here.
            if (progress){
                rxProgress(length(.x2));
                on.exit({rxProgressAbort("Stopped optimizing duplicate expressions")});
                .ret <- unlist(lapply(.x2, function(x){
                    rxTick()
                    ..rxOpt(x)
                }))
                rxProgressStop();
                return(.ret)
            } else {
                return(unlist(lapply(.x2, ..rxOpt, onlyRet=onlyRet)))
            }
        } else if (identical(x[[1]], quote(`*`)) ||
                   identical(x[[1]], quote(`^`)) ||
                   identical(x[[1]], quote(`+`)) ||
                   identical(x[[1]], quote(`-`)) ||
                   identical(x[[1]], quote(`/`)) ||
                   identical(x[[1]], quote(`==`)) ||
                   identical(x[[1]], quote(`>=`)) ||
                   identical(x[[1]], quote(`<=`)) ||
                   identical(x[[1]], quote(`>`)) ||
                   identical(x[[1]], quote(`<`)) ||
                   identical(x[[1]], quote(`!=`)) ||
                   identical(x[[1]], quote(`&&`)) ||
                   identical(x[[1]], quote(`||`)) ||
                   identical(x[[1]], quote(`&`)) ||
                   identical(x[[1]], quote(`|`))){
            if (length(x) == 3){
                return(paste0(..rxOpt(x[[2]], onlyRet=onlyRet), as.character(x[[1]]),
                              ..rxOpt(x[[3]], onlyRet=onlyRet)));
            } else {
                ## Unary Operators
                return(paste(as.character(x[[1]]),
                             ..rxOpt(x[[2]], onlyRet=onlyRet)))
            }
        } else if (identical(x[[1]], quote(`~`)) ||
                   identical(x[[1]], quote(`=`)) ||
                   identical(x[[1]], quote(`<-`))){
            if (onlyRet){
                return(.rxOptExpr(x[[3]]));
            } else {
                return(paste0(..rxOpt(x[[2]]),
                              ifelse(identical(x[[1]], quote(`<-`)), "=", as.character(x[[1]])),
                              .rxOptExpr(x[[3]])))
            }
        } else if (identical(x[[1]], quote(`[`))){
            return(paste0(..rxOpt(x[[2]]), "[", ..rxOpt(x[[3]]), "]"));
        } else {
            .ret0 <- lapply(x, ..rxOpt);
            .ret <- paste0(.ret0[[1]], "(")
            if (.ret == "((") .ret <- "("
            .ret0 <- .ret0[-1];
            .ret <- paste0(.ret, paste(unlist(.ret0), collapse=", "), ")");
            return(.ret)
        }
    }
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
rxOptExpr <- function(x, msg="model"){
    .rxOptEnv$.list <- list();
    .rxOptEnv$.rep <- list();
    .rxOptEnv$.exclude <- "";
    message(paste0("Finding duplicate expressions in ", msg, "..."))
    .p <- eval(parse(text=paste0("quote({", x, "})")))
    .lines <- ..rxOpt(.p, progress=TRUE)
    .rxOptEnv$.list <- .rxOptEnv$.list[which(unlist(.rxOptEnv$.list) > 1L)];
    .exprs <- names(.rxOptEnv$.list)[order(nchar(names(.rxOptEnv$.list)))];
    .exprs <- .exprs[regexpr(rex::rex(start, regNum, end), .exprs, perl=TRUE) == -1]
    .exprs <- .exprs[regexpr(rex::rex(start, or("THETA[", "ETA["), any_numbers, "]", end), .exprs, perl=TRUE) == -1]
    if (length(.exprs) > 0){
        .rp <- rxOptRep_(.exprs)
        .rxOptEnv$.rep <- as.list(.rp[[1]]);
        .rxOptEnv$.exclude <- ""
        message(paste0("Optimizing duplicate expressions in ", msg, "..."))
        return(paste(c(.rp[[2]], ..rxOpt(.p, progress=TRUE)), collapse="\n"))
    } else {
        return(x)
    }
}
