##' Prune branches from RxODE
##'
##' This prunse branches from RxODE.
##'
##' @param x RxODE model that can be accessed by rxNorm
##' @return Pruned RxODE model text.  All logical expressions are
##'     collapsed to rx_lgl_# and are at the beginning of the file.
##' @author Matthew L. Fidler
##' @export
rxPrune <- function(x){
    .f <- function(x, envir=parent.frame()){
        if (is.name(x) || is.atomic(x)){
            return(as.character(x))
        } else if (is.call(x)){
            if (identical(x[[1]], quote(`if`))){
                .if <- envir$.if;
                .if[length(.if) + 1] <- .f(x[[2]], envir=envir);
                envir$.if <- .if
                .x2 <- x[-(1:2)]
                envir$.else <- c();
                if (length(.x2) == 2){
                    .ret1 <- .f(.x2[[1]], envir=envir);
                    .if[length(.if)] <- paste0("1-(", .if[length(.if)], ")");
                    envir$.if <- .if
                    envir$.else <- findLhs(eval(parse(text=paste0("quote({", .ret1, "})"))));
                    .ret2 <- .f(.x2[[2]], envir=envir);
                    .ret <- paste0(.ret1, "\n", .ret2);
                } else if (length(.x2) == 1){
                    .ret <- .f(.x2[[1]], envir=envir);
                }
                .if <- .if[-length(.if)];
                return(.ret)
            } else if (identical(x[[1]], quote(`(`))){
                return(paste0("(", .f(x[[2]], envir=envir), ")"))
            } else if (identical(x[[1]], quote(`{`))){
                .x2 <- x[-1];
                return(paste(lapply(.x2, .f, envir=envir), collapse="\n"));
            } else if (identical(x[[1]], quote(`==`)) ||
                identical(x[[1]], quote(`>=`)) ||
                identical(x[[1]], quote(`<=`)) ||
                identical(x[[1]], quote(`>`)) ||
                identical(x[[1]], quote(`<`)) ||
                identical(x[[1]], quote(`!=`)) ||
                identical(x[[1]], quote(`&&`)) ||
                identical(x[[1]], quote(`||`)) ||
                identical(x[[1]], quote(`&`)) ||
                identical(x[[1]], quote(`|`))){
                .ret <- paste0(.f(x[[2]], envir=envir), as.character(x[[1]]),
                               .f(x[[3]], envir=envir))
                return(.ret)
                ## if (!any(envir$.logic == .ret)){
                ##     envir$.logic[length(envir$.logic) + 1] <- .ret
                ## }
                ## return(paste0("rx_lgl_", which(envir$.logic == .ret)));
            } else if (identical(x[[1]], quote(`=`)) ||
                       identical(x[[1]], quote(`<-`)) ||
                       identical(x[[1]], quote(`~`))){
                if (length(envir$.if > 0)){
                    .if <- paste(paste0("(", envir$.if, ")"), collapse="*");
                    .ret0 <- paste0(.f(x[[2]], envir=envir), as.character(x[[1]]), .if, "*(",
                              .f(x[[3]], envir=envir), ")", ifelse(any(envir$.else == as.character(x[[2]])), paste0("+", as.character(x[[2]])), ""))
                    .ret <- "";
                    ## if (length(envir$.logic) != length(envir$.logic2)){
                    ##     if (length(envir$.logic2) == 0){
                    ##         .ret <- paste0(paste(paste0("rx_lgl_", seq_along(envir$.logic), "~(", envir$.logic, ")"), collapse="\n"), "\n");
                    ##         envir$.logic2 <- envir$.logic;
                    ##     } else {
                    ##         .logic <- envir$.logic
                    ##         .logic <- .logic[-seq_along(envir$.logic2)];
                    ##         .ret <- paste0(paste(paste0("rx_lgl_", seq_along(.logic) + length(envir$.logic2), "~(", .logic, ")"), collapse="\n"), "\n");
                    ##         envir$.logic2 <- envir$.logic;
                    ##     }
                    ## }
                    return(paste0(.ret, .ret0));
                } else {
                    return(paste0(.f(x[[2]], envir=envir), as.character(x[[1]]),
                              .f(x[[3]], envir=envir)));
                }
            } else if (identical(x[[1]], quote(`*`)) ||
                       identical(x[[1]], quote(`^`)) ||
                       identical(x[[1]], quote(`+`)) ||
                       identical(x[[1]], quote(`-`)) ||
                       identical(x[[1]], quote(`/`))){
                if (length(x == 3)){
                    return(paste0(.f(x[[2]], envir=envir), as.character(x[[1]]),
                              .f(x[[3]], envir=envir)));
                } else {
                    ## Unary Operators
                }
            } else {
                return(as.call(lapply(x, .f, envir=envir)))
            }
        } else {
            ## is.pairlist OR is.atomic OR unknown...
            stop("Unsupported expression.");
        }
    }
    .env <- new.env(parent=emptyenv())
    .env$.logic <- c();
    .env$.logic2 <- c();
    .env$.if <- c();
    .ret <- .f(eval(parse(text=paste0("quote({", rxNorm(x), "})"))), envir=.env);
    return(.ret)
}
