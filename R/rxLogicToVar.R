##' Change binary logical expressions to variables
##'
##' This is meant to be used with non-branching code.
##'
##' @param x RxODE model that can be accessed by rxNorm
##' @return list of logical expressions and replaced variables
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxLogicToVar <- function(x){
    .f <- function(x, envir=parent.frame()){
        if (is.name(x) || is.atomic(x)){
            return(as.character(x))
        } else if (is.call(x)){
            if (identical(x[[1]], quote(`(`))){
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
                .w <- which(.ret == envir$.logic);
                if (length(.w) != 1){
                    .w <- length(envir$.logic) + 1;
                    envir$.logic[.w] <- .ret
                }
                return(paste0("rx_lgl_", .w))
            } else if (identical(x[[1]], quote(`=`)) ||
                       identical(x[[1]], quote(`<-`)) ||
                       identical(x[[1]], quote(`~`))){
                return(paste0(.f(x[[2]], envir=envir), as.character(x[[1]]),
                              .f(x[[3]], envir=envir)));
            } else if (identical(x[[1]], quote(`*`)) ||
                       identical(x[[1]], quote(`^`)) ||
                       identical(x[[1]], quote(`+`)) ||
                       identical(x[[1]], quote(`-`)) ||
                       identical(x[[1]], quote(`/`))){
                if (length(x) == 3){
                    return(paste0(.f(x[[2]], envir=envir), as.character(x[[1]]),
                                  .f(x[[3]], envir=envir)));
                } else {
                    ## Unary Operators
                    return(paste(as.character(x[[1]]),
                                 .f(x[[2]], envir=envir)))
                }
            } else {
                .ret0 <- lapply(x, .f, envir=envir);
                .ret <- paste0(.ret0[[1]], "(")
                .ret0 <- .ret0[-1];
                .ret <- paste0(.ret, paste(unlist(.ret0), collapse=", "), ")");
                return(.ret)
            }
        } else { ## nocov start
            ## is.pairlist OR is.atomic OR unknown...
            stop("Unsupported expression.");
        } ## nocov end
    }
    .env <- new.env(parent=emptyenv())
    .env$.logic <- c();
    .ret <- .f(eval(parse(text=paste0("quote({", rxNorm(x), "})"))), envir=.env);
    return(list(.env$.logic, .ret))
}
##' Convert logical variable list and model to logical model
##'
##' @param variables Logical variables that will replace the rx_lgl_#
##'     in the model
##' @param model A model with rx_lgl_# information contained within
##' @return Model with the rx_lgl_# variables replaced with the
##'     logical variables in the data
##' @author Matthew Fidler
##' @keywords internal
##' @export
rxVarToLogic <- function(variables, model){
    .lgl <- paste0("rx_lgl_", seq_along(variables));
    .f <- function(x){
        if (is.name(x) || is.atomic(x)){
            .ret <- as.character(x);
            .w <- which(.ret == .lgl)
            if (length(.w) == 1){
                return(variables[.w]);
            } else {
                return(.ret)
            }
        } else if (is.call(x)){
            if (identical(x[[1]], quote(`(`))){
                return(paste0("(", .f(x[[2]]), ")"))
            } else if (identical(x[[1]], quote(`{`))){
                .x2 <- x[-1];
                return(paste(lapply(.x2, .f), collapse="\n"));
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

                return(paste0(.f(x[[2]]), as.character(x[[1]]),
                              .f(x[[3]])))

            } else if (identical(x[[1]], quote(`=`)) ||
                       identical(x[[1]], quote(`<-`)) ||
                       identical(x[[1]], quote(`~`))){
                return(paste0(.f(x[[2]]), as.character(x[[1]]),
                              .f(x[[3]])));
            } else if (identical(x[[1]], quote(`*`)) ||
                       identical(x[[1]], quote(`^`)) ||
                       identical(x[[1]], quote(`+`)) ||
                       identical(x[[1]], quote(`-`)) ||
                       identical(x[[1]], quote(`/`))){
                if (length(x) == 3){
                    return(paste0(.f(x[[2]]), as.character(x[[1]]),
                                  .f(x[[3]])));
                } else {
                    ## Unary Operators
                    return(paste(as.character(x[[1]]),
                                 .f(x[[2]])))
                }
            } else {
                .ret0 <- lapply(x, .f);
                .ret <- paste0(.ret0[[1]], "(")
                .ret0 <- .ret0[-1];
                .ret <- paste0(.ret, paste(unlist(.ret0), collapse=", "), ")");
                return(.ret)
            }
        } else { ## nocov start
            ## is.pairlist OR is.atomic OR unknown...
            stop("Unsupported expression.");
        } ## nocov end
    }
    return(.f(eval(parse(text=paste0("quote({", rxNorm(model), "})")))));
}
