##' Event Table Function
##'
##' @export
et <- function(...){
    UseMethod("et");
}

##'@rdname et
##'@export
et.default <- function(...){
    .Call(`_RxODE_et_`, list(...), list())
}

##' @export
`$.rxEt` <-  function(obj, arg, exact = FALSE){
    return(.Call(`_RxODE_etUpdate`, obj, arg, NULL))
}


##'@export
print.rxEt <- function(x,...){
    cat(sprintf("EventTable with %s records:\n", x$nobs+x$ndose))
    units <- x$.units;
    cat(sprintf("   %d dosing%s records\n",
                x$ndose,
                ifelse(is.na(units["dosing"]), "",
                       sprintf(" (in %s)", units["dosing"]))))
    cat(sprintf("   %d observation times%s\n",
                x$nobs,
                ifelse(is.na(units["time"]), "",
                       sprintf(" (in %s)", units["time"]))))
    if (x$nobs!=0 | x$ndose!=0){
        print(dplyr::as.tbl(data.frame(x)[, x$.show, drop = FALSE]));
    }
    invisible(x)
}

##'@export
str.rxEt <- function(object, ...){
    cat("rxEt methods and properties:\n");
    cat(" $ get.EventTable   :function ()\n");
    cat(" $ get.obs.rec      :function ()  \n");
    cat(" $ get.nobs         :function ()  \n");
    cat(" $ add.dosing       :function ()  \n");
    cat(" $ clear.dosing     :function ()  \n");
    cat(" $ get.dosing       :function ()  \n");
    cat(" $ add.sampling     :function ()  \n");
    cat(" $ clear.sampling   :function ()  \n");
    cat(" $ get.sampling     :function ()  \n");
    cat(" $ get.units        :function ()  \n");
    cat(" $ import.EventTable:function ()  \n");
    cat(" $ copy             :function ()  \n");
    return(invisible(NextMethod("str", ...)))
}

##'@export
print.rxHidden <- function(x,...){
    cat("\r");
}

##'@export
str.rxHidden <- function(object,...){
    cat("\r");
}

