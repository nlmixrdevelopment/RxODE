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
    cat(sprintf("EventTable with %s records%s:\n", x$nobs+x$ndose,
                ifelse(x$maxId==1, "", sprintf(" (%d IDs)", x$maxId))))
    .units <- x$.units;
    cat(sprintf("   %s dosing records\n",
                x$ndose))
    cat(sprintf("   %s observation times\n",
                x$nobs))
    if (x$nobs!=0 | x$ndose!=0){
        print(dplyr::as.tbl(data.frame(x)[, x$show, drop = FALSE]));
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

##'@export
set_units.rxEt <- function(x, value, ..., mode = units::units_options("set_units_mode")){
    if (missing(value))
        value <- units::unitless
    else if (mode == "symbols") {
        value <- substitute(value)
        if (is.numeric(value) && !identical(value, 1) && !identical(value, 1L))
            stop("The only valid number defining a unit is '1', signifying a unitless unit")
    }
    if (identical(value, units::unitless)){
        warning("Clearing both amount and time units; For more precise control use et(amountUnits=\"\") or et(timeUnits=\"\")")
        return(suppressWarnings({.Call(`_RxODE_et_`, list(amountUnits="", timeUnits=""), x)}))
    } else {
        if (!rxIs(value, "character")) value <- deparse(value);
        .tUnit <- units::set_units(1, "sec", mode="standard");
        .isTime <- try(units::set_units(units::set_units(1, value, mode="standard"), "sec"), silent=TRUE);
        if (inherits(.isTime, "try-error")){
            ## Amount
            return(.Call(`_RxODE_et_`, list(amountUnits=value), x))
        } else {
            ##
            return(.Call(`_RxODE_et_`, list(timeUnits=value), x));
        }
    }
}
