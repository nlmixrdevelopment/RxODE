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
    bound <- .getBound(x, parent.frame(2));
    cat(cli::rule(center=crayon::bold(paste0("EventTable with ",x$nobs+x$ndose, " records"))), "\n");
    ## sprintf(" with %s records%s:\n", x$nobs+x$ndose,
    ##                                           ifelse(x$maxId==1, "", sprintf(" (%d IDs)", abs(x$maxId))))
    .units <- x$.units;
    .maxId <- length(x$IDs)
    if (.maxId !=1){
        cat(sprintf("   %s individuals\n",
                .maxId))
    }
    cat(sprintf("   %s dosing records (see %s$%s())\n",
                x$ndose, crayon::yellow(bound), crayon::blue("get.dosing")))
    cat(sprintf("   %s observation times (see %s$%s())\n",
                x$nobs, crayon::yellow(bound), crayon::blue("get.sampling")))
    if (x$nobs!=0 | x$ndose!=0){
        cat(cli::rule(crayon::bold(paste0("First part of ",crayon::yellow(bound),":"))), "\n");
        print(dplyr::as.tbl(data.frame(x)[, x$show, drop = FALSE]));
    }
    cat(cli::rule(), "\n");
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

##' Add dosing to eventTable
##'
##' This adds a dosing event to the event table.  This is provided for
##' piping syntax through magrittr
##'
##' @param eventTable eventTable object
##' @param dose numeric scalar, dose amount in \code{amount.units};
##' @param nbr.doses integer, number of doses;
##' @param dosing.interval required numeric scalar, time between doses
##'      in \code{time.units}, defaults to 24 of \code{time.units="hours"};
##' @param dosing.to integer, compartment the dose goes into
##'        (first compartment by default);
##' @param rate for infusions, the rate of infusion (default
##'            is \code{NULL}, for bolus dosing;
##' @param amount.units optional string indicating the dosing units.
##'           Defaults to \code{NA} to indicate as per the original \code{EventTable}
##'           definition.
##' @param start.time required dosing start time;
##' @param do.sampling logical, should observation sampling records
##'            be added at the dosing times? Defaults to \code{FALSE}.
##' @param time.units optional string indicating the time units.
##'           Defaults to \code{"hours"} to indicate as per the original \code{EventTable} definition.
##' @param ... Other parameters (ignored)
##' @return eventTable with updated dosing (note the event table will be updated anyway)
##' @author Matthew L. Fidler
##' @seealso \code{\link{eventTable}}, \code{\link{RxODE}}
##' @export
add.dosing <- function(eventTable, dose, nbr.doses = 1L, dosing.interval = 24, dosing.to = 1L, rate = NULL, amount.units = NA_character_, start.time = 0.0, do.sampling = FALSE, time.units = NA_character_, ...) {
    .lst <- list(dose=dose,
                 nbr.doses=nbr.doses,
                 start.time=start.time,
                 do.sampling=do.sampling);
    if (!is.na(amount.units)) .lst$amount.units <- amount.units;
    if (!is.na(time.units)) .lst$time.units <- time.units;
    if (dosing.to != 1) .lst$dosing.to <- dosing.to
    if (!is.null(rate)) .lst$rate <- rate;
    .lst$dosing.interval <- dosing.interval;
    .Call(`_RxODE_et_`, .lst, eventTable);
}

##' Add sampling to eventTable
##'
##' This adds a dosing event to the event table.  This is provided for
##' piping syntax through magrittr
##'
##' @param eventTable An eventTable object
##' @param time a vector of time values (in \code{time.units}).
###' @param time.units an optional string specifying the time
##'     units. Defaults to the units specified when the
##'     \code{EventTable} was initialized.
##' @return eventTable with updated sampling.  (Note the event table
##'     will be updated even if you don't reassign the eventTable)
##' @author Matthew L. Fidler
##' @seealso \code{\link{eventTable}}, \code{\link{RxODE}}
##' @export
add.sampling <- function(eventTable, time, time.units = NA){
    .lst <- list(time=time);
    if (!is.na(time.units)) .lst$time.units <- time.units;
    .Call(`_RxODE_et_`, .lst, eventTable);
}


##' Create an event table object
##'
##' Initializes an object of class \sQuote{EventTable} with methods for
##' adding and querying dosing and observation records
##'
##' @param amount.units string denoting the amount dosing units, e.g.,
##'     \dQuote{mg}, \dQuote{ug}. Default to \code{NA} to denote
##'     unspecified units.  It could also be a solved RxODE object.  In
##'     that case, eventTable(obj) returns the eventTable that was used
##'     to solve the RxODE object.
##'
##' @param time.units string denoting the time units, e.g.,
##'     \dQuote{hours}, \dQuote{days}. Default to \code{"hours"}.
##'
##' An \code{eventTable} is an object that consists of a data.frame
##' storing ordered time-stamped events of an (unspecified) PK/PD
##' dynamic system, units (strings) for dosing and time records, plus a
##' list of functions to add and extract event records.
##'
##' Currently, events can be of two types: dosing events that represent
##' inputs to the system and sampling time events that represent
##' observations of the system with \sQuote{amount.units} and
##' \sQuote{time.units}, respectively. In the future, additional events
##' may include resetting of state variables (compartments), for
##' instance, to indicate time after \dQuote{wash-out}, etc.
##'
##' @return A modified data.frame with the following accessible functions:
##'
##' \item{get.EventTable}{returns the current event table.}
##'
##' \item{add.dosing}{adds dosing records to the event table.
##'
##' Its arguments are
##'
##'   \code{dose}: numeric scalar, dose amount in \code{amount.units};
##'
##'   \code{nbr.doses}: integer, number of doses;
##'
##'   \code{dosing.interval}: required numeric scalar, time between doses
##'      in \code{time.units}, defaults to 24 of \code{time.units="hours"};
##'
##'        \code{dosing.to}: integer, compartment the dose goes into
##'        (first compartment by default);
##'
##'        \code{rate}: for infusions, the rate of infusion (default
##'            is \code{NULL}, for bolus dosing;
##'
##'        \code{start.time}: required dosing start time;
##'
##'        \code{do.sampling}: logical, should observation sampling records
##'            be added at the dosing times? Defaults to \code{FALSE}.
##'
##'        \code{amount.units}: optional string indicating the dosing units.
##'           Defaults to \code{NA} to indicate as per the original \code{EventTable}
##'           definition.
##'
##'        \code{time.units}: optional string indicating the time units.
##'           Defaults to \code{"hours"} to indicate as per the original \code{EventTable} definition.
##'     }
##'
##'    \item{get.dosing}{returns a data.frame of dosing records.}
##'
##'    \item{clear.dosing}{clears or deletes all dosing from event table}
##'
##'    \item{add.sampling}{adds sampling time observation records to the
##'        event table. Its arguments are
##'
##'        \code{time} a vector of time values (in \code{time.units}).
##'
##'        \code{time.units} an optional string specifying the time
##'        units. Defaults to the units specified when the \code{EventTable}
##'        was initialized.
##'
##'        % TODO: should add.sampling() have similar calling sequence
##'        % as add.dosing()?
##'        %\code{sampling.interval}: scalar, time between samples.
##'        %\code{start.time}: scalar, starting observation time.
##'        %\code{end.time}: scalar, end observation time.
##'    }
##'
##'    \item{get.sampling}{returns a data.frame of sampled observation
##'        records.}
##'
##'    \item{clear.sampling}{removes all sampling from event table.}
##'
##'    \item{get.obs.rec}{returns a logical vector indicating
##'        whether each event record represents an observation or not.}
##'
##'    \item{get.nobs}{returns the number of observation (not dosing) records.}
##'
##'    \item{get.units}{returns a two-element character vector with the
##'        dosing and time units, respectively.}
##'
##'    \item{copy}{makes a copy of the current event table. To create
##'        a copy of an event table object use \code{qd2 <- qd$copy()}.}
##'
##'    \item{expand}{Expands the event table for multi-subject solving.
##'    This is done by qd$expand(400) for a 400 subject data expansion}
##'
##' @author Matthew Fidler, Melissa Hallow and Wenping Wang
##'
##' @seealso \code{\link{et}}, \code{\link{RxODE}}
##'
##' @examples
##' # create dosing and observation (sampling) events
##' # QD 50mg dosing, 5 days followed by 25mg 5 days
##' #
##' qd <- eventTable(amount.units = "mg", time.units = "days")
##' #
##' qd$add.dosing(dose=50, nbr.doses=5, dosing.interval = 1, do.sampling=FALSE)
##' #
##' # sample the system's drug amounts hourly the first day, then every 12 hours
##' # for the next 4 days
##' qd$add.sampling(seq(from = 0, to = 1, by = 1/24))
##' qd$add.sampling(seq(from = 1, to = 5, by = 12/24))
##' #
##' #print(qd$get.dosing())     # table of dosing records
##' print(qd$get.nobs())   # number of observation (not dosing) records
##' #
##' # BID dosing, 5 days
##' bid <- eventTable("mg", "days")  # only dosing
##' bid$add.dosing(dose=10000, nbr.doses=2*5,
##'                dosing.interval = 12, do.sampling=FALSE)
##' #
##' # Use the copy() method to create a copy (clone) of an existing
##' # event table (simple assignments just create a new reference to
##' # the same event table object (closure)).
##' #
##' bid.ext <- bid$copy()      # three-day extension for a 2nd cohort
##' bid.ext$add.dosing(dose = 5000, nbr.doses = 2*3,
##'                    start.time = 120, dosing.interval = 12, do.sampling = FALSE)
##'
##' # You can also use the Piping operator to create a table
##'
##' qd2 <- eventTable(amount.units="mg", time.units="days") %>%
##'     add.dosing(dose=50, nbr.doses=5, dosing.interval=1, do.sampling=FALSE) %>%
##'     add.sampling(seq(from=0, to=1, by=1 / 24)) %>%
##'     add.sampling(seq(from=1, to=5, by=12 / 24))
##' #print(qd2$get.dosing())     # table of dosing records
##' print(qd2$get.nobs())   # number of observation (not dosing) records
##'
##' # Note that piping with %>% will update the original table.
##'
##' qd3 <- qd2 %>% add.sampling(seq(from=5, to=10, by=6 / 24))
##' print(qd2$get.nobs())
##' print(qd3$get.nobs())
##' @keywords models data
##' @concept ordinary differential equations
##' @concept Nonlinear regression
##' @concept Pharmacokinetics (PK)
#' @concept Pharmacodynamics (PD)
#' @export
eventTable <- function(amount.units = NA, time.units = "hours"){
    et(amount.units = ifelse(is.na(amount.units), "", amount.units),
       time.units = ifelse(is.na(time.units), "",time.units))
}

##' Sequence of event tables
##'
##' @param ...
##' @param handleSamples can be "clear", "use"
##' @details
##'
##' @return A new event table
##' @author Matthew L Fidler
##' @export
etSeq <- function(...,handleSamples=c("clear", "use"),handleWait=c("smartAddIi", "alwaysAddII")){
    ## etSeq_(List ets, bool clearSampling=clearSampling);
    .sampleIx <- c(clear=0L,use=1L);
    .waitIx <- c(smartAddIi=0L, alwaysAddII=1L)
    .Call(`_RxODE_etSeq_`, list(...), setNames(.sampleIx[match.arg(handleSamples)],NULL),
          setNames(.waitIx[match.arg(handleWait)],NULL),
          0L, TRUE, character(0),logical(0),FALSE);
}

##'@export
seq.rxEt <- function(...){
    do.call(etSeq,list(...));
}

##'@export
c.rxEt <- function(...){
    do.call(etSeq,list(...));
}

##' Repeat an RxODE event table
##'
##'
##' @param times Number of times to repeat the event table
##' @param length.out Invalid with RxODE event tables, will throw an
##'     error if used.
##' @param each Invalid with RxODE event tables, will throw an error
##'     if used.
##' @param n The number of times to repeat the event table.  Overrides
##'     \code{times}.
##' @param wait Waiting time between each repeated event table.  By
##'     default there is no waiting, or wait=0
##' @param id IDs to expand/remove in the event table before repeating.
##' @inheritParams etSeq
##' @param ...
##' @return An event table of repeated events
##' @author Matthew Fidler
##' @export
etRep <- function(x, times=1, length.out=NA, each=NA, n=NULL, wait=0, handleSamples=c("clear", "use"), id=integer(0)){
    if (!is.null(n)){
        times <- n;
    }
    .sampleIx <- c(clear=0L,use=1L);
    if (!is.na(length.out)) stop("'length.out' makes no sense with event tables");
    if (!is.na(each)) stop("'each' makes no sense with event tables");
    .Call(`_RxODE_etRep_`, x, as.integer(times),
          as.double(wait), as.integer(id), setNames(.sampleIx[match.arg(handleSamples)],NULL),
          setNames(.waitIx[match.arg(handleWait)],NULL))
}

##'@export
rep.rxEt <- function(x, ...){
    do.call(etRep,list(x=x,...));
}


##' @importFrom magrittr %>%
##' @export
magrittr::`%>%`
