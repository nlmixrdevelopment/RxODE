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
    return(.Call(`_RxODE_etUpdate`, obj, arg, NULL, exact))
}


##'@export
print.rxEt <- function(x,...){
    if (rxIs(x, "rxEt")){
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
    } else {
        print.data.frame(x)
    }
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
##' @param ... Other parameters passed to \code{\link{et}}.
##' @return eventTable with updated dosing (note the event table will be updated anyway)
##' @author Matthew L. Fidler
##' @seealso \code{\link{eventTable}}, \code{\link{RxODE}}
##' @export
add.dosing <- function(eventTable, dose, nbr.doses = 1L, dosing.interval = 24, dosing.to = 1L, rate = NULL, amount.units = NA_character_, start.time = 0.0, do.sampling = FALSE, time.units = NA_character_, ...) {
    .lst <- list(dose=dose,
                 nbr.doses=nbr.doses,
                 start.time=start.time,
                 do.sampling=do.sampling,
                 ...);
    if (!is.na(amount.units)) .lst$amount.units <- amount.units;
    if (!is.na(time.units)) .lst$time.units <- time.units;
    if (dosing.to != 1) .lst$dosing.to <- dosing.to
    if (!is.null(rate)) .lst$rate <- rate;
    if (nbr.doses > 1){
        .lst$dosing.interval <- dosing.interval;
    } else {
        .lst$dosing.interval <- 0.0;
    }

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
##'      \dQuote{mg}, \dQuote{ug}. Default to \code{NA} to denote
##'      unspecified units.  It could also be a solved RxODE object.  In
##'      that case, eventTable(obj) returns the eventTable that was used
##'      to solve the RxODE object.
##'
##' @param time.units string denoting the time units, e.g.,
##'      \dQuote{hours}, \dQuote{days}. Default to \code{"hours"}.
##'
##'  An \code{eventTable} is an object that consists of a data.frame
##'  storing ordered time-stamped events of an (unspecified) PK/PD
##'  dynamic system, units (strings) for dosing and time records, plus a
##'  list of functions to add and extract event records.
##'
##'  Currently, events can be of two types: dosing events that represent
##'  inputs to the system and sampling time events that represent
##'  observations of the system with \sQuote{amount.units} and
##'  \sQuote{time.units}, respectively. In the future, additional events
##'  may include resetting of state variables (compartments), for
##'  instance, to indicate time after \dQuote{wash-out}, etc.
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
##'
##' @keywords models data
##' @concept ordinary differential equations
##' @concept Nonlinear regression
##' @concept Pharmacokinetics (PK)
##' @concept Pharmacodynamics (PD)
##' @export
eventTable <- function(amount.units = NA, time.units = NA){
    .lst <- list()
    if (!missing(amount.units)) .lst$amount.units <- amount.units;
    if (!missing(time.units)) .lst$time.units <- time.units
    .Call(`_RxODE_et_`, .lst, list())
}

##' Sequence of event tables
##'
##' This combines a sequence of event tables.
##'
##' @param ... The event tables and optionally time between event
##'     tables, called waiting times in this help document.
##'
##' @param samples How to handle samples when repeating an event
##'     table.  The options are: \itemize{ \item{"clear"} Clear
##'     sampling records before combining the datasets \item{"use"}
##'     Use the sampling records when combining the datasets }
##'
##' @param waitII This determines how waiting times between events are
##'     handled. The options are:
##'
##' \itemize{
##'
##' \item \code{"smart"} This "smart" handling of waiting times is the
##' default option.  In this case, if the waiting time is above the
##' last observed inter-dose interval in the first combined event
##' table, then the actual time between doses is given by the wait
##' time.  If it is smaller than the last observed inter-dose
##' interval, the time between event tables is given by the inter-dose
##' interval + the waiting time between event tables.
##'
##' \item \code{"+ii"} In this case, the wait time is added to the
##' inter-dose interval no matter the length of the wait time or
##' inter-dose interval
##'
##' }
##' @param ii If there was no inter-dose intervals found in the event
##'     table, assume that the interdose interval is given by this
##'     \code{ii} value.  By default this is \code{24}.
##'
##' @details
##'
##' This \code{seq}uences all the event tables in added in the
##' argument list \code{...}.  By default when combining the event
##' tables the offset is at least by the last inter-dose interval in
##' the prior event table (or \code{ii}).  If you separate any of the
##' event tables by a number, the event tables will be separated at
##' least the wait time defined by that number or the last inter-dose
##' interval.
##'
##' @examples
##'
##' ## Model from RxODE tutorial
##' mod1 <-RxODE({
##'     KA=2.94E-01;
##'     CL=1.86E+01;
##'     V2=4.02E+01;
##'     Q=1.05E+01;
##'     V3=2.97E+02;
##'     Kin=1;
##'     Kout=1;
##'     EC50=200;
##'     C2 = centr/V2;
##'     C3 = peri/V3;
##'     d/dt(depot) =-KA*depot;
##'     d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
##'     d/dt(peri)  =                    Q*C2 - Q*C3;
##'     d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
##' });
##'
##' ## These are making the more complex regimens of the RxODE tutorial
##'
##' ## bid for 5 days
##' bid <- et(timeUnits="hr") %>%
##'        et(amt=10000,ii=12,until=set_units(5, "days"))
##'
##' ## qd for 5 days
##' qd <- et(timeUnits="hr") %>%
##'       et(amt=20000,ii=24,until=set_units(5, "days"))
##'
##' ## bid for 5 days followed by qd for 5 days
##'
##' et <- seq(bid,qd) %>% et(seq(0,11*24,length.out=100));
##'
##' bidQd <- rxSolve(mod1, et)
##'
##' plot(bidQd, C2)
##'
##'
##' ## Now Infusion for 5 days followed by oral for 5 days
##'
##' ##  note you can dose to a named compartment instead of using the compartment number
##' infusion <- et(timeUnits = "hr") %>%
##'       et(amt=10000, rate=5000, ii=24, until=set_units(5, "days"), cmt="centr")
##'
##'
##' qd <- et(timeUnits = "hr") %>% et(amt=10000, ii=24, until=set_units(5, "days"), cmt="depot")
##'
##' et <- seq(infusion,qd)
##'
##' infusionQd <- rxSolve(mod1, et)
##'
##' plot(infusionQd, C2)
##'
##' ## 2wk-on, 1wk-off
##'
##' qd <- et(timeUnits = "hr") %>% et(amt=10000, ii=24, until=set_units(2, "weeks"), cmt="depot")
##'
##' et <- seq(qd, set_units(1,"weeks"), qd) %>%
##'      add.sampling(set_units(seq(0, 5.5,length.out=200),weeks))
##'
##' wkOnOff <- rxSolve(mod1, et)
##'
##' plot(wkOnOff, C2)
##'
##' ## You can also repeat the cycle easily with the rep function
##'
##' qd <-et(timeUnits = "hr") %>% et(amt=10000, ii=24, until=set_units(2, "weeks"), cmt="depot")
##'
##' et <- etRep(qd, times=4, wait=set_units(1,"weeks")) %>%
##'      add.sampling(set_units(seq(0, 12.5,by=0.1),weeks))
##'
##' repCycle4 <- rxSolve(mod1, et)
##'
##' plot(repCycle4, C2)
##'
##' @return A new event table
##'
##' @author Matthew L Fidler
##'
##' @seealso \code{\link{eventTable}}, \code{\link{et}}, \code{\link{etRep}}, \code{\link{etRbind}},
##'    \code{\link{RxODE}}
##'
##' @references
##'
##' Wang W, Hallow K, James D (2015). "A Tutorial on RxODE: Simulating
##' Differential Equation Pharmacometric Models in R." CPT:
##' Pharmacometrics \& Systems Pharmacology, 5(1), 3-10. ISSN 2163-8306,
##' <URL: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728294/>.
##'
##' @export
etSeq <- function(...,samples=c("clear", "use"), waitII=c("smart", "+ii"), ii=24){
    ## etSeq_(List ets, bool clearSampling=clearSampling);
    .sampleIx <- c(clear=0L,use=1L);
    .waitIx <- c(smart=0L, `+ii`=1L)
    .Call(`_RxODE_etSeq_`, list(...), setNames(.sampleIx[match.arg(samples)],NULL),
          setNames(.waitIx[match.arg(waitII)],NULL), as.double(ii), FALSE, 0L,
          0L, TRUE, character(0),logical(0),FALSE);
}
##' Combining event tables
##'
##' @return A new event table
##' @inheritParams etSeq
##' @author Matthew L Fidler
##' @export
etRbind <- function(...,samples=c("use", "clear"),waitII=c("smart", "+ii"),
                    id=c("merge", "unique")){
    .sampleIx <- c(clear=0L,use=1L);
    .waitIx <- c(smart=0L, `+ii`=1L);
    .idIx <- c(merge=0L,unique=1L);
    .Call(`_RxODE_etSeq_`, list(...), setNames(.sampleIx[match.arg(samples)],NULL),
          setNames(.waitIx[match.arg(waitII)],NULL), as.double(0), TRUE,
          setNames(.idIx[match.arg(id)],NULL),
          0L, TRUE, character(0),logical(0),FALSE);
}

##'@export
rbind.rxEt <- function(..., deparse.level = 1){
    if (!missing(deparse.level)) warning("deparse.level not used with RxODE event tables");
    do.call(etRbind,list(...));
}

##'@rdname etSeq
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
##' @return An event table of repeated events
##' @author Matthew Fidler
##' @export
etRep <- function(x, times=1, length.out=NA, each=NA, n=NULL, wait=0, id=integer(0),
                  samples=c("clear", "use"),
                  waitII=c("smart", "+ii"), ii=24){
    if (!is.null(n)){
        times <- n;
    }
    .sampleIx <- c(clear=0L,use=1L);
    .waitIx <- c(smart=0L, `+ii`=1L)
    if (!is.na(length.out)) stop("'length.out' makes no sense with event tables");
    if (!is.na(each)) stop("'each' makes no sense with event tables");
    .Call(`_RxODE_etRep_`, x, as.integer(times),
          wait, as.integer(id), setNames(.sampleIx[match.arg(samples)],NULL),
          setNames(.waitIx[match.arg(waitII)],NULL), as.double(ii))
}

##'@export
rep.rxEt <- function(x, ...){
    do.call(etRep,list(x=x,...));
}
##' Coerce object to data.frame
##'
##' @param x Object to coerce to et.
##' @param ... Other parameters
##'@export
as.et <- function(x,...){
    UseMethod("as.et");
}
##'@rdname as.et
##'@export
as.et.default <- function(x,...){
    .e <- et();
    .e$import.EventTable(as.data.frame(x));
    return(.e);
}

##' Check to see if this is an rxEt object.
##'
##' @param x object to check to see if it is rxEt
##'
##' If this is an rxEt object that has expired strip all rxEt
##' information.
##'
##' @author Matthew L.Fidler
##' @export
is.rxEt <- function(x){
    .Call(`_RxODE_rxIs`, x, "rxEt");
}


##' @importFrom magrittr %>%
##' @export
magrittr::`%>%`
