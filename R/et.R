##' Event Table Function
##'
##' @param ... Times or event tables.
##'
##' @param time Time is the time of the dose or the sampling times.
##'     This can also be unspecified and is determined by the object
##'     type (list or numeric/integer).
##'
##' @param amt Amount of the dose. If specified, this assumes a dosing
##'     record, instead of a sampling record.
##'
##' @param evid Event ID; This can be:
##'
##' \itemize{
##'
##' \item{0} An observation. This can also be specified as
##'    \code{evid=obs}
##'
##' \item{1} A dose observation.  This can also be specified as
##'    \code{evid=dose}
##'
##' \item{2} A non-dose event. This can also be specified as
##'    \code{evid=other}.
##'
##' \item{3} A reset event.  A reset event resets all the compartment
##' values to zero and turns off all infusions.  This can also be
##' specified as \code{evid=reset}.
##'
##' \item{4} Dose and reset event.  This can also be specified as
##' \code{evid=doseReset} or \code{evid=resetDose}
##'
##' }
##'
##' @param cmt Compartment name or number.  If a number, this is an
##'     integer starting at 1.  Negative compartments are not
##'     supported (there is no way to turn off a compartment
##'     currently). If the compartment is a name, the compartment name
##'     is changed to the correct state/compartment number before
##'     running the simulation.
##'
##'     Can also specify \code{cmt} as \code{dosing.to},
##'     \code{dose.to}, \code{doseTo}, \code{dosingTo}, and
##'     \code{state}.
##'
##' @param ii When specifying a dose, this is the inter-dose interval
##'     for \code{ss}, \code{addl} and \code{until} options (described below).
##'
##' @param addl The number of additional doses at a inter-dose
##'     interval after one dose.
##'
##' @param ss Steady state flag;  It can be one of:
##' \itemize{
##'
##' \item{0} This dose is not a steady state dose
##'
##' \item{1} This dose is a steady state dose with the between/inter
##' dose interval of \code{ii}
##'
##' \item{2} This is a steady state dose that uses the super-position
##' principle to allow more complex steady states, like 10 mg in the
##' morning and 20 mg at night, or dosing at 8 am 12 pm and 8 pm
##' instead of every 12 hours.  Since it uses the super positioning
##' principle, it only makes sense when you know the kinetics are
##' linear.
##' }
##'
##' All other values of \code{SS} are currently invalid.
##'
##' @param rate When positive, this is the rate of infusion.  Otherwise:
##'
##' \itemize{
##'
##' \item{0} No infusion is on this record.
##'
##' \item{-1} Rate of this record is modeled by \code{rate(cmt) =} in
##' the RxODE model.  You may also specify type or rate by
##' \code{rate=model}
##'
##' \item{-2} Duration of this record is modeled by \code{dur(cmt) =}
##' in the RxODE model. You may also specify this type of rate by
##' \code{dur=model} or \code{rate=dur}.
##'
##' }
##'
##' When a modeled bioavailability is applied to positive rates
##' (\code{rate} > 0), the duration of infusion is changed. This is
##' because the data specify the rate and amount, the only think ghat
##' modeled bioavailability can affect is duration.
##'
##' If instead you want the modeled bioavailability to increase the
##' rate of infusion instead of the duration of infusion, specify the
##' \code{dur} instead or model the duration with \code{rate=2}.
##'
##' @param dur Duration of infusion.  When \code{amt} and \code{dur}
##'     are specified the rate is calculated from the two data items.
##'     When \code{dur} is specified instead of \code{rate}, the
##'     bioavailability changes will increase rate instead of
##'     duration.
##'
##' @param until This is the time until the dosing should end.  It can
##'     be an easier way to figure out how many additional doses are
##'     needed over your sampling period.
##'
##' @param id
##'
##' @param amountUnits The units for the dosing records (\code{amt})
##'
##' @param timeUnits The units for the time records (\code{time})
##'
##' @param addSampling This is a boolean indicating if a sampling time
##'     should be added at the same time as a dosing time.  By default
##'     this is \code{FALSE}.
##'
##' @template etExamples
##'
##' @export
et <- function(...){
    .lst <- as.list(match.call()[-1]);
    if (rxIs(.lst[1], "numeric") || rxIs(.lst[1], "integer") ||
        rxIs(.lst[1], "list") || rxIs(.lst[1], "rxEt")){
        ## Use do call on a match.call() to preserve lazy evaluation
        ## By doing this evid=obs will work as well as evid="obs" and evid=1
        do.call(et.default, .lst)
    } else {
        UseMethod("et");
    }
}

##'@rdname
##'@export
et.default <- function(...,time, amt, evid, cmt, addl, ss, rate, dur, until, id,
                       amountUnits, timeUnits, addSampling){
    .lst <- as.list(match.call()[-1]);
    if (!missing(time)){
        .lst$time <- time;
    }
    if (!missing(amt)){
        .lst$amt <- amt;
    }
    if (!missing(evid)){
        .evid <- as.character(substitute(evid))
        if (.evid=="obs" || .evid=="0"){
            .lst$evid <- 0L;
        } else if (.evid=="dose" || .evid=="1") {
            .lst$evid <- 1L;
        } else if (.evid=="other" || .evid=="2") {
            .lst$evid <- 2L;
        } else if (.evid=="reset" || .evid=="3") {
            .lst$evid <- 3L;
        } else if (.evid=="doseReset" || .evid=="resetDose" || .evid=="4") {
            .lst$evid <- 4L;
        } else {
            .lst$evid <- as.integer(evid);
        }
    }
    if (!missing(cmt)){
        .cmt <- as.character(substitute(cmt));
        .cmt1 <- try(suppressWarnings(as.integer(cmt)), silent=TRUE);
        if (inherits(.cmt1, "try-error")){
            .lst$cmt <- .cmt
        } else {
            if (is.na(.cmt1)){
                .lst$cmt <- .cmt
            } else {
                .lst$cmt <- .cmt1
            }
        }
    }
    if (!missing(rate)){
        .rate <- as.character(substitute(rate));
        if (.rate=="model" || .rate=="modeled" ||
            .rate=="modelled" || .rate=="rate"){
            .lst$rate <- -1.0
        } else if (.rate=="dur" || .rate=="duration"){
            .lst$rate <- -2.0;
        } else {
            .lst$rate <- rate;
        }
    }
    if (!missing(dur)){
        .dur <- as.character(substitute(dur));
        if (.dur=="model" || .dur=="modeled" ||
            .dur=="modelled" || .dur=="dur" ||
            .dur=="duration"){
            .lst$rate <- -2.0
            .lst <- .lst[names(.lst) != "dur"]
        } else if (.dur=="rate"){
            .lst$rate <- -1.0;
            .lst <- .lst[names(.lst) != "dur"];
        } else {
            .lst$dur <- dur;
        }
    }
    .unitNames <- names(.lst);
    .unitNames <- .unitNames[regexpr("^(amount|time)",.unitNames) != -1];
    .unitNames <- .unitNames[.unitNames != "time"];
    for (.u in .unitNames){
        if (class(.lst[[.u]])=="name"){
            .tmp <- .lst[[.u]];
            .tmp <- deparse(substitute(.tmp))
            .lst[[.u]] <- .tmp
        }
    }
    .lst <- lapply(.lst,function(x){
        eval(x,parent.frame(4L))
    });
    .Call(`_RxODE_et_`, .lst, list())
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
##'     in \code{time.units}, defaults to 24 of
##'     \code{time.units="hours"};
##' @param dosing.to integer, compartment the dose goes into (first
##'     compartment by default);
##' @param rate for infusions, the rate of infusion (default is
##'     \code{NULL}, for bolus dosing;
##' @param amount.units optional string indicating the dosing units.
##'     Defaults to \code{NA} to indicate as per the original
##'     \code{EventTable} definition.
##' @param start.time required dosing start time;
##' @param do.sampling logical, should observation sampling records be
##'     added at the dosing times? Defaults to \code{FALSE}.
##' @param time.units optional string indicating the time units.
##'     Defaults to \code{"hours"} to indicate as per the original
##'     \code{EventTable} definition.
##' @param ... Other parameters passed to \code{\link{et}}.
##' @return eventTable with updated dosing (note the event table will
##'     be updated anyway)
##' @author Matthew L. Fidler
##' @seealso \code{\link{eventTable}}, \code{\link{RxODE}},
##'     \code{\link{et}}, \code{\link{add.sampling}},
##'     \code{\link{etRep}}, \code{\link{etSeq}},
##'     \code{\link{etRbind}}
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
##' @param time.units an optional string specifying the time
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
##' @template etExamples
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
##' @examples
##'
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

##'@rdname etRbind
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
