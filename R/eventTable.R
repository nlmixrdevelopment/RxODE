                                        # event table (dosing + sampling obs from the system)
# An eventTable object contains a numeric matrix with
# a time vector, an event id  describing two types
# of timed records, doses (input) and sampling times
# (state variables); in the future there could be
                                        # other events (e.g., re-setting after "washoout"
# periods, resetting of compartments (e.g., urine),
# etc.
# TODO:
#   (1) Other events (steady state, resetting compartments, etc.)
#   (2) Covariates (age, sex, weight, biomarkers, etc.)
#   (3) A more comprehensive handling of units for time, amounts,
#       and covariates.


#' Create an event table object
#'
#' Initializes an object of class \sQuote{EventTable} with methods for
#' adding and querying dosing and observation records
#'
#' @param amount.units string denoting the amount dosing units, e.g.,
#'     \dQuote{mg}, \dQuote{ug}. Default to \code{NA} to denote
#'     unspecified units.  It could also be a solved RxODE object.  In
#'     that case, eventTable(obj) returns the eventTable that was used
#'     to solve the RxODE object.
#'
#' @param time.units string denoting the time units, e.g.,
#'     \dQuote{hours}, \dQuote{days}. Default to \code{"hours"}.
#'
#' An \code{eventTable} is an object that consists of a data.frame
#' storing ordered time-stamped events of an (unspecified) PK/PD
#' dynamic system, units (strings) for dosing and time records, plus a
#' list of functions to add and extract event records.
#'
#' Currently, events can be of two types: dosing events that represent
#' inputs to the system and sampling time events that represent
#' observations of the system with \sQuote{amount.units} and
#' \sQuote{time.units}, respectively. In the future, additional events
#' may include resetting of state variables (compartments), for
#' instance, to indicate time after \dQuote{wash-out}, etc.
#'
#' @return A closure with the following list of functions:
#'
#' \item{get.EventTable}{returns the current event table.}
#'
#' \item{add.dosing}{adds dosing records to the event table.
#'
#' Its arguments are
#'
#'   \code{dose}: numeric scalar, dose amount in \code{amount.units};
#'
#'   \code{nbr.doses}: integer, number of doses;
#'
#'   \code{dosing.interval}: required numeric scalar, time between doses
#'      in \code{time.units}, defaults to 24 of \code{time.units="hours"};
#'
#'        \code{dosing.to}: integer, compartment the dose goes into
#'        (first compartment by default);
#'
#'        \code{rate}: for infusions, the rate of infusion (default
#'            is \code{NULL}, for bolus dosing;
#'
#'        \code{start.time}: required dosing start time;
#'
#'        \code{do.sampling}: logical, should observation sampling records
#'            be added at the dosing times? Defaults to \code{FALSE}.
#'
#'        \code{amount.units}: optional string indicating the dosing units.
#'           Defaults to \code{NA} to indicate as per the original \code{EventTable}
#'           definition.
#'
#'        \code{time.units}: optional string indicating the time units.
#'           Defaults to \code{"hours"} to indicate as per the original \code{EventTable} definition.
#'     }
#'
#'    \item{get.dosing}{returns a data.frame of dosing records.}
#'
#'    \item{clear.dosing}{clears or deletes all dosing from event table}
#'
#'    \item{add.sampling}{adds sampling time observation records to the
#'        event table. Its arguments are
#'
#'        \code{time} a vector of time values (in \code{time.units}).
#'
#'        \code{time.units} an optional string specifying the time
#'        units. Defaults to the units specified when the \code{EventTable}
#'        was initialized.
#'
#'        % TODO: should add.sampling() have similar calling sequence
#'        % as add.dosing()?
#'        %\code{sampling.interval}: scalar, time between samples.
#'        %\code{start.time}: scalar, starting observation time.
#'        %\code{end.time}: scalar, end observation time.
#'    }
#'
#'    \item{get.sampling}{returns a data.frame of sampled observation
#'        records.}
#'
#'    \item{clear.sampling}{removes all sampling from event table.}
#'
#'    \item{get.obs.rec}{returns a logical vector indicating
#'        whether each event record represents an observation or not.}
#'
#'    \item{get.nobs}{returns the number of observation (not dosing) records.}
#'
#'    \item{get.units}{returns a two-element character vector with the
#'        dosing and time units, respectively.}
#'
#'    \item{copy}{makes a copy of the current event table. To create
#'        a copy of an event table object use \code{qd2 <- qd$copy()}.}
#'
#'    \item{expand}{Expands the event table for multi-subject solving.
#'    This is done by qd$expand(400) for a 400 subject data expansion}
#'
#' @author Melissa Hallow and Wenping Wang
#'
#' @seealso \code{\link{RxODE}}
#'
#' @examples
#' # create dosing and observation (sampling) events
#' # QD 50mg dosing, 5 days followed by 25mg 5 days
#' #
#' qd <- eventTable(amount.units = "mg", time.units = "days")
#' #
#' qd$add.dosing(dose=50, nbr.doses=5, dosing.interval = 1, do.sampling=FALSE)
#' #
#' # sample the system's drug amounts hourly the first day, then every 12 hours
#' # for the next 4 days
#' qd$add.sampling(seq(from = 0, to = 1, by = 1/24))
#' qd$add.sampling(seq(from = 1, to = 5, by = 12/24))
#' #
#' #print(qd$get.dosing())     # table of dosing records
#' print(qd$get.nobs())   # number of observation (not dosing) records
#' #
#' # BID dosing, 5 days
#' bid <- eventTable("mg", "days")  # only dosing
#' bid$add.dosing(dose=10000, nbr.doses=2*5,
#'                dosing.interval = 12, do.sampling=FALSE)
#' #
#' # Use the copy() method to create a copy (clone) of an existing
#' # event table (simple assignments just create a new reference to
#' # the same event table object (closure)).
#' #
#' bid.ext <- bid$copy()      # three-day extension for a 2nd cohort
#' bid.ext$add.dosing(dose = 5000, nbr.doses = 2*3,
#'                    start.time = 120, dosing.interval = 12, do.sampling = FALSE)
#'
#' # You can also use the Piping operator to create a table
#'
#' qd2 <- eventTable(amount.units="mg", time.units="days") %>%
#'     add.dosing(dose=50, nbr.doses=5, dosing.interval=1, do.sampling=FALSE) %>%
#'     add.sampling(seq(from=0, to=1, by=1 / 24)) %>%
#'     add.sampling(seq(from=1, to=5, by=12 / 24))
#' #print(qd2$get.dosing())     # table of dosing records
#' print(qd2$get.nobs())   # number of observation (not dosing) records
#'
#' # Note that piping with %>% will update the original table.
#'
#' qd3 <- qd2 %>% add.sampling(seq(from=5, to=10, by=6 / 24))
#' print(qd2$get.nobs())
#' print(qd3$get.nobs())
#' @keywords models data
#' @concept ordinary differential equations
#' @concept Nonlinear regression
#' @concept Pharmacokinetics (PK)
#' @concept Pharmacodynamics (PD)
#' @export
eventTable <- function(amount.units = NA, time.units = "hours")
{
    .EventTable <- NULL
    .obs.rec <- logical(0)     # flag for observation records

    .amount.units <- as.vector(amount.units)  # preferred units
    .time.units <- as.vector(time.units)

    "add.dosing" <-
        function(dose,      # amount per dose,
                 nbr.doses = 1,      # single dose default
                 dosing.interval = 24,
                 dosing.to=1,         #to which cmt dosing is admin'ed
                 rate=NULL,            #infusion rate if infusion
                 amount.units = NA,
                 start.time,
                 do.sampling=FALSE,
                 time.units = NA, ...)
        {
            if(nbr.doses < 1L){
                stop("Number of Doses must be at least one.");
            }
            tmp <- as.integer(nbr.doses)
            if(nbr.doses != tmp){
                stop("Number of doses must be a integer >= 1.")
            }
            if(!is.na(amount.units)){
                if(is.na(.amount.units))
                    .amount.units <<- as.vector(amount.units)   # initialize
                else if(tolower(.amount.units)!=tolower(amount.units)){
                    stop("dosing units differ from EventTable's")
                }
            } # else assume amount.units as per eventTable() definition

            if(!is.na(time.units)){
                if(is.na(.time.units))
                    .time.units <<- as.vector(time.units)   # initialize
                else if(tolower(.time.units)!=tolower(time.units)){
                    stop("time units differ from EventTable's")
                }
            } # else assume time.units as per eventTable() definition

            if(missing(dosing.interval) && nbr.doses>1 )
                stop("must specify 'dosing.interval' with multiple doses")

            if(missing(start.time)){
                if(is.null(.EventTable) || all(.obs.rec))
                    start.time <- 0
                else {
                    warning("imputing start.time", immediate = TRUE)
                    start.time <- max(.EventTable$time) + dosing.interval
                }
            }
            time <- start.time+(1:nbr.doses-1)*dosing.interval

            ## TODO: should we code individual flags (infusion vs bolus, etc)
            ## in the table and convert to a mask integer just prior to
            ## invoking the C code?
            ## TODO: Handle units. Check that add.dosing() units don't conflict
            ## with the eventTable definition (preferred units)
            if (is.null(rate)) {#-- bolus
                wh <- floor(dosing.to / 100) * 1e5 + 100*(dosing.to %% 100) + 1
                inp <- data.frame(time=time, evid=wh, amt=dose)
            } else {         #-- infusion
                wh <- floor(dosing.to / 100) * 1e5 + 1e4 +100 * (dosing.to %% 100) + 1
                toff <- dose/rate
                if (rate<=0) {
                    inp <- NULL
                } else {
                    inp <- rbind(
                        data.frame(time=time,      evid=wh, amt=rate),
                        data.frame(time=time+toff, evid=wh, amt=-rate)
                    )
                }
            }

            inp <- rbind(.EventTable, inp)
            inp <- inp[order(inp$time, -inp$evid), ]
            .EventTable <<- inp
            .obs.rec <<- inp$evid==0

            s <- as.list(match.call(expand.dots = TRUE))
            if ("sampling.interval" %in% names(s))
                sampling.interval <- s$sampling.interval
            else sampling.interval <- 1

            if (do.sampling)
                add.sampling(0:(nbr.doses*dosing.interval), time.units = as.vector(time.units))
            invisible()
        }

    "add.sampling" <-
        function(time, time.units = NA)
        {
            if(!is.na(time.units)){
                if(is.na(.time.units))
                    .time.units <<- as.vector(time.units)   # initialize
                else if(tolower(.time.units)!=tolower(time.units)){
                    stop("time units differ from EventTable's")
                }
            } # else assume time.units as per eventTable() definition
            inp <- data.frame(time=time, evid=0, amt=NA)
            inp <- rbind(.EventTable, inp)
            inp <- inp[order(inp$time, -inp$evid), ]
            .EventTable <<- inp
            .obs.rec <<- inp$evid==0

            invisible()
        }

    "clear.sampling" <- function(){
        ## Clears all sampling.
        .EventTable <<- .EventTable[!.obs.rec, ,drop = TRUE] ;
        if (is(.EventTable,"list")){
            .EventTable <<- as.data.frame(.EventTable);
        }
        .obs.rec <<- .EventTable$evid == 0
        invisible()
    }

    "clear.dosing" <- function(){
        .EventTable <<- .EventTable[.obs.rec, ,drop = TRUE] ;
        if (is(.EventTable,"list")){
            .EventTable <<- as.data.frame(.EventTable);
        }
        .obs.rec <<- .EventTable$evid==0
        invisible()
    }

    "import.EventTable" <-
        function(inp)
        {
            if (!is.data.frame(inp))
                stop("input table is not a data.frame")
            vars <- setdiff(c("time", "evid", "amt"), names(inp))
            if (length(vars)) {
                msg <- paste("var(s) not found in input table.\n",
                             paste(vars, collapse=" "))
                stop(msg)
            }
            inp <- inp[, c("time", "evid", "amt")]
            inp <- inp[order(inp$time, -inp$evid), ]
            .EventTable <<- inp        # should we append inp to current EventTable?
            .obs.rec <<- inp$evid==0
            invisible()
        }

    "copy" <-
        function()
        {
                                        # Make a copy (clone) of the current event table.
                                        # Can test the output with identical(old, new, ignore.closure=TRUE).
            self <- environment(add.dosing)     # current environment
            nms <- objects(all.names = TRUE, envir = self)

            out <- eventTable()                 # new, pristine eventTable
            env2 <- environment(out$add.dosing)
            for(name in nms)
                assign(name, get(name, envir = self), envir = env2)
            out
        }

    out <-
        list(
            get.EventTable = function() .EventTable,
            get.obs.rec = function() .obs.rec,
            get.nobs = function() sum(.obs.rec),
            add.dosing = add.dosing,
            clear.dosing = clear.dosing,
            get.dosing = function() .EventTable[!.obs.rec, ,drop = FALSE],
            add.sampling = add.sampling,
            clear.sampling = clear.sampling,
            get.sampling = function() .EventTable[.obs.rec, ,drop = FALSE],
            get.units = function() c(dosing = .amount.units, time = .time.units),
            import.EventTable = import.EventTable,
            copy = copy
            )
    class(out) <- "EventTable"
    out
}

#' @export
"print.EventTable" <-
function(x, ...)
{
   nobs <- x$get.nobs()
   dr <- sum(!x$get.obs.rec())
   nr <- nrow(x$get.EventTable())
   nr <- if(is.null(nr)) 0 else nr
   unts <- x$get.units()
   cat(
      sprintf("EventTable with %d records:\n", nr),
      sprintf("  %d dosing (in %s) records\n", dr, unts[1]),
      sprintf("  %d observation time (in %s) records\n", nobs, unts[2])
   )
   invisible(x)
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
    .Call(`_RxODE_add_dosing_`, eventTable, dose, nbr.doses, dosing.interval, dosing.to, rate, amount.units, start.time, do.sampling, time.units)
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
##' @param ... Other parameters (ignored)
##' @return eventTable with updated sampling.  (Note the event table
##'     will be updated even if you don't reassign the eventTable)
##' @author Matthew L. Fidler
##' @seealso \code{\link{eventTable}}, \code{\link{RxODE}}
##' @export
add.sampling <- function(eventTable, time, time.units = NA, ...){
    .Call(`_RxODE_add_sampling_`, eventTable, time, time.units)
}


##' @importFrom magrittr %>%
##' @export
magrittr::`%>%`


##' @export
print.RxODE.multi.data <- function(x, ...){
    message("RxODE multi-subject data:")
    message(sprintf("  Number of Subjects: %s", x$nSub))
    message(sprintf("  Number of Observations: %s (t=%s to %s%s)", x$nObs, x$min.time, x$max.time, ifelse(is.na(x$time.units), "", paste0(" ", x$time.units))))
    if (x$nDose == 0){
        message("  No Dosing Records.")
    } else {
        message(sprintf("  Number of Dosing Records: %s%s", x$nDose, ifelse(is.na(x$amount.units), "", paste0(" (in ", x$amount.units, ")"))))
    }
    obs.cov <- x$cov.names;
    sim.cov <- x$simulated.vars;
    if (!is.null(obs.cov)){
        message(sprintf("  Covariates: %s", paste(obs.cov, collapse=", ")))
    }
    if (!is.null(sim.cov)){
        message(sprintf("  Simulated Variables: %s", paste(sim.cov, collapse=", ")))
    }
}
