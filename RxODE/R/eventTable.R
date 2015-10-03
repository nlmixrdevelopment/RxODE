# event table (dosing + sampling obs from the system)
# An EventTable object contains a numeric matrix with
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

"eventTable" <-
function(amount.units = NA, time.units = "hours")
{
   .EventTable <- NULL
   .obs.rec <- logical(0)     # flag for observation records

   .amount.units <- amount.units  # preferred units
   .time.units <- time.units

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
      if(!is.na(amount.units)){
         if(is.na(.amount.units)) 
            .amount.units <<- amount.units   # initialize
         else if(tolower(.amount.units)!=tolower(amount.units)){
            stop("dosing units differ from EventTable's")
        }
      } # else assume amount.units as per eventTable() definition

      if(!is.na(time.units)){
         if(is.na(.time.units)) 
            .time.units <<- time.units   # initialize
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
      
      # TODO: should we code individual flags (infusion vs bolus, etc)
      # in the table and convert to a mask integer just prior to
      # invoking the C code?
      # TODO: Handle units. Check that add.dosing() units don't conflict 
      # with the eventTable definition (preferred units)
      if (is.null(rate)) {#-- bolus
         wh <- 100*dosing.to+1
         inp <- data.frame(time=time, evid=wh, amt=dose)
      } else {         #-- infusion
         wh <- 10000+100*dosing.to+1
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
         add.sampling(0:(nbr.doses*dosing.interval), time.units = time.units)
      invisible()
   }

   "add.sampling" <- 
   function(time, time.units = NA)
   {
      if(!is.na(time.units)){
         if(is.na(.time.units)) 
            .time.units <<- time.units   # initialize
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
         get.dosing = function() .EventTable[!.obs.rec, ,drop = FALSE],
         add.sampling = add.sampling,
         get.sampling = function() .EventTable[.obs.rec, ,drop = FALSE],
         get.units = function() c(dosing = .amount.units, time = .time.units),
         import.EventTable = import.EventTable,
         copy = copy
      )
   class(out) <- "EventTable"
   out
}

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
