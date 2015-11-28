\name{eventTable}
\alias{eventTable}
\title{
   Create an event table object
}
\description{
   Initializes an object of class \sQuote{EventTable} with methods for
   adding and querying dosing and observation records
}
\usage{
eventTable(amount.units = NA, time.units = "hours")
}
\arguments{
   \item{amount.units}{
     string denoting the amount dosing units, e.g., \dQuote{mg}, \dQuote{ug}.
     Default to \code{NA} to denote unspecified units.
   }
   \item{time.units}{
     string denoting the time units, e.g., \dQuote{hours}, \dQuote{days}.
     Default to \code{"hours"}.
   }
}
\details{
An \code{EventTable} is an object that consists of a data.frame 
storing ordered time-stamped events of an (unspecified) PK/PD dynamic system, 
units (strings) for dosing and time records, 
plus a list of functions to add and extract event records.

Currently, events can be of two types: dosing events that represent
inputs to the system and sampling time events that 
represent observations of the system with \sQuote{amount.units} and
\sQuote{time.units}, respectively. In the future, additional events
may include resetting of state variables (compartments), for instance,
to indicate time after \dQuote{wash-out}, etc.
}
\value{
A closure with the following list of functions:
   \item{get.EventTable}{returns the current event table.}
   \item{add.dosing}{adds dosing records to the event table.
   Its arguments are

      \code{dose}: numeric scalar, dose amount in \code{amount.units};

      \code{nbr.doses}: integer, number of doses;

      \code{dosing.interval}: required numeric scalar, time between doses
      in \code{time.units}, defaults to 24 of \code{time.units="hours"};

      \code{dosing.to}: integer, compartment the dose goes into 
      (first compartment by default);

      \code{rate}: for infusions, the rate of infusion (default 
      is \code{NULL}, for bolus dosing;
   
      \code{start.time}: required dosing start time;

      \code{do.sampling}: logical, should observation sampling records
      be added at the dosing times? Defaults to \code{FALSE}.

      \code{amount.units}: optional string indicating the dosing units. 
      Defaults to \code{NA} to indicate as per the original \code{EventTable}
      definition.

      \code{time.units}: optional string indicating the time units. 
      Defaults to \code{"hours"} to indicate as per the original \code{EventTable}
      definition.
   }

   \item{get.dosing}{returns a data.frame of dosing records.}
  
   \item{add.sampling}{adds sampling time observation records to the
   event table. Its arguments are

      \code{time} a vector of time values (in \code{time.units}).
  
      \code{time.units} an optional string specifying the time 
      units. Defaults to the units specified when the \code{EventTable}
      was initialized.

      % TODO: should add.sampling() have similar calling sequence
      % as add.dosing()?
      %\code{sampling.interval}: scalar, time between samples.
      %\code{start.time}: scalar, starting observation time.
      %\code{end.time}: scalar, end observation time.
   }

   \item{get.sampling}{returns a data.frame of sampled observation
   records.}

   \item{get.obs.rec}{returns a logical vector indicating 
      whether each event record represents an observation or not.}

   \item{get.nobs}{returns the number of observation (not dosing) records.}

   \item{get.units}{returns a two-element character vector with the
   dosing and time units, respectively.}

   \item{copy}{makes a copy of the current event table. To create
   a copy of an event table object use \code{qd2 <- qd$copy()}.}
   
}
\author{Melissa Hallow and Wenping Wang}

\seealso{
   \code{\link{RxODE}}   
}
\examples{
# create dosing and observation (sampling) events
# QD 50mg dosing, 5 days followed by 25mg 5 days

qd <- eventTable(amount.units = "mg", time.units = "days")

qd$add.dosing(dose=50, nbr.doses=5, dosing.interval = 1, do.sampling=FALSE)

# sample the system's drug amounts hourly the first day, then every 12 hours
# for the next 4 days
qd$add.sampling(seq(from = 0, to = 1, by = 1/24))
qd$add.sampling(seq(from = 1, to = 5, by = 12/24))

#print(qd$get.dosing())     # table of dosing records
print(qd$get.nobs())   # number of observation (not dosing) records

# BID dosing, 5 days

bid <- eventTable("mg", "days")  # only dosing
bid$add.dosing(dose=10000, nbr.doses=2*5, 
   dosing.interval = 12, do.sampling=FALSE)

# Use the copy() method to create a copy (clone) of an existing
# event table (simple assignments just create a new reference to
# the same event table object (closure)).

bid.ext <- bid$copy()      # three-day extension for a 2nd cohort
bid.ext$add.dosing(dose = 5000, nbr.doses = 2*3,
   start.time = 120, dosing.interval = 12, do.sampling = FALSE)
}

\keyword{models}
\keyword{data}
\concept{ordinary differential equations}
\concept{Nonlinear regression}
\concept{Pharmacokinetics (PK)}
\concept{Pharmacodynamics (PD)}
