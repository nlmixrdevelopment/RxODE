\name{RxODE}
\alias{RxODE}
\title{
  Create an ODE-based model specification 
}
\description{
  Create a dynamic ODE-based model object suitably for 
  translation into fast C code
}
\usage{
RxODE(model, modName = basename(wd), wd = getwd(), filename = NULL, do.compile = TRUE, ...)
}

\arguments{
   \item{model}{
     a string containing the set of ordinary differential equations (ODE)
     and other expressions defining the changes in the dynamic 
     system (see also the \code{filename} argument).
     For details, see the section \dQuote{\code{RxODE} Syntax} below.
     }
   \item{modName}{
     a string to be used as the model name. This string is
     used for naming various aspects of the computations,
     including generating C symbol names, dynamic libraries, etc.
     Therefore, it is necessary that \code{modName} consists of
     simple ASCII alphanumeric characters starting with a letter.
   }
   \item{wd}{
     character string with a working directory where to create
     a subdirectory according to \code{modName} (defaults to the
     current working directory \code{getwd()}. A subdirectoy 
     named after the \dQuote{\code{modName.d}} will be created and populated
     with a C file, a dynamic loading library, plus various other 
     working files.
   }
   \item{filename}{
     A file name or connection object where the ODE-based model
     specification resides. Only one of \code{model} or \code{filename}
     may be specified.
   }
   \item{do.compile}{
     logical specifying whether to carry out the parsing of the ODE
     into C code and its compilation. Default is \code{TRUE}.
   }
   \item{...}{
   any other arguments are passed to the function \code{\link{readLines}},
   (e.g., encoding).
   }
}
\section{RxODE Syntax}{
   An \code{RxODE} model specification consists of one or more 
   statements terminated by semi-colons, \sQuote{\code{;}}, and
   optional comments (comments are delimited by \code{#} and an 
	end-of-line marker).  \strong{NB:} Comments are not allowed 
	inside statements.

   A block of statements is a set of statements delimeted by
	curly braces, \sQuote{\code{\{ ... \}}}.
   Statements can be either assignments or conditional \code{if}
   statements. Assignment statements can be either \dQuote{simple}
	assignmets, where the left hand is an identifier (i.e., variable), or 
	special \dQuote{time-derivative} assignments, where the left hand
   specifies the change of that variable with respect to time  
	e.g., \code{d/dt(depot)}.

   Expressions in assignment and \sQuote{\code{if}} statements can be 
	numeric or logical (no character expressions are currently supported).
   Numeric expressions can include the following numeric operators 
	(\sQuote{\code{+}}, \sQuote{\code{-}}, \sQuote{\code{*}},
   \sQuote{\code{/}}, \sQuote{\code{^}}),	and
	those mathematical functions defined in the C or the
	R math libraries (e.g., \code{fabs}, \code{exp}, \code{log}, \code{sin}).
   (Notice that the modulo operator \sQuote{\code{\%}} is currently
   not supported.)

   Identifiers in an \code{RxODE} model specification can refer to:
	\itemize{
      \item state variables in the dynamic system (e.g., compartments in a
	   pharmacokinetics/pharamcodynamics model);
		\item implied input variable, \code{t} (time), 
      \code{podo} (oral dose, for absorption models), and
      \code{tlast} (last time point);
		\item model parameters, (\code{ka} rate of absorption, \code{CL} 
      clearance, etc.);
		\item others, as created by assignments as part of the model
		specification.
   }

	The values of these variables at pre-specified time points are
	saved as part of the fitted/integrated/solved model (see
	\code{\link{eventTable}}, in particular its member function 
	\code{add.sampling} that defines a set of time points at which
	to capture a snapshot of the syste via the values of these variables).
}
\value{
An object (closure) of class \dQuote{\code{RxODE}} consisting of
the following list of strings and functions:

   \item{modName}{the name of the model (a copy of the input argument).}
   \item{model}{a character string holding the source model specification.}
   \item{get.modelVars}{a function that returns a list with 3 character 
      vectors, \code{params}, 
      \code{state}, and \code{lhs} of variable names used in the model
      specification. These will be output when the model is computed (i.e.,
      the ODE solved by integration).}
   \item{solve}{this function solves (integrates) the ODE. The following
      are its arguments:

      \code{params}: a numeric named vector of starting estimates for
      every parameter in the ODE system; 
 
      \code{events}: an \code{EventTable} object describing the 
      input (e.g., doses) to the dynamic system and observation sampling time points
      (see  \code{\link{eventTable}});

      \code{inits}: a vector of initial values of the state variables
      (e.g., amounts in each compartment);

      \code{stiff}: a logical (\code{TRUE} by default) indicating whether
      the ODE system is stifff or not;

      \code{trans_abs}: a logical (\code{FALSE} by default) indicating
      whether to fit a transit absorption term 
      (TODO: need further documentation and example);

      \code{atol}: a numeric absolute tolerance (1e-08 by default);

      \code{rtol}: a numeric relative tolerance (1e-06 by default).

      The output of \dQuote{solve} is a matrix with as many rows as there
      are sampled time points and as many columns as system variables 
      (as defined by the ODEs and additional assigments in the RxODE model 
      code).
   }
   \item{isValid}{a function that (naively) checks for model validity,
      namely that the C object code reflects the latest model 
      specification.}
   \item{dynLoad}{a function with one \code{force=FALSE} argument
      that dynamically loads the object code if needed.}
   \item{cmpMgr}{a \dQuote{compilation manager} object, see 
      \code{\link{rx.initCmpMgr}}.}
   \item{delete}{removes all created model files, including C and DDL files.
      The model object is no longer valid and should be removed, e.g., 
      \code{rm(m1)}.}
   \item{run}{deprecated, use \code{solve}.}
   \item{parse}{deprecated.}
   \item{compile}{deprecated.}
   \item{get.index}{helper function to extract the index of one or
      more system variables (state, parameter, or other).
   }
   \item{getObj}{internal (not user callable) function.}
}
\author{Melissa Hallow, Wenping Wang, David A. James}
\seealso{
   \code{\link{eventTable}}
}
\examples{
# Step 1 - Create a model specification
ode <- "
   # A 4-compartment model, 3 PK and a PD (effect) compartment
   # (notice state variable names 'depot', 'centr', 'peri', 'eff')

   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"

m1 <- RxODE(model = ode, modName = "m1")
print(m1)

# Step 2 - Create the model input as an EventTable,
# including dosing and observation (sampling) events

# QD (once daily) dosing for 5 days.

qd <- eventTable(amount.units="ug", time.units = "hours")
qd$add.dosing(dose=10000, nbr.doses=5, dosing.interval = 24)

# Sample the system hourly during the first day, every 8 hours
# then after

qd$add.sampling(0:24)
qd$add.sampling(seq(from = 24+8, to = 5*24, by = 8))

# Step 3 - set starting parameter estimates and initial
# values of the state

theta <- 
   c(KA=.291, CL=18.6, 
     V2=40.2, Q=10.5, V3=297.0,
     Kin=1.0, Kout=1.0, EC50=200.0)

# init state variable
inits <- c(0, 0, 0, 1)      

# Step 4 - Fit the model to the data

qd.cp <- m1$solve(theta, events = qd, inits)

head(qd.cp)

}
\keyword{models}
\keyword{nonlinear}
\concept{Nonlinear regression}
\concept{ODE models}
\concept{Ordinary differential equations}
\concept{Pharmacokinetics (PK)}
\concept{Pharmacodynamics (PD)}
