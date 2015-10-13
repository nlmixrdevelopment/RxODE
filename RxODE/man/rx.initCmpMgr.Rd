\name{rx.initCmpMgr}
\alias{rx.initCmpMgr}
\title{
  A compilation manager for RxODE models
}
\description{
   This function parses, compiles, links, and loads the shared 
   object that implements an RxODE model. 
}
\usage{
   rx.initCmpMgr(model, modName, wd)
}
\arguments{
   \item{model}{
      a string containing the source \code{RxODE} model;
   }
   \item{modName}{
      a string with a model identifier (this string will be used to
      construct filenames and C symbols, therefore it must be suitable
      for creating those identifiers);
   }
   \item{wd}{
      a string with a file path to a working directory where to create 
      various C and object files.
   }
}
\details{
   The function parses and compiles (if needed) the \code{RxODE} model 
   specified in the string \code{model} into a dynamic link library (DLL
   on Windows) or a shared object (\code{*.so} on Unix-like systems).
   It then dynamically loads this code into the current R session. 
   (Models previously parsed and compiled in previous R sessions 
   only need to be dynamically loaded into the current R session.)
}
\value{
   An object (closure) with the following member functions:
   \item{parse}{
      this function parses (translates) the ODE-based model 
      specification and generates a C file to implement the model.}
   \item{compile}{
      compiles the generated C file for the ODE system
      and dynamically loads the machine code in the shared object.}
   \item{dynLoad}{
      if needed, dynamically loads the dynamic library
      produced in the \code{compile()} step.  Note that the shared
      object persists across R sessions, thus the \code{dynLoad} needs
      to be issued as needed.}
   \item{dynUnload}{
      this function unloads the previously dynamically loaded
      model object code.  Mainly for internal use.}
   \item{ode_solver}{
      a string with the name of the C symbol for this model solver.}
   \item{dllfile}{
      a string with the name of the dynamic link (or shared object) file.}
   \item{get.modelVars}{
      function that returns a list with 3 character
      vectors, \code{params}, \code{state}, and \code{lhs} of variable 
      names (identifiers) used in the model specification. 
      These will be output when the model is computed (i.e., the ODE 
      solved).}
   \item{isValid}{
      a function that (naively) checks for model validity,
      namely that the C object code reflects the latest model 
      specification.}
   \item{get.index}{
      helper function to extract the index of one or
      more system variables (state, parameter, or other).}
   \item{getObj}{
      internal (not user callable) function.}
}
\seealso{
  \code{\link{RxODE}}
}
\examples{\dontrun{
   cmpMgt <- rx.initCmpMgr(model, "tst1", wd = ".")
}
}
\keyword{models}
\keyword{ODE}
\concept{ordinary differential equations}
