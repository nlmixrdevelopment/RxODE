% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symengine.R
\name{rxFun}
\alias{rxFun}
\alias{rxRmFun}
\title{Add user function to RxODE}
\usage{
rxFun(name, args, cCode)

rxRmFun(name)
}
\arguments{
\item{name}{This gives the name of the user function}

\item{args}{This gives the arguments of the user function}

\item{cCode}{This is the C-code for the new function}
}
\value{
nothing
}
\description{
This adds a user function to RxODE that can be called.  If needed,
these functions can be differentiated by numerical differences or
by adding the derivatives to RxODE's internal derivative table
with \code{\link[=rxD]{rxD()}}
}
\examples{
\donttest{
## Right now RxODE is not aware of the function f
## Therefore it cannot translate it to symengine or
## Compile a model with it.

try(RxODE("a=fun(a,b,c)"))

## Note for this approach to work, it cannot interfere with C
## function names or reserved RxODE specical terms.  Therefore
## f(x) would not work since f is an alias for bioaviability.

fun <- "
double fun(double a, double b, double c) {
  return a*a+b*a+c;
}
" ## C-code for function

rxFun("fun", c("a", "b", "c"), fun) ## Added function

## Now RxODE knows how to translate this function to symengine

rxToSE("fun(a,b,c)")

## And will take a central difference when calculating derivatives

rxFromSE("Derivative(fun(a,b,c),a)")

## Of course, you could specify the derivative table manually
rxD("fun", list(
  function(a, b, c) {
    paste0("2*", a, "+", b)
  },
  function(a, b, c) {
    return(a)
  },
  function(a, b, c) {
    return("0.0")
  }
))

rxFromSE("Derivative(fun(a,b,c),a)")

# You can also remove the functions by `rxRmFun`

rxRmFun("fun")
}
}
\author{
Matthew L. Fidler
}
