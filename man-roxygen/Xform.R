##' \code{thetaMat} simulations (using the separation strategy for
##' covariance simulation), how should the `thetaMat` values be
##' turned int standard deviation values:
##'
##' \itemize{
##' \item \code{identity} This is when standard deviation values are
##' directly modeled by the \code{params} and \code{thetaMat} matrix
##'
##' \item \code{variance} This is when the \code{params} and \code{thetaMat}
##' simulates the variance that are directly modeled by the
##' \code{thetaMat} matrix
##'
##' \item \code{log} This is when the \code{params} and \code{thetaMat}
##' simulates \code{log(sd)}
##'
##' \item \code{nlmixrSqrt} This is when the \code{params} and
##' \code{thetaMat} simulates the inverse cholesky decomposed matrix
##' with the \code{x^2} modeled along the diagonal.  This only works
##' with a diagonal matrix.
##'
##' \item \code{nlmixrLog} This is when the \code{params} and
##' \code{thetaMat} simulates the inverse cholesky decomposed matrix
##' with the \code{exp(x^2)} along the diagonal.  This only works
##' with a diagonal matrix.
##'
##' \item \code{nlmixrIdentity} This is when the \code{params} and
##' \code{thetaMat} simulates the inverse cholesky decomposed matrrix.
##' This only works with a diagonal matrix.
##'
