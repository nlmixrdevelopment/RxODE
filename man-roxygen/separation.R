##' Tells the type of separation strategy when
##' simulating covariance with parameter uncertainty with standard
##' deviations modeled in the \code{thetaMat} matrix.
##'
##' \code{"lkj"} simulates the correlation matrix from the
##' \code{rLKJ1} matrix with the distribution parameter \code{eta}
##' equal to the degrees of freedom \code{nu} by \code{(nu-1)/2}
##'
##' \code{"separation"} simulates from the identity inverse Wishart
##' covariance matrix with \code{nu} degrees of freedom.  This is then
##' converted to a covariance matrix and augmented with the modeled
##' standard deviations.  While computationally more complex than the
##' \code{"lkj"} prior, it performs better when the covariance matrix
##' size is greater or equal to 10
##'
##' \code{"auto"} chooses \code{"lkj"} when the dimension of the
##' matrix is less than 10 and \code{"separation"} when greater
##' than equal to 10.
