##' Simulate random normal variable from threefry generator
##'
##' @inheritParams stats::rnorm
##'
##' @param n number of observations
##'
##' @param ncores Number of cores for the simulation
##'
##' Note, in RxODE models \code{rxnorm()} does not accept \code{n} or
##' \code{ncores}.  These are determined by the \code{rxSolve()}
##' values.
##'
##' Care should be taken with this method not to encounter the
##' birthday problem, described
##' \url{https://www.johndcook.com/blog/2016/01/29/random-number-generator-seed-mistakes/}.
##' Since the \code{sitmo} \code{threefry}, this currently generates
##' one random deviate from the uniform distribution to seed the
##' engine \code{threefry} and then run the code.
##'
##' Therefore, a simple call to the random number generated followed by a second
##' call to random number generated may have identical seeds.  As the number of
##' random number generator calls are increased the probability that the
##' birthday problem will increase.
##'
##' The key to avoid this problem is to run all simulations in the
##' `RxODE` environment once (therefore one seed or series of seeds
##' for the whole simulation) or to pre-generate all random variables
##' used for the simulation.
##'
##' @examples
##'
##' rxnorm(n=10) # with rxnorm you have to explicitly state n
##' rxnorm(n=10,ncores=2) # You can parallelize the simulation using openMP
##'
##' rxnorm(2,3) ## The first 2 arguments are the mean and standard deviation
##'
##'
##' ## This example uses `rxnorm` directly in the model
##'
##' rx <- RxODE({
##'   a = rxnorm()
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxnorm <- function(mean = 0, sd = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(mean, len=1);
    checkmate::assertNumeric(sd, lower=0, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    rxnorm_(mean, sd, n, ncores)
}



rxnormV <- function(mean = 0, sd = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(mean, len=1);
    checkmate::assertNumeric(sd, lower=0, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    rxnorm_(mean, sd, n, ncores)
}
