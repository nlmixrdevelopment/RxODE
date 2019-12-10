##' Simulate random normal variable from threefry generator
##'
##' @inheritParams stats::rnorm
##'
##' @param n number of observations
##'
##' @param ncores Number of cores for the simulation
##'
##' @template birthdayProblem
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
    .Call(`_RxODE_rxnorm_`, mean, sd, n, ncores)
}



##' Simulate random normal variable from vandercorput generator
##'
##' @inheritParams stats::rnorm
##'
##' @param n number of observations
##'
##' @param ncores Number of cores for the simulation
##'
##' @template birthdayProblem
##'
##' @examples
##'
##' rxnormV(n=10) # with rxnorm you have to explicitly state n
##' rxnormV(n=10,ncores=2) # You can parallelize the simulation using openMP
##'
##' rxnormV(2,3) ## The first 2 arguments are the mean and standard deviation
##'
##'
##' ## This example uses `rxnorm` directly in the model
##'
##' rx <- RxODE({
##'   a = rxnormV()
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxnormV <- function(mean = 0, sd = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(mean, len=1);
    checkmate::assertNumeric(sd, lower=0, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxnormV_`, mean, sd, n, ncores)
}
