##' @export
##' @rdname rxnormV
rxnorm <- function(mean = 0, sd = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(mean, len=1);
    checkmate::assertNumeric(sd, lower=0, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxnorm_`, mean, sd, n, ncores)
}

##' Simulate random normal variable from threefry/vandercorput generator
##'
##' @inheritParams stats::rnorm
##'
##' @param n number of observations
##'
##' @param ncores Number of cores for the simulation
##'
##' `rxnorm` simulates using the threefry sitmo generator; `rxnormV`
##' uses the vandercorput generator
##'
##' @examples
##'
##' ## Use threefry engine
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
##' ## Use vandercorput generator
##'
##' rxnormV(n=10) # with rxnorm you have to explicitly state n
##' rxnormV(n=10,ncores=2) # You can parallelize the simulation using openMP
##'
##' rxnormV(2,3) ## The first 2 arguments are the mean and standard deviation
##'
##'
##' ## This example uses `rxnormV` directly in the model
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

##' Simulate random poisson variable from threefry generator
##'
##' @inheritParams stats::rpois
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxpois(lambda=3, n=10) # with rxnorm you have to explicitly state n
##' rxpois(lambda=3, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxpois(4) ## The first arguments are the lambda parameter
##'
##'
##' ## This example uses `rxpois` directly in the model
##'
##' rx <- RxODE({
##'   a = rxpois(3)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxpois <- function(lambda, n=1L, ncores=1L){
    checkmate::assertNumeric(lambda, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxpois_`, lambda, n, ncores)
}


##' Simulate student t variable from threefry generator
##'
##' @inheritParams stats::rt
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxt(df=3, n=10) # with rxnorm you have to explicitly state n
##' rxt(df=3, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxt(4) ## The first argument is the df parameter
##'
##'
##' ## This example uses `rxt` directly in the model
##'
##' rx <- RxODE({
##'   a = rxt(3)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxt <- function(df, n=1L, ncores=1L){
    checkmate::assertNumeric(df, len=1, lower=0);
    if (df == 0) stop("'df' must be greater than 0")
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxt__`, df, n, ncores)
}
