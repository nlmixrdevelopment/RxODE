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
##' rxpois(lambda=3, n=10) # with rxpois you have to explicitly state n
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
##' rxt(df=3, n=10) # with rxt you have to explicitly state n
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

##' Simulate uniform variable from threefry generator
##'
##' @inheritParams stats::runif
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxunif(min=0, max=4, n=10) # with rxunif you have to explicitly state n
##' rxunif(min=0, man=4, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxunif()
##'
##'
##' ## This example uses `rxunif` directly in the model
##'
##' rx <- RxODE({
##'   a = rxunif(0,3)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxunif <- function(min = 0, max = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(min, len=1);
    checkmate::assertNumeric(max, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxunif_`, min, max, n, ncores)
}


##' Simulate weibull variable from threefry generator
##'
##' @inheritParams stats::rweibull
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxweibull(shape=1, scale=4, n=10) # with rxweibull you have to explicitly state n
##' rxweibull(shape=1, scale=4, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxweibull(3)
##'
##'
##' ## This example uses `rxweibull` directly in the model
##'
##' rx <- RxODE({
##'   a = rxweibull(1,3)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxweibull <- function(shape, scale = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(shape, len=1);
    checkmate::assertNumeric(scale, len=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxweibull_`, shape, scale, n, ncores)
}


##' Simulate geometric variable from threefry generator
##'
##' @inheritParams stats::rgeom
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxgeom(0.5, n=10) # with rxgeom you have to explicitly state n
##' rxgeom(0.25, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxgeom(0.75)
##'
##'
##' ## This example uses `rxgeom` directly in the model
##'
##' rx <- RxODE({
##'   a = rxgeom(0.24)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxgeom <- function(prob, n=1L, ncores=1L){
    checkmate::assertNumeric(prob, len=1, lower=0, upper=1);
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxgeom_`, prob, n, ncores)
}


##' Simulate beta variable from threefry generator
##'
##' @inheritParams stats::rbeta
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxbeta(0.5, 0.5, n=10) # with rxbeta you have to explicitly state n
##' rxbeta(5, 1, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxbeta(1, 3)
##'
##'
##' ## This example uses `rxbeta` directly in the model
##'
##' rx <- RxODE({
##'   a = rxbeta(2, 2)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxbeta <- function(shape1, shape2, n=1L, ncores=1L){
    checkmate::assertNumeric(shape1, len=1, lower=0);
    if (shape1 == 0) stop("'shape1' cannot be 0");
    checkmate::assertNumeric(shape2, len=1, lower=0);
    if (shape2 == 0) stop("'shape2' cannot be 0");
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxbeta_`, shape1, shape2, n, ncores)
}

##' Simulate gamma variable from threefry generator
##'
##' @inheritParams stats::rgamma
##' @inheritParams rxnormV
##'
##' @template birthdayProblem
##' @examples
##'
##' ## Use threefry engine
##'
##' rxbeta(0.5, 0.5, n=10) # with rxbeta you have to explicitly state n
##' rxbeta(5, 1, n=10, ncores=2) # You can parallelize the simulation using openMP
##'
##' rxbeta(1, 3)
##'
##'
##' ## This example uses `rxbeta` directly in the model
##'
##' rx <- RxODE({
##'   a = rxbeta(2, 2)
##' })
##'
##' et <- et(1,id=1:2)
##'
##' s <- rxSolve(rx,et)
##'
##' @export
rxgamma <- function(shape, rate = 1/scale, scale = 1, n=1L, ncores=1L){
    checkmate::assertNumeric(shape, len=1, lower=0);
    if (shape == 0) stop("'shape' cannot be 0");
    checkmate::assertNumeric(rate, len=1, lower=0);
    if (rate == 0 || scale == 0) stop("'rate'/'scale' cannot be 0");
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15)
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
    checkmate::assertCount(n)
    checkmate::assertCount(ncores)
    rxSeedEng(ncores)
    .Call(`_RxODE_rxgamma_`, shape, rate, n, ncores)
}
