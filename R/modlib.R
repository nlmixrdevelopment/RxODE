##' RxODE one compartment model (solved)
##' @examples
##' oral1cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot
"oral1cmt"

##' RxODE two compartment model (solved)
##' @examples
##' oral2cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot
"oral2cmt"

##' RxODE three compartment model (solved)
##' @examples
##' oral3cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot
"oral3cmt"

.rxPkgDll <- function(obj){
    obj$mdir <- file.path(system.file(package=obj$package), "include")
    file.path(system.file(package=obj$package), "rx", basename(obj$rxDll$dll))
}

.rxModLib <- function(){
    ## First Create Model
    message("oral1cmt")
    oral1cmt <- RxODE({
        popCl <- 1
        popV <- 20
        popKa <- 1
        bsvCl <- 0
        bsvV  <- 0
        bsvKa <- 0
        cl ~ popCl * exp(bsvCl)
        v ~ popV * exp(bsvV)
        ka ~ popKa * exp(bsvKa)
        popLagDepot <- 0
        popLagCentral <- 0
        popRateCentral <- 0
        popDurCentral <- 0
        bsvLagDepot <- 0
        bsvLagCentral <- 0
        bsvRateCentral <- 0
        bsvDurCentral <- 0
        lag(depot) <- popLagDepot * exp(bsvLagDepot)
        lag(central) <- popLagCentral * exp(bsvLagCentral)
        rate(central) <- popRateCentral *  exp(bsvRateCentral)
        dur(central) <- popDurCentral * exp(bsvDurCentral)
        cp <- linCmt()
    });
    rxDelete(oral1cmt);
    ## Second optimize expressions and recompile using the package= option
    oral1cmt <- RxODE(rxOptExpr(rxNorm(oral1cmt)), package="RxODE", modName="oral1cmt");
    usethis::use_data(oral1cmt, overwrite = TRUE);

    message("oral2cmt")
    oral2cmt <- RxODE({
        popCl <- 1
        popV <- 20
        popKa <- 1
        popVp <- 10
        popQ <- 2
        bsvCl <-0
        bsvV <- 0
        bsvKa <-0
        bsvVp <- 0
        bsvQ <-0
        cl ~ popCl * exp(bsvCl)
        v ~ popV * exp(bsvV)
        ka ~ popKa * exp(bsvKa)
        q ~ popQ * exp(bsvQ)
        vp ~ popVp * exp(bsvVp)
        popLagDepot <- 0
        popLagCentral <- 0
        popRateCentral <- 0
        popDurCentral <- 0
        bsvLagDepot <- 0
        bsvLagCentral <- 0
        bsvRateCentral <- 0
        bsvDurCentral <- 0
        lag(depot) <- popLagDepot * exp(bsvLagDepot)
        lag(central) <- popLagCentral * exp(bsvLagCentral)
        rate(central) <- popRateCentral * exp(bsvRateCentral)
        dur(central) <- popDurCentral * exp(bsvDurCentral)
        cp <- linCmt()
    });
    oral2cmt <- RxODE(rxOptExpr(rxNorm(oral2cmt)), package="RxODE", modName="oral2cmt")
    usethis::use_data(oral2cmt, overwrite = TRUE);

    message("oral3cmt")
    oral3cmt <- RxODE({
        cl <- 1
        v <- 20
        ka <- 1
        vp <- 10
        q <- 2
        q2 <- 2
        vp2 <- 100
        lagDepot <- 0
        lagCentral <- 0
        rateCentral <- 0
        durCentral <- 0
        lag(depot) <- lagDepot
        lag(central) <- lagCentral
        rate(central) <- rateCentral
        dur(central) <- durCentral
        cp <- linCmt()
    });
    oral3cmt <- RxODE(rxOptExpr(rxNorm(oral3cmt)), package="RxODE", modName="oral3cmt")
    usethis::use_data(oral3cmt, overwrite = TRUE);
}
