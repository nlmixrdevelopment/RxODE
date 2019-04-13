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
    obj$mdir <- file.path(system.file(package=obj$package), "rx")
    file.path(system.file(package=obj$package), "rx", basename(obj$rxDll$dll))
}

##' Use model object in your package
##' @param obj model to save.
##' @export
rxUse <- function(obj, internal = FALSE, overwrite = TRUE, compress = "bzip2"){
    rxReq("usethis")
    rxReq("devtools")
    if (missing(obj)){
        .env <- new.env();
        assign("internal", internal, .env)
        assign("overwrite", overwrite, .env)
        assign("compress", compress, .env)
        for (.f in list.files(path=devtools::package_file("data"),pattern="\\.rda$",full.names=TRUE)){
            load(.f, envir=.env)
            .f2 <- basename(.f);
            .f2 <- substr(.f2, 0, nchar(.f2) - 4)
            if (rxIs(.env[[.f2]], "RxODE")){
                message(sprintf("Recompile %s (if needed)", .f2))
                eval(parse(text=sprintf("RxODE::rxUse(%s, internal=internal, overwrite=overwrite, compress=compress)", .f2)),
                     envir=.env)
            }
        }
    } else {
        .modName <- as.character(substitute(obj));
        .pkg <- basename(usethis::proj_get())
        .env <- new.env();
        assign(.modName, RxODE(rxNorm(obj), package=.pkg, modName=.modName), .env);
        assign("internal", internal, .env)
        assign("overwrite", overwrite, .env)
        assign("compress", compress, .env)
        eval(parse(text=sprintf("usethis::use_data(%s, internal=internal, overwrite=overwrite, compress=compress)", .modName)),
             envir=.env)

    }
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
        popCl <- 1
        popV <- 20
        popKa <- 1
        popVp <- 10
        popQ <- 2
        popQ2 <- 2
        popVp2 <- 100
        bsvCl <- 0
        bsvV <- 0
        bsvKa <- 0
        bsvVp <- 0
        bsvQ <- 0
        bsvQ2 <- 0
        bsvVp2 <- 0
        cl ~ popCl * exp(bsvCl)
        v ~ popV * exp(bsvV)
        ka ~ popKa * exp(bsvKa)
        q ~ popQ * exp(bsvQ)
        vp ~ popVp * exp(bsvVp)
        q2 ~ popQ2 * exp(bsvQ2)
        vp2 ~ popVp2 * exp(bsvVp2)
        popLagDepot <- 0
        popLagCentral <- 0
        popRateCentral <- 0
        popDurCentral <- 0
        bsvLagDepot <- 0
        bsvLagCentral <- 0
        bsvRateCentral <- 0
        bsvDurCentral <- 0
        lag(depot)    <- popLagDepot * exp(bsvLagDepot)
        lag(central)  <- popLagCentral * exp(bsvLagCentral)
        rate(central) <- popRateCentral * exp(bsvRateCentral)
        dur(central)  <- popDurCentral * exp(bsvDurCentral)
        cp <- linCmt()
    });
    oral3cmt <- RxODE(rxOptExpr(rxNorm(oral3cmt)), package="RxODE", modName="oral3cmt")
    usethis::use_data(oral3cmt, overwrite = TRUE);
}
