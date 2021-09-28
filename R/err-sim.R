# This handles the errors for simulations

##' This is a S3 method for getting the distribution lines for a simulatin
##'
##' @param env Parsed RxODE model environment
##' @param pred1 Predictions 1
##' @return Lines for the simulation of `ipred` and `dv`. This is based on the idea that the focei parameters are defined
##' @author Matthew Fidler
##' @noRd
rxGetDistributionSimulationLines <- function(env, pred1) {
  UseMethod("rxGetDistributionSimulationLines")
}

.simulationFun <- list(
  "t"="rt"
)

.getQuotedDistributionAndSimulationArgs <- function(env, pred1) {
  .dist <- pred1$distribution
  .nargs <- max(.errDist[[.dist]])
  .cnd <- pred1$cond
  .argName <- .namedArgumentsToPredDf[[.dist]]
  .args <- vapply(seq(1:.nargs), function(.i){
    .curDist <- .argName[.i]
    if (!is.na(pred1[[.curDist]])) {
      return(pred1[[.curDist]])
    } else {
      .curDist <- paste0(.dist, ifelse(.i == 1, "", .i))
      .w <- which(env$iniDf$err == .curDist & env$iniDf$condition == .cnd)
      if (length(.w) == 1) {
        return(env$iniDf$name[.w])
      } else {
        return("")
      }
    }
  }, character(1))
  as.call(lapply(c(.simulationFun[[.dist]], .args[.args != ""]), .enQuote))
}

rxGetDistributionSimulationLines.norm <- function(env, pred1) {
  .ret <- vector("list", 2)
  .ret[[1]] <- bquote(ipred=rxTBSi(rx_pred_, rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  .ret[[2]] <- bquote(sim=rxTBSi(rx_pred_+sqrt(rx_r_) * rnorm(), rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  .ret
}


rxGetDistributionSimulationLines.t <- function(env, pred1) {
  .ret <- vector("list", 2)
  .ret[[1]] <- bquote(ipred=rxTBSi(rx_pred_, rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  .ret[[2]] <- bquote(sim=rxTBSi(rx_pred_+sqrt(rx_r_) * .(.getQuotedDistributionAndSimulationArgs(env, pred1)), rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  .ret
}


