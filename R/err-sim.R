# This handles the errors for simulations
.createSimLineObject <- function(x, line) {
  .predDf <- get("predDf", x)
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine)
  class(.ret) <- c(paste(.predLine$distribution), "rxGetDistributionSimulationLines")
  .ret
}

#' This is a S3 method for getting the distribution lines for a RxODE simulation
#'
#' @param line  Parsed RxODE model environment
#' @return Lines for the simulation of `ipred` and `dv`. This is based on the idea that the focei parameters are defined
#' @author Matthew Fidler
#' @keywords internal
#' @export
rxGetDistributionSimulationLines <- function(line) {
  UseMethod("rxGetDistributionSimulationLines")
}

.simulationFun <- list(
  "t"="rt",
  "pois"="rpois",
  "binom"="rbinom",
  "beta"="rbeta",
  "chisq"="rchisq",
  "dexp"="rexp",
  "f"="rf",
  "geom"="rgeom",
  "hyper"="rhyper",
  "unif"="runif",
  "weibull"="rweibull",
  "ordinal"="rordinal"
)

.getQuotedDistributionAndSimulationArgs <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
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

#' @rdname rxGetDistributionSimulationLines
#' @export
rxGetDistributionSimulationLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .ret <- vector("list", 2)
  .ret[[1]] <- bquote(ipredSim <- rxTBSi(rx_pred_, rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  .ret[[2]] <- bquote(sim <- rxTBSi(rx_pred_+sqrt(rx_r_) * rnorm(), rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  c(.handleSingleErrTypeNormOrTFoceiBase(env, pred1), .ret)
}

#' @rdname rxGetDistributionSimulationLines
#' @export
rxGetDistributionSimulationLines.t <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .ret <- vector("list", 2)
  .ret[[1]] <- bquote(ipredSim <- rxTBSi(rx_pred_, rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  .ret[[2]] <- bquote(sim <- rxTBSi(rx_pred_+sqrt(rx_r_) * .(.getQuotedDistributionAndSimulationArgs(line)), rx_lambda_, rx_yj_, rx_low_, rx_hi_))
  c(.handleSingleErrTypeNormOrTFoceiBase(env, pred1), .ret)
}

#' @rdname rxGetDistributionSimulationLines
#' @export
rxGetDistributionSimulationLines.default <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  .ret <- vector("list", 1)
  .ret[[1]] <- bquote(sim <- .(.getQuotedDistributionAndSimulationArgs(line)))
  .ret
}

#' @rdname rxGetDistributionSimulationLines
#' @export
rxGetDistributionSimulationLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createSimLineObject(line, c)
    rxGetDistributionSimulationLines(.mod)
  })
}

#' Combine Error Lines and create RxODE expression
#'
#' @param uiModel UI model
#' @param errLines Error lines; If missing, get the error lines from
#'   `rxGetDistributionSimulationLines()`
#' @return quoted extression that can be evaluated to compiled RxODE
#'   model
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
#' @examples
#'
#'
#' one.cmt <- function() {
#'    ini({
#'      ## You may label each parameter with a comment
#'      tka <- 0.45 # Log Ka
#'      tcl <- log(c(0, 2.7, 100)) # Log Cl
#'      ## This works with interactive models
#'      ## You may also label the preceding line with label("label text")
#'      tv <- 3.45; label("log V")
#'      ## the label("Label name") works with all models
#'      eta.ka ~ 0.6
#'      eta.cl ~ 0.3
#'      eta.v ~ 0.1
#'      add.sd <- 0.7
#'    })
#'    model({
#'      ka <- exp(tka + eta.ka)
#'      cl <- exp(tcl + eta.cl)
#'       v <- exp(tv + eta.v)
#'       linCmt() ~ add(add.sd)
#'    })
#' }
#'
#' f <- RxODE(one.cmt)
#'
#' # You can get the simulation model easily by
#' rxCombineErrorLines(f)
#'
#' # You can then get the compiled model by simply evaluting the model:
#' r <- eval(rxCombineErrorLines(f))
#'
rxCombineErrorLines <- function(uiModel, errLines=NULL) {
  if(!inherits(uiModel, "rxUi")) {
    stop("uiModel must be a evaluated UI model by RxODE(modelFunction) or modelFunction()",
         call.=FALSE)
  }
  if (is.null(errLines)) {
    errLines <- rxGetDistributionSimulationLines(uiModel)
  }
  .lenLines <- sum(vapply(seq_along(errLines), function(i){
    length(errLines[[i]])
  }, integer(1)))
  .expr <- uiModel$lstExpr
  .predDf <- uiModel$predDf
  .lenLines <- .lenLines + length(uiModel$lstExpr) - length(.predDf$line)
  .ret <- vector("list", .lenLines + 1)
  .curErrLine <- 1
  .k <- 2
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.expr)) {
    if (.i %in% .predDf$line) {
      .curErr <- errLines[[.curErrLine]]
      for (.j in seq_along(.curErr)) {
        .ret[[.k]] <- .curErr[[.j]]
        .k <- .k + 1
      }
      .curErrLine <- .curErrLine + 1
    } else {
      .ret[[.k]] <- .expr[[.i]]
      .k <- .k + 1
    }
  }
  as.call(list(quote(`RxODE`), as.call(.ret)))
}
