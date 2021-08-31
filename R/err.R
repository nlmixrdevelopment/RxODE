## This is a list of supported distributions with the number of arguments they currently support.
.errDist <- list(
  "dpois" = 0,
  "dbinom" = 0:1,
  "dbern" = 0,
  "bern" = 0,
  "dbeta" = 2:3,
  ##
  ## "dnbinom"=2:3,  ## dnbinom is in R; FIXME: how does ot compare to dneg_binomial
  ## "dneg_binomial", ## not in base R (but in glnmm2)
  ##
  ## Available as external package http://ugrad.stat.ubc.ca/R/library/rmutil/html/BetaBinom.html
  ## "dbetabinomial", ## not in base R (but in glnmm2)
  "dt" = 1:2,
  "pois" = 0,
  "binom" = 0:1,
  "beta" = 2:3,
  "t" = 1:2,
  "add" = 1,
  "norm" = 1,
  "dnorm" = 1,
  "prop" = 1,
  "propT" = 1,
  "pow" = 2,
  "powT" = 2,
  "tbs" = 1,
  "boxCox" = 1,
  "tbsYj" = 1,
  "yeoJohnson" = 1,
  "logn" = 1,
  "dlogn" = 1,
  "logitNorm" = 1:3,
  "probitNorm" = 1:3,
  "lnorm" = 1,
  "dlnorm" = 1
)

.errDistsPositive <- c("add", "norm", "dnorm", "prop", "propT", "pow", "powT", "logn", "dlogn", "lnorm", "dlnorm", "logitNorm", "probitNorm")


.errUnsupportedDists <- c(
  "dchisq", "chisq", "dexp", "df", "f", "dgeom", "geom",
  "dhyper", "hyper", "dunif", "unif",
  "dweibull", "weibull",
  ## for testing...
  "nlmixrDist"
)

.errAddDists <- c("add", "prop", "propT", "norm", "pow", "powT", "dnorm", "logn", "lnorm", "dlnorm", "tbs", "tbsYj", "boxCox",
                  "yeoJohnson", "logitNorm", "probitNorm")


## the desired outcome for each expression is to capture the condition
## when the multiple endpoint occurs, the lower and upper for the
## transformation (if any) and the error type of the model.  The
## errors should be labeled "add" for an additive error and pow, pow2
## for the first and second argument of the power distribution.  The
## ini block should have already been processed by lotri, so we have
## some information about the model.  The predictions will be replaced
## by rxPred__# for easy replacement with the correct information for
## simulating or estimating in different nlmixr algorithms and
## different RxODE simulation routines

## Currently you can support the following types of expressions:
##
## f ~ add(add.sd)
##
## This assumes a named DVID of f
##
## or
##
## g ~ add(add.sd) | namedDvid
##
## In this case the named DVID is namedDvid
##
## In the final RxODE model expression these expressions will look like:
##
## rxPred__1 = f
## rxPred__2 = g
##
## With this change, there is sufficient information to send to the
## mu-referencing routine


.errHandleSingleDistributionTerm <- function(funName, expression, env) {

}

.errHandleSingleTerm <- function(funName, expression, env) {

}

##' Handle the error structure term
##'
##' @param expression The error structure term
##'
##' @param env The environment that holds information about the error
##'   structure
##'
##' @return
##' @author Matthew Fidler
.errHandleErrorStructure <- function(expression, env) {
  if (identical(expression[[1]], quote(`+`))) {
    env$isAnAdditiveExpression <- TRUE
    .errHandleErrorStructure(expression[[2]], env)
    .errHandleErrorStructure(expression[[3]], env)
  } else if (env$isAnAdditiveExpression) {
    .currErr <- deparse1(expression[[1]])
    if (.currErr %in% .errAddDists) {
      .errHandleSingleDistributionTerm(.currErr, expression, env)
    } else {
      .errHandleSingleTerm(.currErr, expression, env)
    }
  } else {
    .currErr <- deparse1(expression[[1]])
    if (.currErr %in% names(.errDist)) {
      .errHandleSingleDistributionTerm(.currErr, expression, env)
    } else {
      .errHandleSingleTerm(.currErr, expression, env)
    }
  }
}

##' Handle the right hand side conditional expression (if it exists)
##'
##' @param expression A right hand side of the tilde equation
##' @param env Environment for storing information about the expression
##' @return An expression without a conditional statement.
##'
##' @details
##'
##' In addition to stripping the conditional statement out of the
##' expression, the environment is modified when a conditional
##' expression is present.  First, the environment variable
##' `needsToBeAnErrorExpression` is changed to `TRUE`.  Second, the
##' expression `curCondition` is modified to match the information
##' within the conditional statement.
##'
##' @author Matthew Fidler
.errHandleCondition <- function(expression, env) {
  if (identical(expression[[1]], quote(`|`))) {
    env$needsToBeAnErrorExpression  <- TRUE
    env$curCondition <- deparse1(expression[[3]])
    return(expression[[2]])
  }
  expression
}

##' Handle the error expressions
##'
##' @param expression Single tilde error expression
##' @param env Environment with initial estimate data.frame
##' @return
##' @author Matthew Fidler
.errHandleTilde <- function(expression, env) {
  .left <- expression[[2]]
  env$curCondition <- deparse1(.left)
  env$needsToBeAnErrorExpression <- FALSE
  .right <- .errHandleCondition(expression[[3]], env)
  env$isAnAdditiveExpression <- FALSE
  .errHandleErrorStructure(.right, env)
}

##' Process the errors in the quoted expression
##'
##' @param x
##' @param df
##' @return
##' @author Matthew Fidler
##' @examples
##' lmat <- lotri({
##'  ## You may label each parameter with a comment
##'  tka <- 0.45 # Log Ka
##'  tcl <- log(c(0, 2.7, 100)) # Log Cl
##'  ## This works with interactive models
##'  ## You may also label the preceding line with label("label text")
##'  tv <- 3.45; label("log V")
##'  tvp <- 3.45; label("log V")
##'  cl.wt <- 0.1
##'  v.wt <- 0.1
##'  cl.sex <- 0.1
##'  v.sex <- 0.1
##'  cl.age <- 0.1
##'  v.age <- 0.1
##'  vp.wt <- 1
##'  vp.sex <- 1
##'  vp.age <- 1
##'  ## the label("Label name") works with all models
##'  eta.ka ~ 0.6
##'  eta.cl ~ 0.3
##'  eta.v ~ 0.1
##'  add.sd <- 0.7
##'})
##'
##'
##' iniDf <- as.data.frame(lmat)
##'
##' .errProcessExpression()
##' @noRd
.errProcessExpression <- function(x, df) {
  # ntheta neta1 neta2   name lower       est   upper   fix  err  label
  # backTransform condition trLow trHi
  .env <- new.env(parent=emptyenv())
  .env$top <- TRUE
  .env$df <- df
  # Add error structure like nlmixr ui had before transitioning to RxODE
  .env$df$err <- NA_character_
  .env$df$trLow <- .env$df$trHi <- NA_real_
  if (is.call(x)) {
    if (.env$top && identical(x[[1]], quote(`{`))) {
      .env$top <- FALSE
      .y <- x[-1]
      .ret <- character(length(.y))
      for (.i in seq_along(.y)) {
        if (identical(.y[[.i]][[1]], quote(`~`))) {
          .errHandleTilde(.y[[.i]], .env)
        } else {
          .ret[[.i]] <- deparse1(.y[[.i]])
        }
      }
      return(.ret)
    }
  }
  return(invisible(NULL))
}
