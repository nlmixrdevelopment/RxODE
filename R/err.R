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



##' Process a single expression and add information to the provided environment
##'
##' .. content for \details{} ..
##' @title
##' @param expression
##' @param env
##' @return
##' @author Matthew Fidler
##' @noRd
.errProcessExpression <- function(expression, env) {

}
