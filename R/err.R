## This is a list of supported distributions with the number of arguments they currently support.
.errDist <- list(
  "dpois" = 0,
  "pois" = 0,
  "dbinom" = 0:1,
  "binom"=0:1,
  "dbern" = 0,
  "bern" = 0,
  "dbeta" = 2:3,
  "beta" = 2:3,
  "dt" = 1:2,
  "t" = 1:2,
  ##
  ## "dnbinom"=2:3,  ## dnbinom is in R; FIXME: how does ot compare to dneg_binomial
  ## "dneg_binomial", ## not in base R (but in glnmm2)
  ##
  ## Available as external package http://ugrad.stat.ubc.ca/R/library/rmutil/html/BetaBinom.html
  ## "dbetabinomial", ## not in base R (but in glnmm2)
  "add" = 1,
  "norm" = 1,
  "dnorm" = 1,
  "prop" = 1,
  "propT" = 1,
  "propF" = 2,
  "pow" = 2,
  "powT" = 2,
  "powF"=3,
  "tbs" = 1,
  "boxCox" = 1,
  "tbsYj" = 1,
  "yeoJohnson" = 1,
  "logn" = 1,
  "lnorm" = 1,
  "dlnorm" = 1,
  "dlogn" = 1,
  "logitNorm" = 1:3,
  "probitNorm" = 1:3,
  "combined1"=0,
  "combined2"=0,
  "comb1"=0,
  "comb2"=0
)

.errDistsPositive <- c("add", "norm", "dnorm", "prop", "propT", "pow", "powT", "logn", "dlogn", "lnorm", "dlnorm", "logitNorm", "probitNorm")


.errUnsupportedDists <- c(
  "dchisq", "chisq", "dexp", "df", "f", "dgeom", "geom",
  "dhyper", "hyper", "dunif", "unif",
  "dweibull", "weibull",
  ## for testing...
  "nlmixrDist"
)

.errAddDists <- c("add", "prop", "propT", "propF", "norm", "pow", "powT", "powF", "dnorm", "logn", "lnorm", "dlnorm", "tbs", "tbsYj", "boxCox",
                  "yeoJohnson", "logitNorm", "probitNorm", "combined1", "combined2", "comb1", "comb2")

.errIdenticalDists <- list(
  "add"=c("norm", "dnorm"),
  "lnorm"=c("logn", "dlogn", "dlnorm"),
  "boxCox"="tbs",
  "yeoJohnson"="tbsYj",
  "pois"="dpois",
  "binom"="dbinom",
  "bern"="dbern",
  "beta"="dbeta",
  "t"="dt",
  "combined1"="comb1",
  "combined2"="comb2"
)


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


##' Change distribution name to the preferred distribution name term
##'
##' This is determined by the internal preferred condition name list
##' `.errIdenticalDists`
##'
##' @param dist This is the input distribution
##' @return Preferred distribution term
##' @author Matthew Fidler
##'
##' @examples
##'
##' rxPreferredDistributionName("dnorm")
##'
##' rxPreferredDistributionName("add")
##'
##' # can be vectorized
##'
##' rxPreferredDistributionName(c("add","dnorm"))
##'
##' @export
rxPreferredDistributionName <- function(dist) {
  if (length(dist) == 1) {
    .names <- names(.errIdenticalDists)
    for(.n in .names) {
      if (dist == .n) return(.n)
      else if (dist %in% .errIdenticalDists[[.n]]) return(.n)
    }
    dist
  } else {
    vapply(dist, rxPreferredDistributionName, character(1))
  }
}

.incompatibleTransformations <-
  list(boxCox=c("yeoJohnson", "lnorm"),
       yeoJohnson=c("boxCox", "lnorm"),
       lnorm=c("yeoJohnson", "boxCox", "logit", "logitNorm", "probit", "probitNorm", "logit + yeoJohnson"),
       logit=c("probit", "probitNorm", "lnorm"),
       probit=c("logit", "lnorm"),
       "logit + yeoJohnson"=c("boxCox", "probit", "probitNorm", "lnorm"),
       "probit + yeoJohnson"=c("boxCox", "logit", "logitNorm", "lnorm"),
       "logit + boxCox"=c("yeoJohnson", "probit", "probitNorm", "lnorm"),
       "probit + boxCox"=c("yeoJohnson", "logit", "logitNorm", "lnorm"))

.rxTransformCombineLevels <- c(
  "boxCox", # 1
  "yeoJohnson", # 2
  "untransformed", # 3
  "lnorm", # 4
  "logit", # 5
  "logit + yeoJohnson", # 6
  "probit", # 7
  "probit + yeoJohnson", #8
  "logit + boxCox", # 9
  "probit + boxCox" # 10
)

.rxAddPropLevels <- c(
  "combined1", # 1
  "combined2", # 2
  "default" # 3
)

.incompatibleAddProp <- list(
  combined1="combined2",
  combined2="combined1"
)

.rxErrType <- c(
  "add", # 1
  "prop", # 2
  "pow", # 3
  "add + prop", # 4
  "add + pow", # 5
  "none" # 6
)
##' Demote the error type
##'
##' @param errType Error type factor
##' @return Demoted Error Type
##' @author Matthew Fidler
##' @export
##' @examples
##' rxDistributionCombine("add") %>%
##'   rxDistributionCombine("prop")
##'
##' # This removes the internal additive error
##' rxDistributionCombine("add") %>%
##'   rxDistributionCombine("prop") %>%
##'   rxDemoteAddErr()
##'
##' # This is used for logitNorm(NA), the additive portion is stripped
##' @keywords internal
rxDemoteAddErr <- function(errType) {
  if (inherits(errType, "factor")) {
    return(structure(switch(as.character(errType),
                            add=6L,
                            "add + prop"=2L,
                            "add + pow"=3L,
                            as.integer(errType)),
                     .Label =.rxErrType,
                     class="factor"))
  } else if (inherits(errType, "rxCombinedErrorList")) {
    return(.rxTransformCombineListOrChar(list(transform=errType$transform,
                                              errType=rxDemoteAddErr(errType$errType),
                                              errTypeF=errType$errTypeF,
                                              addProp=errType$addProp)))
  }
}

.incompatibleErrType <- list(prop=c("pow", "powT", "powF"),
                             propT=c("pow", "powT", "powF"),
                             propF=c("pow", "powT", "powF"),
                             pow=c("prop", "propT", "propF"),
                             powT=c("prop", "propT", "propF"),
                             powF=c("prop", "propT", "propF"))


.rxErrTypeF <- c(
  "untransformed", # 1
  "transformed", # 2
  "f", # 3
  "none")

.incompatibleErrTypeF <- list(prop=c("propT", "propF", "powT", "powF"),
                              propT=c("prop", "propF", "pow", "powF"),
                              propF=c("prop", "propT", "pow", "powT"),
                              pow=c("propF", "propT", "powF", "powT"),
                              powT=c("prop", "propF", "pow", "powF"),
                              powF=c("prop", "propT", "pow", "powT"))

.incompatibleErr <- function(err1, err2) {
  .errs <- sort(c(err1, err2))
   paste0("`", .errs[1], "` and `", .errs[2], "` are incompatible")
}
##' Combine error types to get the model F type
##'
##' @param newAddProp New error type
##' @param oldAddProp old error type
##' @return factor of the error type function OR a string indicating a syntax error
##' @author Matthew Fidler
##' @noRd
.rxCombineAddProp <- function(newAddProp, oldAddProp="default") {
  .tmp <- as.character(oldAddProp)
  .w <- which(names(.incompatibleAddProp) == .tmp)
  if (length(.w) == 1L) {
    if (newTransform %in% .incompatibleAddProp[[.w]]) {
      return(.incompatibleErr(.tmp, newTransform))
    }
  }
  structure(switch(newAddProp,
                   combined1 =1L,
                   combined2=2L,
                   default=3L,
                   ifelse(inherits(oldAddProp, "character"),
                          switch(oldAddProp,
                                 combined1 =1L,
                                 combined2=2L,
                                 default=3L,
                                 3L),
                          as.integer(oldAddProp))),
            .Label=.rxAddPropLevels,
            class="factor")
}

##' Combine error types to get the model F type
##'
##' @param newErrTypeF New error type
##' @param oldErrTypeF old error type
##' @return factor of the error type function OR a string indicating a syntax error
##' @author Matthew Fidler
##' @noRd
.rxCombineErrTypeF <- function(newErrTypeF, oldErrTypeF="none") {
  .tmp <- as.character(oldErrTypeF)
  .w <- which(names(.incompatibleErrTypeF) == .tmp)
  if (length(.w) == 1L) {
    if (newTransform %in% .incompatibleErrTypeF[[.w]]) {
      return(.incompatibleErr(.tmp, newTransform))
    }
  }
  structure(switch(newErrTypeF,
                   prop =1L,
                   propT=2L,
                   propF=3L,
                   pow=1L,
                   powT=2L,
                   powF=3L,
                   ifelse(inherits(oldErrTypeF, "character"),
                          switch(oldErrTypeF,
                                 prop =1L,
                                 propT=2L,
                                 propF=3L,
                                 pow=1L,
                                 powT=2L,
                                 powF=3L,
                                 4L),
                          as.integer(oldErrTypeF))),
            .Label=.rxErrTypeF,
            class="factor")
}


#' Combine error model
#'
#' @param newErrType New error portion combined with current proportion
#' @param oldErrType Old error information
#' @return A factor for error type OR a string with error information
#' @author Matthew Fidler
#' @noRd
.rxCombineErrType <- function(newErrType, oldErrType="none") {
  .tmp <- as.character(oldErrType)
  .w <- which(names(.incompatibleErrType) == .tmp)
  if (length(.w) == 1L) {
    if (newErrType %in% .incompatibleErrType[[.w]]) {
      return(.incompatibleErr(.tmp, newErrType))
    }
  }
  structure(switch(.tmp,
                   add=switch(newErrType,
                              prop=4L,
                              propT=4L,
                              propF=4L,
                              pow=5L,
                              powT=5L,
                              powF=5L,
                              1L),
                   prop=switch(newErrType,
                               add=4L,
                               probit=4L,
                               probitNorm=4L,
                               logit=4L,
                               logitNorm=4L,
                               lnorm=4L,
                               2L),
                   pow=switch(newErrType,
                              add=5L,
                              probit=5L,
                              probitNorm=5L,
                              logit=5L,
                              logitNorm=5L,
                              lnorm=5L,
                              3L),
                   none=structure(switch(newErrType,
                                         add=1L,
                                         lnorm=1L,
                                         logit=1L,
                                         logitNorm=1L,
                                         probit=1L,
                                         probitNorm=1L,
                                         prop=2L,
                                         propT=2L,
                                         propF=2L,
                                         pow=3L,
                                         powT=3L,
                                         powF=3L,
                                         6L)),
                   as.integer(oldErrType)),
            .Label=.rxErrType,
            class="factor")
}
##' Combine transformations
##'
##' @param newTransform New error structure added together
##' @param oldTransform Old transformation added together
##' @return A factor describing the transformation type
##' @author Matthew Fidler
##' @noRd
.rxCombineTransform <- function(newTransform, oldTransform="untransformed") {
  .tmp <- as.character(oldTransform)
  .w <- which(names(.incompatibleTransformations) == .tmp)
  if (length(.w) == 1L) {
    if (newTransform %in% .incompatibleTransformations[[.w]]) {
      return(.incompatibleErr(.tmp, newTransform))
    }
  }
  structure(switch(.tmp,
                   boxCox=switch(newTransform,
                                 logitNorm=9L,
                                 probitNorm=10L,
                                 1L),
                   yeoJohnson=switch(newTransform,
                                     logitNorm=6L,
                                     probitNorm=8L,
                                     2L),
                   logit=switch(newTransform,
                                boxCox=9L,
                                yeoJohnson=6L,
                                5L),
                   probit=switch(newTransform,
                                 boxCox=10L,
                                 yeoJohnson=8L,
                                 7L),
                   untransformed=structure(switch(newTransform,
                                                  boxCox=1L,
                                                  yeoJohnson=2L,
                                                  lnorm=4L,
                                                  logitNorm=5L,
                                                  probitNorm=7L,
                                                  3L),
                                           .Label=.rxTransformCombineLevels,
                                           class="factor"),
                   as.integer(oldTransform)),
            .Label=.rxTransformCombineLevels,
            class="factor")
}

.rxTransformCombineListOrChar <- function(inputList) {
  if (inherits(inputList, "character")) return(inputList)
  .err  <- NULL
  for (i in names(inputList)) {
    if (inherits(inputList[[i]], "character")) {
      .err <- c(.err, inputList[[i]])
    }
  }
  if (is.null(.err)) {
    .ret <- inputList
    class(.ret) <- "rxCombinedErrorList"
    return(.ret)
  } else {
    return(paste(.err, collapse="\n"))
  }
}


##' Combine transformations and error structures
##'
##' Combine error information to figure out what transformation is
##' being applied for the current endpoint
##'
##'
##' @param oldDistribution This is the old transformation, by default is
##'   zero representing no prior transformation. This parameter is
##'   first to allow piping. When the parameter `addTransform` is
##'   missing and `oldDistribution` is a character value, this functions
##'   swaps `oldDistribution` and `addTransform` and assigns
##'   `oldDistribution` to zero assuming that there is no prior
##'   distribution.
##'
##' @param newTransform This is the new distribution that is being
##'   "added" to the current transformation.  These assumes the inputs
##'   are in the preferred distribution name, as determined by
##'   `rxPreferredDistributionName()`
##'
##'
##' @return The new transformation as a factor
##'
##' @author Matthew Fidler
##'
##' @examples
##'
##' rxDistributionCombine("probitNorm")
##'
##' rxDistributionCombine("probitNorm") %>%
##'   rxDistributionCombine("boxCox")
##'
##'
##' @export
##' @keywords internal
rxDistributionCombine <- function(oldDistribution, newDistribution) {
  if (missing(newDistribution) && inherits(oldDistribution, "character")) {
    return(.rxTransformCombineListOrChar(list(transform=.rxCombineTransform(oldDistribution),
                                              errType=.rxCombineErrType(oldDistribution),
                                              errTypeF=.rxCombineErrTypeF(oldDistribution),
                                              addProp=.rxCombineAddProp(oldDistribution))))
  } else if (inherits(oldDistribution, "rxCombinedErrorList")) {
    return(.rxTransformCombineListOrChar(list(transform=.rxCombineTransform(newDistribution, oldDistribution$transform),
                                              errType=.rxCombineErrType(newDistribution, oldDistribution$errType),
                                              errTypeF=.rxCombineErrTypeF(newDistribution, oldDistribution$errTypeF),
                                              addProp=.rxCombineAddProp(newDistribution, oldDistribution$addProp))))
  } else {
    stop("old transform not in the proper format", call.=FALSE)
  }
}

##' Checks to see if an expression is numeric
##'
##' @param expression quoted expression
##' @param env Environment to store result in `env$.numeric`
##' @return TRUE if this is an expression containing a positive or negative expression or FALSE if it is an expression that doesn't contain an expression.
##' @author Matthew Fidler
##' @noRd
.is.numeric <- function(expression, env) {
  if (is.numeric(expression)) {
    env$.numeric <- expression
    return(TRUE)
  } else if (length(expression) == 2L) {
    if (identical(expression[[1]], quote(`-`)) &&
          is.numeric(expression[[2]])) {
      env$.numeric <- -(expression[[2]])
      return(TRUE)
    } else if (identical(expression[[1]], quote(`+`)) &&
          is.numeric(expression[[2]])) {
      env$.numeric <- expression[[2]]
      return(TRUE)
    }
  }
  return(FALSE)
}

##' This handles the error distribution for a single argument.
##'
##' @param argumentNumber The argument number of the distribution being processed
##' @param funName Function name string of the distribution name
##' @param expression Function expression (including the function name)
##' @param env Environment where the names and calculations are made for the user interface.
##' @return None, called for the side effects.
##' @author Matthew Fidler
##' @noRd
.errHandleSingleDistributionArgument <- function(argumentNumber, funName, expression, env) {
  .cur <- expression[[argumentNumber + 1]]
  .isLogitOrProbit <- (funName %in% c("logitNorm", "probitNorm"))
  if (is.name(.cur)) {
    .curName <- as.character(.cur)
    .w <- which(env$df$name == .curName)
    if (.isLogitOrProbit && argumentNumber > 2) {
      env$err <- c(env$err,
                   paste0("`", funName, "()` requires numeric bounds"))
    } else if (length(.w) == 1L) {
      .df  <- env$df
      .df$err[.w] <- ifelse(argumentNumber == 1, funName, paste0(funName, argumentNumber))
      .df$condition[.w] <- rxPreferredDistributionName(env$curCondition)
      assign("df", .df, envir=env)
      assign("lastDistAssign", .curName, envir=env)
    }
  } else if (.is.numeric(.cur, env)) {
    if (.isLogitOrProbit) {
      .curName <- env$lastDistAssign
      .w <- which(env$df$name == .curName)
      if (length(.w) == 1L) {
        env$trLimit[argumentNumber - 1] <- env$.numeric
        #env$df[.w, c("trHi", "trLow")[argumentNumber - 1]] <- env$.numeric
      }
    }
  } else {
  }
}

##' Checks and modifies the error term as needed in the user interface function.
##'
##' @param funName Name of the distributional error
##' @param expression The R expression (including the function name)
##' @param env The environment that holds information about the error structure
##' @return Nothing; Called for the side effects
##' @author Matthew Fidler
##' @noRd
.errHandleSingleDistributionTerm <- function(funName, expression, env) {
  .nargs <- length(expression) - 1
  .errDistArgs <- .errDist[[funName]]
  if (.nargs %in% .errDistArgs) {
    if (.nargs > 0) {
      lapply(seq(1, .nargs), .errHandleSingleDistributionArgument, funName=funName, expression=expression, env=env)
    }
    env$distInfo <- rxDistributionCombine(env$distInfo, funName)
    env$needsToBeAnErrorExpression <- TRUE
  } else {
    .min <- range(.errDistArgs)
    .max <- .min[2]
    .min <- .min[1]
    if (.min == .max) {
      assign("err", c(env$err,
                      paste0("`", funName, "` requires ",
                             .max, " argument(s), you specified ", .nargs)), envir=env)
    } else {
      assign("err", c(env$err,
                      paste0("`", funName, "` requires ",
                             .min, " to ", .max, " argument(s), you specified ", .nargs)),
             envir=env)
    }
  }
}
##' This handles a function that is not an error term.
##'
##' @param funName Function that is being called, as a character
##' @param expression Expression (including function)
##' @param env Environment that stores the information about errors
##' @return Nothing, called for the side effects
##' @author Matthew Fidler
##' @noRd
.errHandleSingleTerm <- function(funName, expression, env) {
  env$hasNonErrorTerm <- TRUE
}

##' Handle the error structure term
##'
##' @param expression The error structure term
##'
##' @param env The environment that holds information about the error
##'   structure
##'
##' @return Nothing, called for the side efects
##' @author Matthew Fidler
##' @noRd
.errHandleErrorStructure <- function(expression, env) {
  if (identical(expression[[1]], quote(`+`))) {
    env$isAnAdditiveExpression <- TRUE
    .errHandleErrorStructure(expression[[2]], env)
    .errHandleErrorStructure(expression[[3]], env)
  } else if (env$isAnAdditiveExpression) {
    .currErr <- deparse1(expression[[1]])
    if (.currErr %in% .errAddDists) {
      .errHandleSingleDistributionTerm(.currErr, expression, env)
    } else if (.currErr %in% names(.errDist)) {
      assign("err", c(env$err,
                      paste0("`", .currErr, "` is incorrectly added to an error expression")), envir=env)
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
##' @noRd
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
##' @noRd
.errHandleTilde <- function(expression, env) {
  .left <- expression[[2]]
  env$trLimit <- c(-Inf, Inf)
  env$curCondition <- env$curVar <- deparse1(.left)
  env$hasNonErrorTerm <- FALSE
  env$needsToBeAnErrorExpression <- FALSE
  .right <- .errHandleCondition(expression[[3]], env)
  env$isAnAdditiveExpression <- FALSE
  env$distInfo <- rxDistributionCombine("")
  .errHandleErrorStructure(.right, env)
  if (inherits(env$distInfo, "character")) {
    env$err <- c(env$err, env$distInfo)
    env$distInfo <- rxDistributionCombine("")
  }
  env$predDf <- rbind(env$predDf,
                      data.frame(cond=env$curCondition, var=env$curVar, dvid=env$curDvid,
                                 trHi=env$trLimit[1], trLow=env$trLimit[2],
                                 transform=env$distInfo$transform,
                                 errType=env$distInfo$errType,
                                 errTypeF=env$distInfo$errTypeF,
                                 addProp=env$distInfo$addProp))
  env$curDvid <- env$curDvid + 1
  if (env$hasNonErrorTerm & env$needsToBeAnErrorExpression) {
    assign("err", c(env$err, "syntax error: cannot mix additive expression with algebraic expressions"),
           envir=env)
  }
}

##' Process the errors in the quoted expression
##'
##' @param x Quoted expression for parsing
##' @param df lotri data.frame of estimates
##' @return Environment with error information setup.
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
  .env$err <- NULL
  # Add error structure like nlmixr ui had before transitioning to RxODE
  .env$df$err <- NA_character_
  #.env$df$trLow <- .env$df$trHi <- NA_real_
  .env$curDvid <- 1
  # Pred df needs to be finalized with compartment information from parsing the raw RxODE model
  .env$predDf  <- NULL
  .env$lastDistAssign <- ""
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
      if (!is.null(.env$err)) {
        stop(paste(c("Syntax Errors:", paste(" ", .env$err)), collapse="\n"),
             call.=FALSE)
      }
      ## print(.env$df)
      ## print(.env$predDf)
      return(.ret)
    }
  }
  return(invisible(NULL))
}
