## This is a list of supported distributions with the number of arguments they currently support.
.errDist <- list(
  "dpois" = 1,
  "pois" = 1,
  "dbinom" = 1:2,
  "binom"=1:2,
  "dbern" = 1,
  "bern" = 1,
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
  "comb2"=0,
  "dchisq"=1:2,
  "chisq"=1:2,
  "dexp"=0:1,
  "df"=2:3,
  "f"=2:3,
  "dgeom"=1,
  "geom"=1,
  "dhyper"=3,
  "hyper"=3,
  "dunif"=0:2,
  "unif"=0:2,
  "dweibull"=1:2,
  "weibull"=1:2
)

.errDistsPositive <- c("add", "norm", "dnorm", "prop", "propT", "pow", "powT", "logn", "dlogn", "lnorm", "dlnorm", "logitNorm", "probitNorm")


.errUnsupportedDists <- c(
  ## for testing...
  "nlmixrDist"
)

.errAddDists <- c("add", "prop", "propT", "propF", "norm", "pow", "powT", "powF", "dnorm", "logn", "lnorm", "dlnorm", "tbs", "tbsYj", "boxCox",
                  "yeoJohnson", "logitNorm", "probitNorm", "combined1", "combined2", "comb1", "comb2", "t")

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
  "combined2"="comb2",
  "chisq"="dchisq", #6
  "f"="df", #8
  "geom"="dgeom", #9
  "hyper"="dhyper", #10
  "unif"="dunif", #11
  "weibull"="dweibull"#12
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


#' Change distribution name to the preferred distribution name term
#'
#' This is determined by the internal preferred condition name list
#' `.errIdenticalDists`
#'
#' @param dist This is the input distribution
#' @return Preferred distribution term
#' @author Matthew Fidler
#'
#' @examples
#'
#' rxPreferredDistributionName("dnorm")
#'
#' rxPreferredDistributionName("add")
#'
#' # can be vectorized
#'
#' rxPreferredDistributionName(c("add","dnorm"))
#'
#' @export
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

.rxTransformHasALambdaParameter <- function(distribution) {
  (as.integer(distribution) %in% c(1L, 2L, 6L, 8L, 9L, 10L))
}

.rxTransformHasBounds <- function(distribution) {
    (as.integer(distribution) %in% 5:10)
}

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

.rxDistributionType <- c(
  "norm", #1
  "pois", #2
  "binom", #3
  "beta", #4
  "t", #5
  "chisq", #6
  "dexp", #7
  "f", #8
  "geom", #9
  "hyper", #10
  "unif", #11
  "weibull", #12
  "ordinal", # 13
  "-2LL" #14
)

#' Demote the error type
#'
#' @param errType Error type factor
#' @return Demoted Error Type
#' @author Matthew Fidler
#' @export
#' @examples
#' rxErrTypeCombine("add") %>%
#'   rxErrTypeCombine("prop")
#'
#' # This removes the internal additive error
#' rxErrTypeCombine("add") %>%
#'   rxErrTypeCombine("prop") %>%
#'   rxDemoteAddErr()
#'
#' # This is used for logitNorm(NA), the additive portion is stripped
#' @keywords internal
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
#' Combine error types to get the model F type
#'
#' @param newAddProp New error type
#' @param oldAddProp old error type
#' @return factor of the error type function OR a string indicating a syntax error
#' @author Matthew Fidler
#' @noRd
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

#' Combine error types to get the model F type
#'
#' @param newErrTypeF New error type
#' @param oldErrTypeF old error type
#' @return factor of the error type function OR a string indicating a syntax error
#' @author Matthew Fidler
#' @noRd
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
#' Combine transformations
#'
#' @param newTransform New error structure added together
#' @param oldTransform Old transformation added together
#' @return A factor describing the transformation type
#' @author Matthew Fidler
#' @noRd
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

#' This is a wrapper to make sure that the transformation combination returns the correct value
#'
#' @param inputList This is either an input list or character vector of length one
#' @return Either a complete list, or a character vector which represents the parsed error that was encountered
#' @author Matthew Fidler
#' @noRd
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


#' Combine transformations and error structures
#'
#' Combine error information to figure out what transformation is
#' being applied for the current endpoint
#'
#'
#' @param oldErrType This is the old transformation, by default is
#'   zero representing no prior transformation. This parameter is
#'   first to allow piping. When the parameter `addTransform` is
#'   missing and `oldErrType` is a character value, this functions
#'   swaps `oldErrType` and `addTransform` and assigns
#'   `oldErrType` to zero assuming that there is no prior
#'   distribution.
#'
#' @param newTransform This is the new distribution that is being
#'   "added" to the current transformation.  These assumes the inputs
#'   are in the preferred distribution name, as determined by
#'   `rxPreferredDistributionName()`
#'
#'
#' @return The new transformation as a factor
#'
#' @author Matthew Fidler
#'
#' @examples
#'
#' rxErrTypeCombine("probitNorm")
#'
#' rxErrTypeCombine("probitNorm") %>%
#'   rxErrTypeCombine("boxCox")
#'
#'
#' @export
#' @keywords internal
rxErrTypeCombine <- function(oldErrType, newErrType) {
  if (missing(newErrType) && inherits(oldErrType, "character")) {
    return(.rxTransformCombineListOrChar(list(transform=.rxCombineTransform(oldErrType),
                                              errType=.rxCombineErrType(oldErrType),
                                              errTypeF=.rxCombineErrTypeF(oldErrType),
                                              addProp=.rxCombineAddProp(oldErrType))))
  } else if (inherits(oldErrType, "rxCombinedErrorList")) {
    return(.rxTransformCombineListOrChar(list(transform=.rxCombineTransform(newErrType, oldErrType$transform),
                                              errType=.rxCombineErrType(newErrType, oldErrType$errType),
                                              errTypeF=.rxCombineErrTypeF(newErrType, oldErrType$errTypeF),
                                              addProp=.rxCombineAddProp(newErrType, oldErrType$addProp))))
  } else {
    stop("old transform not in the proper format", call.=FALSE)
  }
}

#' Checks to see if an expression is numeric
#'
#' @param expression quoted expression
#' @param env Environment to store result in `env$.numeric`
#' @return TRUE if this is an expression containing a positive or negative expression or FALSE if it is an expression that doesn't contain an expression.
#' @author Matthew Fidler
#' @noRd
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



.allowDemoteAddDistributions <- c("lnorm", "probitNorm", "logitNorm")

.namedArgumentsToPredDf <- list(
  add="a",
  lnorm="a",
  boxCox="lambda",
  yeoJohnson="lambda",
  pow=c("b", "c"),
  powT=c("b", "c"),
  powF=c("b", "c", "f"),
  prop="b",
  propT="b",
  propF=c("b", "f"),
  t=c("d", "e")
)

#' This handles the error distribution for a single argument.
#'
#' @param argumentNumber The argument number of the distribution being processed
#' @param funName Function name string of the distribution name
#' @param expression Function expression (including the function name)
#' @param env Environment where the names and calculations are made for the user interface.
#' @return None, called for the side effects.
#' @author Matthew Fidler
#' @noRd
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
    } else {
      .w <- which(names(.namedArgumentsToPredDf) == funName)
      if (length(.w) == 1) {
        .curLst <- .namedArgumentsToPredDf[[.w]]
        if (argumentNumber <= length(.curLst)) {
          assign(.curLst[argumentNumber], .curName, envir=env)
          return(NULL)
        }
      }
      env$err <- c(env$err,
                   paste0("In the error expression, the variable `", .curName, "` must be estimated, not calculated"))
    }
  } else if (.is.numeric(.cur, env)) {
    if (.isLogitOrProbit) {
      .curName <- env$lastDistAssign
      .w <- which(env$df$name == .curName)
      if (length(.w) == 1L) {
        env$trLimit[argumentNumber - 1] <- env$.numeric
      }
    }
  } else if (is.logical(.cur)) {
    if (is.na(.cur)) {
      if (argumentNumber == 1 && funName %in% .allowDemoteAddDistributions) {
        env$needToDemoteAdditiveExpression <- TRUE
      } else {
        env$err <- c(env$err,
                     paste0("NA in `", funName, "()` cannot be used here"))
      }
    }
  }
}

#' Checks and modifies the error term as needed in the user interface function.
#'
#' @param funName Name of the distributional error
#' @param expression The R expression (including the function name)
#' @param env The environment that holds information about the error structure
#' @return Nothing; Called for the side effects
#' @author Matthew Fidler
#' @noRd
.errHandleSingleDistributionTerm <- function(funName, expression, env) {
  .nargs <- length(expression) - 1
  .doIt <- FALSE
  if (funName == "ordinal") {
    if (.nargs == 0) {
      assign("err", c(env$err,
                      paste0("ordinal errors require at least 1 argument (i.e. err ~ c(err1))", envir=env)))
      return(invisible())
    }
    .doIt <- TRUE
  } else {
    .errDistArgs <- .errDist[[funName]]
    .doIt <- .nargs %in% .errDistArgs
  }
  if (.doIt) {
    if (.nargs > 0) {
      lapply(seq(1, .nargs), .errHandleSingleDistributionArgument, funName=funName, expression=expression, env=env)
    }
    env$errTypeInfo <- rxErrTypeCombine(env$errTypeInfo, funName)
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
#' This handles a function that is not an error term.
#'
#' @param funName Function that is being called, as a character
#' @param expression Expression (including function)
#' @param env Environment that stores the information about errors
#' @return Nothing, called for the side effects
#' @author Matthew Fidler
#' @noRd
.errHandleSingleTerm <- function(funName, expression, env) {
  env$hasNonErrorTerm <- TRUE
}

#' Handle the error structure term
#'
#' @param expression The error structure term
#'
#' @param env The environment that holds information about the error
#'   structure
#'
#' @return Nothing, called for the side efects
#' @author Matthew Fidler
#' @noRd
.errHandleErrorStructure <- function(expression, env) {
  if (identical(expression[[1]], quote(`+`))) {
    env$isAnAdditiveExpression <- TRUE
    .errHandleErrorStructure(expression[[2]], env)
    .errHandleErrorStructure(expression[[3]], env)
  } else if (env$isAnAdditiveExpression) {
    .currErr <- rxPreferredDistributionName(deparse1(expression[[1]]))
    if (.currErr %in% .errAddDists) {
      if (.currErr == "t") env$distribution <- "t"
      .errHandleSingleDistributionTerm(.currErr, expression, env)
    } else if (.currErr %in% names(.errDist)) {
      assign("err", c(env$err,
                      paste0("`", .currErr, "` is incorrectly added to an error expression")), envir=env)
    } else {
      .errHandleSingleTerm(.currErr, expression, env)
    }
  } else {
    .currErr <- deparse1(expression[[1]])
    if (.currErr == "c") {
      env$distribution <- "ordinal"
      .errHandleSingleDistributionTerm("ordinal", expression, env)
    } else if (.currErr %in% names(.errDist)) {
      .currErr <- rxPreferredDistributionName(.currErr)
      env$distribution <- .currErr
      .errHandleSingleDistributionTerm(.currErr, expression, env)
    } else {
      .errHandleSingleTerm(.currErr, expression, env)
    }
  }
}

#' Handle the right hand side conditional expression (if it exists)
#'
#' @param expression A right hand side of the tilde equation
#' @param env Environment for storing information about the expression
#' @return An expression without a conditional statement.
#'
#' @details
#'
#' In addition to stripping the conditional statement out of the
#' expression, the environment is modified when a conditional
#' expression is present.  First, the environment variable
#' `needsToBeAnErrorExpression` is changed to `TRUE`.  Second, the
#' expression `curCondition` is modified to match the information
#' within the conditional statement.
#'
#' @author Matthew Fidler
#' @noRd
.errHandleCondition <- function(expression, env) {
  if (identical(expression[[1]], quote(`|`))) {
    env$needsToBeAnErrorExpression  <- TRUE
    env$curCondition <- deparse1(expression[[3]])
    return(expression[[2]])
  }
  expression
}

#' Handle -2LL equivalent for n2ll or `linCmt()` statements for lhs
#'
#' @param expression Left handed side of the equation
#'
#' @param env Environment for storing information about the error
#'   structure
#'
#' @return Nothing called for side effects
#'
#' @details
#'
#' The takes the expression
#'
#' n2ll(var) ~  log(...)
#'
#' And strips the `n2ll` and sets the flag `env$n2ll` to `TRUE`
#'
#' Otherwise it leaves the `expression` alone and returns the value
#'
#' This takes the expression:
#'
#' linCmt() ~ add() ...
#'
#' And returns a `rxLinCmt` as the lhs value and sets the linCmt flag to TRUE
#'
#'
#' @author Matthew Fidler
#' @noRd
.errHandleN2llOrLinCmt <- function(expression, env) {
  if (is.call(expression)) {
    if (identical(expression[[1]], quote(`n2ll`)) &&
          length(expression) == 2L) {
      env$n2ll <- TRUE
      return(expression[[2]])
    } else if (identical(expression[[1]], quote(`linCmt`))) {
      env$linCmt <- TRUE
      return(quote(`rxLinCmt`))
    }
  }
  expression
}


#' Handle the error expressions
#'
#' @param expression Single tilde error expression
#' @param env Environment with initial estimate data.frame
#' @return
#' @author Matthew Fidler
#' @noRd
.errHandleTilde <- function(expression, env) {
  env$n2ll <- FALSE
  env$linCmt <- FALSE
  .left <- .errHandleN2llOrLinCmt(expression[[2]], env)
  env$trLimit <- c(-Inf, Inf)
  env$a <- env$b <- env$c <- env$d <- env$e <- env$f <- env$lambda <- NA_character_
  env$curCondition <- env$curVar <- deparse1(.left)
  env$hasNonErrorTerm <- FALSE
  env$needsToBeAnErrorExpression <- FALSE
  env$needToDemoteAdditiveExpression <- FALSE
  .right <- .errHandleCondition(expression[[3]], env)
  env$isAnAdditiveExpression <- FALSE
  env$errTypeInfo <- rxErrTypeCombine("")
  env$distribution <- "norm"
  .errHandleErrorStructure(.right, env)
  if (inherits(env$errTypeInfo, "character")) {
    env$err <- c(env$err, env$errTypeInfo)
    env$errTypeInfo <- rxErrTypeCombine("")
  } else if (env$needToDemoteAdditiveExpression) {
    env$errTypeInfo <- rxDemoteAddErr(env$errTypeInfo)
  }
  if (env$hasNonErrorTerm & env$needsToBeAnErrorExpression) {
    assign("errGlobal", c(env$errGlobal, "cannot mix error expression with algebraic expressions"),
           envir=env)
  } else if (env$hasNonErrorTerm) {
    if (env$n2ll) {
      env$distribution <- "-2LL"
      env$predDf <- rbind(env$predDf,
                          data.frame(cond=env$curCondition, var=env$curVar, dvid=env$curDvid,
                                     trHi=env$trLimit[1], trLow=env$trLimit[2],
                                     transform=env$errTypeInfo$transform,
                                     errType=env$errTypeInfo$errType,
                                     errTypeF=env$errTypeInfo$errTypeF,
                                     addProp=env$errTypeInfo$addProp,
                                     distribution=factor(env$distribution, levels=.rxDistributionType),
                                     line=env$line,
                                     a=env$a,
                                     b=env$b,
                                     c=env$c,
                                     d=env$d,
                                     e=env$e,
                                     f=env$f,
                                     lambda=env$lambda,
                                     linCmt=env$linCmt))
      env$curDvid <- env$curDvid + 1L

    }
  } else if (!env$hasNonErrorTerm) {
    if (env$n2ll) {
      assign("errGlobal", c(env$errGlobal, "a -2 log-likelihood expression cannot use abbreviated error codes like add() + prop() "),
             envir=env)
    } else {
      if (!(env$distribution %in% .rxDistributionType)) {
        env$distribution <- "norm"
      }
      env$predDf <- rbind(env$predDf,
                          data.frame(cond=env$curCondition, var=env$curVar, dvid=env$curDvid,
                                     trHi=env$trLimit[1], trLow=env$trLimit[2],
                                     transform=env$errTypeInfo$transform,
                                     errType=env$errTypeInfo$errType,
                                     errTypeF=env$errTypeInfo$errTypeF,
                                     addProp=env$errTypeInfo$addProp,
                                     distribution=factor(env$distribution, levels=.rxDistributionType),
                                     line=env$line,
                                     a=env$a,
                                     b=env$b,
                                     c=env$c,
                                     d=env$d,
                                     e=env$e,
                                     f=env$f,
                                     lambda=env$lambda,
                                     linCmt=env$linCmt))
      env$curDvid <- env$curDvid + 1L
    }

  }
}

#' Process the errors in the quoted expression
#'
#' @param x Quoted expression for parsing
#' @param df lotri data.frame of estimates
#' @return Environment with error information setup.
#' @author Matthew Fidler
#' @examples
#' lmat <- lotri({
#'  ## You may label each parameter with a comment
#'  tka <- 0.45 # Log Ka
#'  tcl <- log(c(0, 2.7, 100)) # Log Cl
#'  ## This works with interactive models
#'  ## You may also label the preceding line with label("label text")
#'  tv <- 3.45; label("log V")
#'  tvp <- 3.45; label("log V")
#'  cl.wt <- 0.1
#'  v.wt <- 0.1
#'  cl.sex <- 0.1
#'  v.sex <- 0.1
#'  cl.age <- 0.1
#'  v.age <- 0.1
#'  vp.wt <- 1
#'  vp.sex <- 1
#'  vp.age <- 1
#'  ## the label("Label name") works with all models
#'  eta.ka ~ 0.6
#'  eta.cl ~ 0.3
#'  eta.v ~ 0.1
#'  add.sd <- 0.7
#'})
#'
#'
#' iniDf <- as.data.frame(lmat)
#'
#' .errProcessExpression()
#' @noRd
.errProcessExpression <- function(x, ini) {
  # ntheta neta1 neta2   name lower       est   upper   fix  err  label
  # backTransform condition trLow trHi
  .env <- new.env(parent=emptyenv())
  .env$eta <- dimnames(ini)[[1]]
  .env$top <- TRUE
  .env$df <- as.data.frame(ini)
  .env$err <- NULL
  .env$errGlobal <- NULL
  # Add error structure like nlmixr ui had before transitioning to RxODE
  .env$df$err <- NA_character_
  #.env$df$trLow <- .env$df$trHi <- NA_real_
  .env$curDvid <- 1L
  # Pred df needs to be finalized with compartment information from parsing the raw RxODE model
  .env$predDf  <- NULL
  .env$lastDistAssign <- ""
  if (is.call(x)) {
    if (.env$top && identical(x[[1]], quote(`{`))) {
      .env$top <- FALSE
      .y <- x[-1]
      .env$lstChr <- character(length(.y))
      .env$lstErr <- vector(length(.y), mode="list")
      .env$lstExpr <- vector(length(.y), mode="list")
      .env$hasErrors <- FALSE
      for (.i in seq_along(.y)) {
        .env$line <- .i
        if (identical(.y[[.i]][[1]], quote(`~`))) {
          .errHandleTilde(.y[[.i]], .env)
        }
        .env$lstChr[[.i]] <- deparse1(.y[[.i]])
        .env$lstExpr[[.i]] <- .y[[.i]]
        if (!is.null(.env$err)) {
          .env$lstErr[[.i]] <- paste(paste(" ", .env$err), collapse="\n")
          .env$err <- NULL
          .env$hasErrors <- TRUE
        }
      }
      .env$iniDf <- .env$df
      if (is.null(.env$predDf)) {
        .env$errGlobal <- c(.env$errGlobal,
                            "there must be at least one prediction in the model({}) block.  Use `~` for predictions")
      } else if (length(.env$predDf[, 1]) == 0L) {
        .env$errGlobal <- c(.env$errGlobal,
                            "there must be at least one prediction in the model({}) block.  Use `~` for predictions")
      }
      if (!is.null(.env$errGlobal)) {
        stop(paste(.env$errGlobal, collapse="\n"), call.=FALSE)
      }
      .env$mv0 <- rxModelVars(paste(.env$lstChr[-.env$predDf$line], collapse="\n"))
      .env$curCmt <- length(.env$mv0$state)
      .env$extraCmt <- NULL
      .env$predDf$cmt <- vapply(seq_along(.env$predDf$line),
                                function(i) {
                                  .cmtName <- .env$predDf$cond[i]
                                  .w <- which(.cmtName == .env$mv0$state)
                                  if (length(.w) == 1L) return(.w)
                                  assign("curCmt", .env$curCmt + 1L, envir=.env)
                                  assign("extraCmt", c(.env$extraCmt,
                                                       paste0("cmt(", .cmtName, ")")),
                                         envir=.env)
                                  .env$curCmt
                                }, integer(1))

      .env$extraDvid <- paste0("dvid(", paste(.env$predDf$cmt, collapse = ","), ")")
      # Cleanup the environment
      .rm <- intersect(c("curCondition", "curDvid", "curVar", "df",
                         "errTypeInfo", "err", "hasNonErrorTerm", "isAnAdditiveExpression",
                         "lastDistAssign", "line", "needsToBeAnErrorExpression", "needToDemoteAdditiveExpression",
                         "top", "trLimit", ".numeric", "a", "b", "c", "d", "e", "f",  "lambda",
                         "curCmt", "errGlobal", "linCmt", "n2ll", "distribution"),
                       ls(envir=.env, all.names=TRUE))
      if (length(.rm) > 0) rm(list=.rm, envir=.env)
      return(.env)
    }
  }
  return(invisible(NULL))
}
#' Turn a character expression into quoted symbol
#'
#' @param chr Character symbol
#' @return Quoted symbol
#' @author Matthew Fidler
#' @noRd
.enQuote <- function(chr) {
  eval(parse(text=paste0("quote(", chr, ")")))
}

#' Get the lambda value based on the pred information
#'
#' @param env Environment that has the environment
#' @param pred1 Single error data frame
#' @return Lambda expression
#' @author Matthew Fidler
#' @noRd
.rxGetLambdaFromPred1AndIni <- function(env, pred1) {
  if (!is.na(pred1$lambda)) {
    return(.enQuote(pred1$lambda))
  }
  if (.rxTransformHasALambdaParameter(pred1$transform)) {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err %in% c("boxCox", "yeoJohnson") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      return(.enQuote(env$iniDf$name[.w]))
    } else {
      stop("cannot find lambda", call.=FALSE)
    }
  }
  return(1)
}
#' Get the lower boundary condition when the transformation requires it
#'
#' @param env  Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return Lower Boundary
#' @author Matthew Fidler
#' @noRd
.rxGetLowBoundaryPred1AndIni <- function(env, pred1) {
  if (.rxTransformHasBounds(pred1$transform)) {
    return(pred1$trLow)
  }
  return(0)
}
#' Get the upper boundary condition when the transformation it
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return Upper Boundary
#' @author Matthew Fidler
#' @noRd
.rxGetHiBoundaryPred1AndIni <- function(env, pred1) {
  if (.rxTransformHasBounds(pred1$transform)) {
    return(pred1$trHi)
  }
  return(1)
}

#' Get the prediction name
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return The prediction symbol
#' @author Matthew Fidler
#' @noRd
.rxGetPredictionF <- function(env, pred1) {
  .f <- pred1$var
  if (.f == "rxLinCmt") {
    return(quote(linCmt()))
  }
  .enQuote(.f)
}

#' Get the prediction transformation
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @param yj The transformation number for the current error
#' @return The transformation expression
#' @author Matthew Fidler
#' @noRd
.rxGetPredictionFTransform <- function(env, pred1, yj) {
  if (yj == 2) {
    return(quote(rx_pred_f_))
  } else if (yj == 3) {
    return(quote(log(rx_pred_f_)))
  } else {
    return(quote(rxTBS(rx_pred_f_, rx_lambda_, rx_yj_, rx_low_, rx_hi_)))
  }
}

#' Get the additive transformation
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return The quoted symbolic name of the additive standard deviation
#' @author Matthew Fidler
#' @noRd
.rxGetVarianceForErrorAdd <- function(env, pred1) {
  if (!is.na(pred1$a)) {
    return(.enQuote(pred1$a))
  }
  .cnd <- pred1$cond
  .w <- which(env$iniDf$err %in% c("add", "lnorm") & env$iniDf$condition == .cnd)
  if (length(.w) == 1L) {
    .p1 <- .enQuote(env$iniDf$name[.w])
  } else {
    stop("cannot find additive standard deviation", call.=FALSE)
  }
  bquote((.(.p1)) ^ 2)
}

#' Based on current error get the F that is used for prop or pow expressions
#'
#' @param env Environment of the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return The f expression
#' @author Matthew Fidler
#' @noRd
.rxGetVarianceForErrorPropOrPowF <- function(env, pred1) {
  switch(as.character(pred1$errTypeF),
         untransformed=quote(rx_pred_f_),
         transformed=quote(rx_pred_),
         f=.enQuote(pred1$f),
         none=quote(rx_pred_f_))
}

#' Get Variance for proportional error
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return The quoted proportional error
#' @author Matthew Fidler
#' @noRd
.rxGetVarianceForErrorProp <- function(env, pred1) {
  .f <- .rxGetVarianceForErrorPropOrPowF(env, pred1)
  if (!is.na(pred1$b)) {
    .p1 <- .enQuote(pred1$b)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err == "prop" & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  return(bquote((.(.f))^2*(.(.p1))^2))
}

#' Get the Variance for the additive + prop error
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return The quoted additive + proportional expression
#' @author Matthew Fidler
#' @noRd
.rxGetVarianceForErrorPow <- function(env, pred1) {
  .f <- .rxGetVarianceForErrorPropOrPowF(env, pred1)
  if (!is.na(pred1$b)) {
    .p1 <- .enQuote(pred1$b)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err == "pow" & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find power standard deviation", call.=FALSE)
    }
  }
  if (!is.na(pred1$c)) {
    .p2 <- .enQuote(pred1$c)
  } else {
    .w <- which(env$iniDf$err == "pow2" & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p2 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find exponent of power expression", call.=FALSE)
    }
  }
  bquote((.(.f))^(2 * .(.p2))*(.(.p1))^2)
}

#' Get Variance for proportional error
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return The add + prop error expression
#' @author Matthew Fidler
#' @noRd
.rxGetVarianceForErrorAddProp <- function(env, pred1) {
  if (!is.na(pred1$a)) {
    .p1 <- .enQuote(pred1$a)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err %in% c("add", "lnorm") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find additive standard deviation", call.=FALSE)
    }
  }
  .f <- .rxGetVarianceForErrorPropOrPowF(env, pred1)
  if (!is.na(pred1$b)) {
    .p2 <- .enQuote(pred1$b)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err == "prop" & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p2 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  if (pred1$addProp == "default") {
    .addProp <- getOption("rxAddProp", "combined2")
  } else {
    .addProp <- pred1$addProp
  }
  if (.addProp == "combined2") {
    return(bquote((.(.p1))^2+ (.(.f))^2*(.(.p2))^2))
  } else {
    return(bquote( ( (.(.p1)) + (.(.f)) * (.(.p2)) ) ^ 2))
  }
}

#' Additive + Power
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return additive + power
#' @author Matthew Fidler
.rxGetVarianceForErrorAddPow <- function(env, pred1) {
  if (!is.na(pred1$a)) {
    .p1 <- .enQuote(pred1$a)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err %in% c("add", "lnorm") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find additive standard deviation", call.=FALSE)
    }
  }
  .f <- .rxGetVarianceForErrorPropOrPowF(env, pred1)
  if (!is.na(pred1$b)) {
    .p2 <- .enQuote(pred1$b)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err == "prop" & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p2 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  if (!is.na(pred1$c)) {
    .p3 <- .enQuote(pred1$c)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err == "prop" & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p3 <- .enQuote(env$iniDf$name[.w])
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  if (pred1$addProp == "default") {
    .addProp <- getOption("rxAddProp", "combined2")
  } else {
    .addProp <- pred1$addProp
  }
  if (.addProp == "combined2") {
    return(bquote((.(.p1))^2 + (.(.f))^(2 * .(.p3)) * (.(.p2))^2))
  } else {
    return(bquote( ( (.(.p1)) + (.(.f)) ^ (.(.p3))* (.(.p2)) ) ^ 2))
  }
}

.rxGetVarianceForErrorType <- function(env, pred1) {
  switch(as.character(pred1$errType),
         "add"=.rxGetVarianceForErrorAdd(env, pred1), # 1
         "prop"=.rxGetVarianceForErrorProp(env, pred1), # 2
         "pow"=.rxGetVarianceForErrorPow(env, pred1), # 3
         "add + prop"=.rxGetVarianceForErrorAddProp(env, pred1),# 4
         "add + pow"=.rxGetVarianceForErrorAddPow(env, pred1) # 5
         )
}

#' Handle the single error for normal or t distributions
#'
#' @param env Environment for the parsed model
#' @param pred1 The `data.frame` of the current error
#' @return A list of the lines added.  The lines will contain
#' - `rx_yj_` which is an integer that corresponds to the transformation type.
#' - `rx_lambda_` is the transformation lambda
#' - `rx_low_` The lower boundary of the transformation
#' - `rx_hi_` The upper boundary of the transformation
#' - `rx_pred_f_` The prediction function
#' - `rx_pred_` The transformed prediction function
#' - `rx_r_` The transformed variance
#' @author Matthew Fidler
#' @noRd
.handleSingleErrTypeNormOrTFoceiBase <- function(env, pred1) {
  .ret <- vector("list", 7)
  .yj <- as.double(pred1$transform) - 1
  .ret[[1]] <- bquote(rx_yj_ ~ .(.yj))
  .ret[[2]] <- bquote(rx_lambda_~.(.rxGetLambdaFromPred1AndIni(env, pred1)))
  .ret[[3]] <- bquote(rx_low_ ~ .(.rxGetLowBoundaryPred1AndIni(env, pred1)))
  .ret[[4]] <- bquote(rx_hi_ ~ .(.rxGetHiBoundaryPred1AndIni(env, pred1)))
  .ret[[5]] <- bquote(rx_pred_f_ ~ .(.rxGetPredictionF(env, pred1)))
  .ret[[6]] <- bquote(rx_pred_ ~ .(.rxGetPredictionFTransform(env, pred1, .yj)))
  .ret[[7]] <- bquote(rx_r_ ~ .(.rxGetVarianceForErrorType(env, pred1)))
  .ret
}

#' Calculate the focei base information.
#'
#' @param env Environment for the parsed model
#' @param i The error row that is being parsed
#' @return quoted R lines for the focei setup
#' @author Matthew Fidler
#' @noRd
.handleSingleErrTypeFoceiBase <- function(env, i) {
  .pred1 <- env$predDf[i, ]
  if (tmp$predDf$distribution %in% c("norm", "t")) {
    .handleSingleErrTypeNormOrTFoceiBase(env, .pred1)
  } else {
    # This is for the non-normal cases
  }
}

