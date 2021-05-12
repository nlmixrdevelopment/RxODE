## "Cl = exp(t.cl + eta.cl + cov.wt*wt + cov.wt2/70*wt)"
## "F = expit(t.f + eta.f  + cov.f*wt  + cov.w2/0.5*wt)"
## "F = probitInv()"

## .rxGetMuRefSumFirst(quote(t.cl+eta.cl+eta.cl2+eta.cl3),theta="t.cl", eta="eta.cl")
## .rxGetMuRefSumFirst(quote(t.cl + eta.cl + cov.wt*wt + cov.wt2/70*wt), theta=c("t.cl", "cov.wt", "cov.wt2"), "eta.cl")
## .rxGetMuRefSumFirst(quote(t.cl + eta.cl + cov.wt*wt + cov.wt2/70*wt), theta=c("t.cl", "cov.wt", "cov.wt2"), "eta.cl")

##' Determine if expression x is "clean"
##'
##' Clean means that it is free from calculated variables (lhs), and state.
##'
##' The expression needs to be clean of these calculated variables and
##' states to safely extract into a mu expression at the top of the
##' model.
##'
##' @param x The expression that is evaluated; must include an `info`
##'   that has a list of defined states and lhs values
##' @param env Environment for assigning information
##' @return A boolean indicating if the environment is clean from
##'   confounding elements
##' @author Matthew Fidler
##' @examples
##'
##' env <- new.env(parent=emptyenv())
##'
##' env$info <- list(state= c("depot", "center"), lhs=c("ka", "cl", "v", "cp"))
##'
##' .rxMuRefIsClean(quote(exp(tka + eta.ka)), env)
##' .rxMuRefIsClean(quote(exp(tka + eta.ka + depot)), env)
##' .rxMuRefIsClean(quote(exp(tka + eta.ka + cp)), env)
##'
##' @noRd
.rxMuRefIsClean <- function(x, env) {
  if (is.name(x) ) {
    .n <- as.character(x)
    if (any(.n == env$info$state)) {
      return(FALSE)
    } else if (any(.n == env$info$lhs)) {
      return(FALSE)
    }
    return(TRUE)
  } else if (is.call(x)) {
    return(all(unlist(lapply(x[-1], .rxMuRefIsClean, env=env))))
  } else {
    return(TRUE)
  }
}
##' Does this have an eta parameter or covariate in the expression?
##'
##' @inheritParams .rxMuRefIsClean
##' @return boolean saying if the expression has an eta parameter or
##'   covariate
##' @author Matthew Fidler
##' @examples
##'
##' env <- new.env(parent=emptyenv())
##'
##' env$info <- list(theta="tka", eta="eta.ka", cov="wt")
##'
##' .rxMuRefHasThetaEtaOrCov(quote(exp(tka + eta.ka)/v), env)
##' .rxMuRefHasThetaEtaOrCov(quote(cp/v), env)
##'
##' @noRd
.rxMuRefHasThetaEtaOrCov <- function(x, env) {
  if (is.name(x)) {
    .n <- as.character(x)
    if (any(.n == env$info$eta)) {
      return(TRUE)
    } else if (any(.n == env$info$theta)) {
      return(TRUE)
    } else if (any(.n == env$info$cov)) {
      return(TRUE)
    }
    return(FALSE)
  } else if (is.call(x)) {
    return(any(unlist(lapply(x[-1], .rxMuRefHasThetaEtaOrCov, env=env))))
  } else {
    return(FALSE)
  }
}

##' This determines if a whole line is "clean"
##'
##' This not only determines if the lhs of the line is independent of
##' any prior declarations (by `.rxMuRefIsClean()`), but it makes sure
##' that the rhs is not a special RxODE expression like `d/dt(depot)`
##' or `rate(depot)` `depot(0)` etc.
##'
##' This also assigns the line status in the enviroment `env` to
##' `.curLineClean` and also resets the current evaluation function to
##' `.curEval`
##'
##' @inheritParams .rxMuRefIsClean
##' @return boolean indicating if the line is clean
##' @author Matthew Fidler
##' @noRd
.rxMuRefLineIsClean <- function(x, env) {
  # First figure out if the mu reference line is clean
  .clean <- FALSE
  if (length(x[[2]]) == 1L && is.name(x[[2]])) {
    env$info$lhs <- c(as.character(x[[2]]), env$info$lhs)
    .clean <- TRUE
  }
  if (.clean) .clean <- .rxMuRefIsClean(x[[3]], env)
  assign(".curLineClean", .clean, env)
  assign(".curEval", "", env)
  return(env$.curLineClean)
}

##' Does this line have an eta?
##'
##' @inheritParams .rxMuRefIsClean
##' @return boolean indicating the line has etas
##' @author Matthew Fidler
##' @noRd
.rxMuRefLineHasEta <- function(x, env) {
  if (is.name(x)) {
    .n <- as.character(x)
    if (any(.n == env$info$eta)) {
      return(TRUE)
    }
    return(FALSE)
  } else if (is.call(x)) {
    return(any(unlist(lapply(x[-1], .rxMuRefLineHasEta, env=env))))
  } else {
    return(FALSE)
  }
}

.rxIsOp <- function(x) {
  (identical(x[[1]], quote(`*`)) ||
     identical(x[[1]], quote(`^`)) ||
     identical(x[[1]], quote(`+`)) ||
     identical(x[[1]], quote(`-`)) ||
     identical(x[[1]], quote(`/`)))
}

.rxIsLogicalOp <- function(x) {
  (identical(x[[1]], quote(`==`)) ||
     identical(x[[1]], quote(`>`)) ||
     identical(x[[1]], quote(`<`)) ||
     identical(x[[1]], quote(`<=`)) ||
     identical(x[[1]], quote(`>=`)) ||
     identical(x[[1]], quote(`!=`)) ||
     identical(x[[1]], quote(`&&`)) ||
     identical(x[[1]], quote(`||`)) ||
     identical(x[[1]], quote(`|`)) ||
     identical(x[[1]], quote(`&`)))
}
##' Handle the limit for logit/expit types of functions
##'
##' @param x This is the parsed expression tree (from R) that has limits
##' @param env environment where parsing information is saved
##' @return Nothing; Called for its side effects
##' @author Matthew Fidler
##' @noRd
.rxMuRefHandleLimits <- function(x, env) {
  if (length(x) == 4) {
    # expit(x, 1, 2)
    if (is.numeric(x[[3]])) {
      env$curLow <- as.numeric(x[[3]])
    } else {
      env$err <- unique(c(env$err, paste0("syntax error '", deparse1(x), "': limits must be numeric")))
      env$curLow <- -Inf
    }
    if (is.numeric(x[[4]])) {
      env$curHi <- as.numeric(x[[4]])
    } else {
      env$err <- unique(c(env$err, paste0("syntax error '", deparse1(x), "': limits must be numeric")))
      env$curHi <- Inf
    }
    x <- x[1:2]
  } else if (length(x) == 3) {
    # expit(x, 1)
    if (is.numeric(x[[3]])) {
      env$curLow <- as.numeric(x[[3]])
    } else {
      env$err <- unique(c(env$err, paste0("syntax error '", deparse1(x), "': limits must be numeric")))
      env$curLow <- -Inf
    }
    env$curHi <- 1
    x <- x[1:2]
  } else {
    env$curLow <- 0
    env$curHi <- 1
  }
  if (env$curLow >= env$curHi) {
    env$err <- unique(c(env$err, paste0("syntax error '", deparse1(x), "': limits must be lower, higher")))
  }
  x
}

##' Extract single variable names from a expression
##'
##' @param x Expression
##' @param names Names to prepend to the final names output
##' @param env Environment that contains a boolean `.found` that
##'   indicates if an interesting expression has been found
##' @return character vector of names that are uncomplicated, like
##'   "a", "b"; If the names are part of a larger expression these
##'   names are skipped.
##' @author Matthew Fidler
##' @noRd
.muRefExtractSingleVariableNames <- function(x, names, env) {
  c(names, do.call(`c`, lapply(x, function(y) {
    if(is.name(y)) {
      env$found <- TRUE
      return(as.character(y))
    }
    return(NULL)
  })))
}
##' Extract mu-style covariates that is theta + eta + cov*theta.cov
##'
##' @param x expression to evaluate
##'
##' @param doubleNames A list of the covariates with the estimates
##'   attached.  This is in a single expression so wt*theta.cov1 +
##'   wt*theta.cov2 will add an error at the end of the expression
##'
##' @param env Environment with $info that has information about the
##'   parsed model
##'
##' @return A list of covariates with estimates attached
##'
##' @author Matthew Fidler
##' @noRd
.muRefExtractMultiplyMuCovariates <- function(x, doubleNames, env) {
  c(doubleNames, do.call(`c`, lapply(x, function(y) {
    if(is.call(y) && identical(y[[1]], quote(`*`))) {
      .y2 <- y[-1]
      if (length(.y2) == 2) {
        if (is.name(.y2[[1]]) && is.name(.y2[[2]])) {
          .y2 <- vapply(.y2, as.character, character(1))
          if (any(.y2[1] == env$info$cov) &&
                any(.y2[2] == env$info$theta)) {
            if (any(.y2[1] == names(doubleNames))) {
              env$err <- unique(c(env$err, paste0("syntax error: covariate '", .y2[1],
                                                  "' is duplicated in mu-referenced expression for '",
                                                  .y2[2], "' and '", doubleNames[[.y2[1]]], "'")))
            }
            env$.found <- TRUE
            return(setNames(list(.y2[2]), .y2[1]))
          } else if (any(.y2[2] == env$info$cov) &&
                       any(.y2[1] == env$info$theta)) {
            if (any(.y2[2] == names(doubleNames))) {
              env$err <- unique(c(env$err, paste0("syntax error: covariate '", .y2[2],
                                                  "' is duplicated in mu-referenced expression for '",
                                                  .y2[1], "' and '", doubleNames[[.y2[2]]], "'")))
            }
            env$.found <- TRUE
            return(setNames(list(.y2[1]), .y2[2]))
          }
        }
      }
      return(NULL)
    }
    return(NULL)
  })))
}

.muRefNextAdditiveExpression <- function(x) {
  .expr <- NULL
  for (i in seq_along(length(x))) {
    if(is.call(x[[i]])) {
      .expr <- x[[i]]
      if (identical(.expr[[1]], quote(`+`))){
        .expr <- .expr[-1]
      } else {
        .expr <- NULL
      }
    }
  }
  .expr
}

# This function handles the extra information in a theta based mu referenced

##' This handles the case where there is a mu referenced single theta with a reference population eta
##'
##' @param .we this represents the item number in `.names` that is the
##'   mu-referenced population eta
##'
##' @param .wt this represents the item number in `.names` that is the
##'   mu-referenced population theta
##'
##' @param .names This is the names of the single item names in the
##'   additive expression that may be a mu-referenced expression
##'
##' @param .doubleNames The double names of the covariate and
##'   estimated theta for the mu referenced covariate.
##'
##' The double names is a list. The elements of this list are the
##' population parameter estimates for a data-based covariate. The
##' names of the list are the covariate that the population parameters
##' are estimating an effect for.
##'
##' @param .extraItems This represents the extra items that are in the
##'   additive expression that may be part of the full additive
##'   expression
##'
##' @param env This is the mu-reference environment where mu reference
##'   information is saved
##'
##' @return Nothing; Called for its side effects and modifications in
##'   the mu reference environment `env`
##'
##' @author Matthew Fidler
##' @noRd
.muRefHandleSingleThetaMuRef <- function(.we, .wt, .names, .doubleNames, .extraItems, env) {
  # Here the mu reference is possible
  #print(.names)
  #message("double names")
  #print(.doubleNames)
  if (length(.we) == 1) {
    # Simple theta/eta mu-referencing
    .curEta <- .names[.we]
    .wEtaInDf <- which(env$nonMuEtas$eta == .curEta)
    if (length(.wEtaInDf) > 0) {
      # This is when the mu referenced eta has already been observed at least once in a prior expression
      #
      # Here we will check to see if this should be downgraded to an additive expression
      if (!all(env$nonMuEtas$curEval[.wEtaInDf] == env$.curEval)) {
        # In this case, there is a specialized expression in more than
        # one form; for example:
        #
        # q1 = exp(tq+eta.q)
        # q2 = expit(tq+eta.q)
        #
        # The mu expression would be an additive expression tq+eta.q
        #
        # Downgrade to additive expression
        env$nonMuEtas$curEval[.wEtaInDf] <- ""
      }
    } else {
      .wEtaInDf <- which(env$muRefDataFrame$eta == .curEta)
      .ce <- env$.curEval
      if (length(.wEtaInDf) > 0) {
        # duplicated ETAs, if everything is not the same then it isn't really mu-referenced

        # The duplicated ETAs can occur for shared etas for 2
        # different population parameters.  This can make sense for 2
        # different treatments with a similar between subject
        # variability that can be pooled, like:
        #
        # emaxA = tv.emaxA + eta.emax
        # emaxB = tv.emaxB + eta.emax
        #
        # In this case, the eta.emax is no longer mu-referenced
        #
        if (!all(env$muRefDataFrame$theta[.wEtaInDf] == .names[.wt]) |
              !all(env$muRefDataFrame$eta[.wEtaInDf] == .curEta)) {
          if (!all(env$muRefDataFrame$curEval[.wEtaInDf] == .ce)) {
            .ce <- ""
          }
          env$nonMuEtas <- rbind(env$nonMuEtas, data.frame(eta=.curEta, curEval=.ce))
          ## Here we check to see
          env$muRefDataFrame <- env$muRefDataFrame[-.wEtaInDf,, drop = FALSE]
        } else {
          if (!all(env$muRefDataFrame$curEval[.wEtaInDf] == env$.curEval)) {
            # Downgrade to additive expression
            env$muRefDataFrame$curEval[.wEtaInDf] <- ""
          }
        }
      } else {
        env$muRefDataFrame <- rbind(env$muRefDataFrame, data.frame(theta=.names[.wt], eta=.names[.we], curEval=env$.curEval))
        if (!is.null(.extraItems)) {
          env$muRefThetaExtra <- rbind(env$muRefThetaExtra, data.frame(theta=.names[.wt], extra=.extraItems))
        }
      }
    }
  } else if (length(.we) != 0) {
    stop("currently do not support IOV etc")
  }
  if (length(.wt) == 1) {
    if (length(.doubleNames) > 0) {
      .doubleNames <- .doubleNames[names(.doubleNames) != ""]
      if (length(.doubleNames) > 0) {
        .w <- which(env$muRefCovariateDataFrame$popParameter == .names[.wt])
        .covariate <- names(.doubleNames)
        .covariateParameter <- setNames(unlist(.doubleNames), NULL)
        if (length(.w) > 0) {
          # Already defined something
          .tmp <- env$muRefCovariateDataFrame[.w, ]
          .tmp2 <- with(.tmp, paste0(covariate, ",", covariateParameter))
          if (!all(.tmp2 %in% paste0(.covariate, ",", .covariateParameter))) {
            stop(sprintf("improper covariate mu-referencing for '%s', more than one parameter is being estimated for a single covariate", .names[.wt]))
          }
        }
        .df <- data.frame(popParameter=.names[.wt], covariate=.covariate, covariateParameter=.covariateParameter)
        env$muRefCovariateDataFrame <- rbind(env$muRefCovariateDataFrame, .df)
      }
    }
  }
}


## To reduce code from nlmixr, the
## Covariate references should have the following structure:
## To reduce code from nlmixr the mu reference:
## -------------------------------
## > f$nmodel$mu.ref
## $eta.ka
## [1] "tka"
##
## $eta.cl
## [1] "tcl"
##
## $eta.v
## [1] "tv"

## > f$cov.ref
## $age
## cov.age
##   "tcl"
##
## $wt
## cov.wt
##  "tcl"


##' Handle the + expressions to determine mu-reference expressions
##'
##' @param x additive Call Expression
##' @param env Environment information
##' @return nothing
##' @author Matthew Fidler
##' @noRd
.muRefHandlePlus <- function(x, env) {
  .x2 <- x[-1]
  .names <- NULL
  .doubleNames <- list(0)
  .extraItems <- NULL
  while (!is.null(.x2)) {
    env$.found <- FALSE
    .names <- .muRefExtractSingleVariableNames(.x2, .names, env)
    .doubleNames <- .muRefExtractMultiplyMuCovariates(.x2, .doubleNames, env)
    if (!env$.found && !is.name(.x2[[2]])) {
      .extraItems <- c(.extraItems, deparse1(.x2[[2]]))
    }
    .x2 <- .muRefNextAdditiveExpression(.x2)
  }
  .wt <- which(.names %in% env$info$theta)
  .we <- which(.names %in% env$info$eta)
  if (length(.wt) >= 2) {
    env$err <- unique(c(env$err,
                        paste0("syntax error: 2+ single population parameters in a single mu-referenced expression: '",
                               paste(env$info$theta[.wt], collapse="', '"), "'")))
  } else if (length(.wt) == 1) {
    #print(.extraItems)
    .ord <- order(vapply(.extraItems, nchar, integer(0)))
    .muRefHandleSingleThetaMuRef(.we, .wt, .names, .doubleNames, .extraItems, env)
  } else if (length(.we) == 1){
    .curEta <- .names[.we]
    .wEtaInDf <- which(env$nonMuEtas$eta == .curEta)
    .ce <- env$.curEval
    if (length(.wEtaInDf) > 0) {
      if (!all(env$muRefDataFrame$curEval[.wEtaInDf] == env$.curEval)) {
        # Downgrade to additive expression
        env$muRefDataFrame$curEval[.wEtaInDf] <- ""
      }
    } else {
      env$nonMuEtas <- rbind(env$nonMuEtas, data.frame(eta=.curEta, curEval=.ce))
    }
  }
  invisible()
}


.handleSingleEtaIfExists <- function(expr, env) {
  .curEta <- deparse1(expr)
  if (any(.curEta == env$info$eta)) {
    ## eta singlet
    .wEtaInDf <- which(env$muRefDataFrame$eta == .curEta)
    .ce <- env$.curEval
    if (length(.wEtaInDf) > 0) {
      # duplicated ETAs, if everything is not the same then it isn't really mu-referenced
      if (!all(env$muRefDataFrame$curEval[.wEtaInDf] == .ce)) {
        .ce <- ""
      }
      env$nonMuEtas <- rbind(env$nonMuEtas, data.frame(eta=.curEta, curEval=.ce))
      env$muRefDataFrame <- env$muRefDataFrame[-.wEtaInDf,, drop = FALSE]
    } else {
      .wEtaInDf <- which(env$nonMuEtas$eta == .curEta)
      .ce <- env$.curEval
      if (length(.wEtaInDf) > 0) {
        if (!all(env$nonMuEtas$curEval[.wEtaInDf] == env$.curEval)) {
          # Downgrade to additive expression
          env$nonMuEtas$curEval[.wEtaInDf] <- ""
        }
      } else {
        env$nonMuEtas <- rbind(env$nonMuEtas, data.frame(eta=.curEta, curEval=.ce))
      }
    }
  }
}


## f$probit.theta.low  f$probit.theta.hi
## f$probit.theta
## f$logit.theta
## f$log.theta
## "tka"     "tcl"     "cov.wt"  "cov.age" "tv"

## f$probit.eta.low  f$probit.eta.hi
## f$probit.eta
## f$probit.eta.low  f$probit.eta.hi
## f$logit.eta
## f$log.eta
## "eta.ka" "eta.cl" "eta.v"
##
## > f$oneTheta
## "tka" "tcl" "tv"
.rxMuRef0 <- function(x, env) {
  if (is.call(x)) {
    if (env$top && identical(x[[1]], quote(`{`))) {
      env$top <- FALSE
      y <- x[-1]
      for (.i in seq_along(y)) {
        assign(".curEval", "", env)
        x <- y[[.i]]
        if (identical(x[[1]], quote(`=`)) ||
              identical(x[[1]], quote(`~`))) {
          .handleSingleEtaIfExists(x[[3]], env)
          if (.rxMuRefHasThetaEtaOrCov(x[[3]], env)){
            # This line has etas or covariates and might need to be
            # separated into mu-referenced line
            .rxMuRefLineIsClean(x, env)
            lapply(x, .rxMuRef0, env=env)
          } else {
            # This line does not depend on etas or covariates
            # simply add to the body
            env$body <- c(env$body, list(x))
          }
        } else {
          ## This line is a special statement, simply add to the body
          env$body <- c(env$body, list(x))
        }
      }
    } else if (identical(x[[1]], quote(`+`))) {
      .muRefHandlePlus(x, env)
    } else {
      assign(".curEval", as.character(x[[1]]), env)
      .handleSingleEtaIfExists(x[[2]], env)
      if (env$.curEval == "probitInv" ||
            env$.curEval == "expit" ||
            env$.curEval == "logit" ||
            env$.curEval == "probit") {
        x <- .rxMuRefHandleLimits(x, env)
      }
      lapply(x[-1], .rxMuRef0, env=env)
    }
  } else if (is.name(x)) {
    .name <- as.character(x)
    if (any(.name == env$info$theta)) {
      env$oneTheta <- unique(.name, env$oneTheta)
    }
  }
}

## 1. $state : states
## 2. $params : params
## 3. $lhs: lhs
## 4. theta: theta from ini
## 5. eta: eta from ini

.rxMuRefSetupInitialEnvironment <- function(mod, ini) {
  .eta <- dimnames(ini)[[1]]
  .iniDf <- as.data.frame(ini)
  .theta <- .iniDf$name[!is.na(.iniDf$ntheta)]
  .mv  <- rxModelVars(mod)
  .expr <- eval(parse(text=paste0("quote({",rxNorm(.mv),"})")))
  .state <- .mv$state

  .params <- .mv$params
  .lhs <- .mv$lhs
  # Covariates are model based parameters not described by theta/eta
  .info <- list(state=.state,
                lhs=NULL,
                theta=.theta,
                eta=.eta,
                cov=setdiff(.params, c(.theta, .eta)))
  .env <- new.env(parent=emptyenv())
  .env$param <- list()
  .env$body <- list()
  .env$info <- .info
  .env$top <- TRUE

  # probit/probitInv
  .env$probit.theta.low <- NULL
  .env$probit.theta.hi <- NULL
  .env$probit.theta <- NULL

  .env$probitInv.theta.low <- NULL
  .env$probitInv.theta.hi <- NULL
  .env$probitInv.theta <- NULL

  # logit/expit
  .env$logit.theta <- NULL
  .env$logit.theta.low <- NULL
  .env$logit.theta.hi <- NULL

  .env$expit.theta <- NULL
  .env$expit.theta.low <- NULL
  .env$expit.theta.hi <- NULL

  .env$log.theta <- NULL
  .env$exp.theta <- NULL

  .env$cov.ref <- NULL
  .env$err <- NULL
  .env$.expr <- .expr
  .env$muRefDataFrame <- data.frame(eta=character(0), theta=character(0), curEval=character(0))
  .env$muRefThetaExtra <- data.frame(theta=character(0), extra=character(0))
  .env$muRefCovariateDataFrame <- data.frame(eta=character(0), theta=character(0), curEval=character(0))
  .env$nonMuEtas <- data.frame(popParameter=character(0), covariate=character(0), covariateParameter=character(0))
  return(.env)
}

##' Get mu-referencing model from model variables
##'
##' The rxMuRef is the core of the nlmixr ui functions
##'
##' This function takes the initialization values from `lotri()` the
##' parsed RxODE model to generate mu-referenced models adding
##' mu-references for etas that do not have them to allow saem to
##' support non mu-referenced models by a parsing trick.
##'
##' @param mod Model
##' @param ini lotri initialization information
##' @return
##' @author Matthew Fidler
##' @examples
##'
##' # First lets get a lotri initialization block:
##'
##' ini <- lotri({
##'    ## You may label each parameter with a comment
##'    tka <- 0.45 # Log Ka
##'    tcl <- log(c(0, 2.7, 100)) # Log Cl
##'    ## This works with interactive models
##'    ## You may also label the preceding line with label("label text")
##'    tv <- 3.45; label("log V")
##'    ## the label("Label name") works with all models
##'    eta.ka ~ 0.6
##'    eta.cl + eta.v ~ c(0.3,
##'                       0.001, 0.1)
##'    add.sd <- 0.7
##' })
##'
##'
##' ini2 <- lotri({
##'        ## You may label each parameter with a comment
##'        tka <- 0.45 # Log Ka
##'        tcl <- log(c(0, 2.7, 100)) # Log Cl
##'        ## This works with interactive models
##'        ## You may also label the preceding line with label("label text")
##'        tv <- 3.45; label("log V")
##'        add.sd <- 0.7
##'     })
##'
##' @export
rxMuRef <- function(mod, ini=NULL) {
  if (!inherits(ini, "lotriFix")) {
    stop("requires a lotri object with at least one fixed effect", call.=FALSE)
  }
  .env <- .rxMuRefSetupInitialEnvironment(mod, ini)
  .rxMuRef0(.env$.expr, env=.env)
  if (length(.env$err) > 0) {
    stop(paste0("syntax/parsing errors:\n",
                paste(.env$err, collapse="\n")),
         call.=FALSE)
  }
  return(invisible(.env))
}
