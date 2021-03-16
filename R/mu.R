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
##' `curLineClean` and also resets the current evaluation function to
##' `curEval`
##'
##' @inheritParams .rxMuRefIsClean
##' @return boolean indicating if the line is clean
##' @author Matthew Fidler
##' @noRd
.rxMuRefLineIsClean <- function(x, env) {
  # First figure out if the
  .clean <- FALSE
  if (length(x[[2]]) == 1L && is.name(x[[2]])){
    env$info$lhs <- c(as.character(x[[2]]), env$info$lhs)
    .clean <- TRUE
  }
  if (.clean) .clean <- .rxMuRefIsClean(x[[3]], env)
  assign("curLineClean", .clean, env)
  assign("curEval", "", env)
  return(env$curLineClean)
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

.rxMuRefHandleLimits <- function(x, env) {
  if (length(x) == 4) {
    # expit(x, 1, 2)
    print(x[[1]])
    print(x[[2]])
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

## To reduce code from nlmixr, the
## Covariate references should have the following structure:
## ------------------------------
## > f$cov.ref
## $age
## cov.age
##   "tcl"
##
## $wt
## cov.wt
##  "tcl"


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
  force(env)
  if (is.call(x)) {
    if (env$top && identical(x[[1]], quote(`{`))) {
      env$top <- FALSE
      y <- x[-1]
      for (.i in seq_along(y)) {
        x <- y[[.i]]
        if (identical(x[[1]], quote(`=`)) ||
              identical(x[[1]], quote(`~`))) {
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
      stop()
    } else if (.rxIsOp(x)) {
      print(x)
    } else if (!.rxIsLogicalOp(x)){
      assign("curEval", as.character(x[[1]]), env)
      if (env$curEval == "probitInv" ||
            env$curEval == "expit") {
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
rxMuRef <- function(mod, theta=NULL, eta=NULL) {
 .mv  <- rxModelVars(mod)
 .expr <- eval(parse(text=paste0("quote({",rxNorm(.mv),"})")))
 .state <- .mv$state
 .params <- .mv$params
 .lhs <- .mv$lhs
 # Covariates are model based parameters not described by theta/eta
 .info <- list(state=.state,
               lhs=NULL,
               theta=theta,
               eta=eta,
               cov=setdiff(.params, c(theta, eta)))
 .env <- new.env(parent=emptyenv())
 .env$param <- list()
 .env$body <- list()
 .env$info <- .info
 .env$top <- TRUE
 .env$probit.theta.low <- NULL
 .env$probit.theta.hi <- NULL
 .env$probit.theta <- NULL
 .env$logit.theta <- NULL
 .env$log.theta <- NULL
 .env$cov.ref <- NULL
 .env$err <- NULL
 .rxMuRef0(.expr, env=.env)
 if (length(.env$err) > 0) {
   stop(paste0("syntax/parsing errors:",
               paste(.env$err, collapse="\n")),
        call.=FALSE)
 }
 return(invisible())
}
