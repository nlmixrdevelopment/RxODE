## "Cl = exp(t.cl + eta.cl + cov.wt*wt + cov.wt2/70*wt)"
## "F = expit(t.f + eta.f  + cov.f*wt  + cov.w2/0.5*wt)"
## "F = probitInv()"


.rxGetMuRefTheta <- NULL
.rxGetMuRefEta <- NULL
rxGetMuRef <- function(model, theta, eta){
  assignInMyNamespace(".rxGetMuRefTheta", theta)
  assignInMyNamespace(".rxGetMuRefEta", eta)
}

.rxGetMuRefSumRecur <- function(expr, theta, eta, env) {
  if (is.name(expr)) {
    .n <- as.character(expr)
    if (any(.n == theta)) {
      assign(".thetaSingle", unique(c(env$.thetaSingle, .n)), envir=env)
    } else if (any(.n == eta)) {
      assign(".etaSingle", unique(c(env$.etaSingle, .n)), envir=env)
    }
  } else if (identical(expr[[1]], quote(`+`)) ||
               identical(expr[[1]], quote(`-`))) {
    lapply(expr[-1], .rxGetMuRefSumRecur, theta=theta, eta=eta, env=env)
  } else if (length(expr) == 3) {
    if (identical(expr[[1]], quote(`*`)) && is.name(expr[[2]]) && is.name(expr[[3]])) {
      .expr <- vapply(expr[-1], as.character, character(1))
      .idx <- .expr %in% theta
      if (sum(.idx) == 1L) {
        # This expression is theta=mu
        .theta <- .expr[.idx]
        .cov <- .expr[!.idx]
        .w <- which(names(env$.covThetaMu) == .theta)
        if (length(.w) == 1) {
          ## In this case
          ## t.cl + eta.cl + t.wt*wt + t.wt*wt2
          ## Perhaps warn, but at least this should not be a mu-referenced covariate
          assign(".covThetaMu", .covThetaMu[-.w], envir=env)
          assign(".badThetaMu", unique(c(.badThetaMu, .theta)), envir=env)
        } else if (!any(.theta == env$.badThetaMu)) { # Make sure this isn't a bad mu-referenced theta
          assign(".covThetaMu", c(stats::setNames(list(.cov), .theta),
                                  env$.covThetaMu), envir=env)
        }
      }
    }
  }
}

## .rxGetMuRefSumFirst(quote(t.cl+eta.cl+eta.cl2+eta.cl3),theta="t.cl", eta="eta.cl")
## .rxGetMuRefSumFirst(quote(t.cl + eta.cl + cov.wt*wt + cov.wt2/70*wt), theta=c("t.cl", "cov.wt", "cov.wt2"), "eta.cl")
## .rxGetMuRefSumFirst(quote(t.cl + eta.cl + cov.wt*wt + cov.wt2/70*wt), theta=c("t.cl", "cov.wt", "cov.wt2"), "eta.cl")
.rxGetMuRefSumFirst <- function(expr, theta=c(), eta=c()) {
  ## Looking for eta + theta + theta*cov  ## This function is the initial function called for any sum
  .env <- new.env(parent=emptyenv())
  .env$.thetaSingle <- NULL
  .env$.etaSingle <- NULL
  .env$.covThetaMu <- NULL
  .env$.badThetaMu <- NULL
  .rxGetMuRefSumRecur(expr, theta, eta, .env)
  ## print(.env$.thetaSingle)
  ## print(.env$.etaSingle)
  ## print(.env$.covThetaMu)

  ## Case 1: exp()
  ## exp(t.cl + eta.cl + theta*cov + extra)
  ## exp(t.cl + eta.cl + theta*cov)*exp(extra)
  ## mu1 = exp(t.cl + eta.cl + theta*cov)

  ## Case 2/3: expit() Can't be as easily decomposed as exp/add
}

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
  if (is.name(x)) {
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

.rxMuRef0 <- function(x, env) {
  if (is.call(x)) {
    if (identical(x[[1]], quote(`=`)) ||
          identical(x[[1]], quote(`~`))) {
      assign("curLine", x, env)
      .clean <- FALSE
      if (length(x[[2]]) == 1L && is.name(x[[2]])){
        env$info$lhs <- c(as.character(x[[2]]), env$info$lhs)
        .clean <- TRUE
      }
      if (.clean) .clean <- .rxMuRefIsClean(x[[3]], env)
      assign("curLineClean", .clean, env)
      assign("curEval", "", env)
    } else if (identical(x[[1]], quote(`+`))) {

    } else {
      assign("curEval", as.character(x[[1]]), env)
    }
    lapply(x[-1], .rxMuRef0, env=env)
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
 .env$info <- .info
 .rxMuRef0(.expr, .env)
 return(invisible())
}
