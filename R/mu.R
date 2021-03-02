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
.rxGetMuRefSumFirst <- function(expr, theta=c(), eta=c()) {
  ## Looking for eta + theta + theta*cov  ## This function is the initial function called for any sum
  .env <- new.env(parent=emptyenv())
  .env$.thetaSingle <- NULL
  .env$.etaSingle <- NULL
  .env$.covThetaMu <- NULL
  .env$.badThetaMu
  .rxGetMuRefSumRecur(expr, theta, eta, .env)
  ## print(.env$.thetaSingle)
  ## print(.env$.etaSingle)
  ## print(.env$.covThetaMu)
}
