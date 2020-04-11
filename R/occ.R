##' Expand the omega matrix to include the IOV matrix
##'
##' @param omega Original (named) omega matrix
##' @param iov Original (named) iov matrix
##' @param iovNames IOV names that are used in the model
##' @param omegaMatNames Omega matrix names (if not in omega matrix)
##' @param iovMatNames IOV matrix names (if not in iov matrix)
##' @return List with three elements; `m` is Named matrix for internal
##'     R simulation. `lower`/`upper` are the upper and lower
##'     boundaries of the truncated multi-variate normal.
##' @author Matthew Fidler
##' @examples
##'
##' omega <- lotri(eta.cl+eta.v ~c(1, 0.5, 1))
##' iov <- lotri(iov.cl+iov.v ~c(1, 0.5, 1))
##' .rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"))
##'
##' .rxExpandOmegaIov(omega)
##' @noRd
.rxExpandOmegaIov <- function(omega, iov=NULL, iovNames=NULL,
                              omegaLower=-Inf, omegaUpper=Inf,
                              iovLower=-Inf, iovUpper=Inf,
                              omegaMatNames=NULL, iovMatNames=NULL) {
  if (!is.null(omegaMatNames)) {
    .dim1 <- omegaMatNames
  } else {
    if (is.null(dimnames(omega))) stop("'omega' needs to be named")
    .dim1 <- dimnames(omega)[[1]]
  }
  .len0 <- length(.dim1)
  if (length(omegaLower) == 1) omegaLower <- rep(omegaLower, .len0)
  else if (length(omegaLower) != .len0)
    stop("'omega' and 'omegaLower' need compatible dimensions")
  if (length(omegaUpper) == 1) omegaUpper <- rep(omegaUpper, .len0)
  else if (length(omegaUpper) != .len0)
    stop("'omega' and 'omegaUpper' need compatible dimensions")
  if (is.null(iov) && is.null(iovNames)){
    return(list(m=omega, lower=omegaLower,
                upper=omegaUpper))
  } else {
    if (!is.null(iovMatNames)) {
      .dim2 <- iovMatNames
    } else {
      if (is.null(dimnames(iov))) stop("'iov' needs to be named")
      .dim2 <- dimnames(iov)[[1]]
    }
    if (is.null(iovNames)) stop("'iovNames' needs to have model-based terms")
    .len1 <- length(.dim2)
    if (length(iovLower) == 1) iovLower <- rep(iovLower, .len1)
    else if (length(iovLower) != .len0)
      stop("'iov' and 'iovLower' need compatible dimensions")
    if (length(iovUpper) == 1) iovUpper <- rep(iovUpper, .len1)
    else if (length(iovUpper) != .len0)
      stop("'iov' and 'iovUpper' need compatible dimensions")
    .len2 <- length(iovNames)
    if (.len2 %% .len1 != 0)
      stop('iov and full iov expansion does not match');
    .len3 <- .len2 / .len1
    iovLower <- rep(iovLower, .len3)
    iovUpper <- rep(iovUpper, .len3)
    .m <- lapply(seq(1, .len3), function(...) iov)
    .m <- as.matrix(Matrix::bdiag(c(list(omega), .m)));
    .dn <- c(.dim1, iovNames)
    dimnames(.m) <- list(.dn, .dn)
    return(list(m=.m,
                lower=c(omegaLower, iovLower),
                upper=c(omegaUpper, iovUpper)))
  }
}


##' This expands the theta parameters to a data frame with one value
##' per study
##'
##' @inheritParams rxSolve
##'
##' @noRd
.expandTheta <- function(theta, thetaMat=NULL,
                         thetaLower= -Inf, thetaUpper=Inf, nStud=1L,
                         nCoresRV=1L) {
  checkmate::assertNumeric(theta, finite=TRUE, any.missing=FALSE,
                           min.len=1, names="strict")
  if (inherits(nStud, "numeric")) nStud <- as.integer(nStud)
  checkmate::assertInteger(nStud, len=1, lower=1)
  if (inherits(nCoresRV, "numeric")) nCoresRV <- as.integer(nCoresRV)
  checkmate::assertInteger(nCoresRV, len=1, lower=1)
  
  if (is.null(thetaMat)) {
    if (inherits(theta, "numeric")){
      return(.Call(`_vecDF`, theta, as.integer(nStud), PACKAGE='RxODE'))#nolint
    } else if (inherits(theta, "matrix")) {
      return(as.data.frame(theta));
    } else if (inherits(theta, "data.frame")) {
      return(theta)
    }
  } else {
    checkmate::assertMatrix(thetaMat, mode="numeric",
                            any.missing=FALSE, col.names="strict",
                            min.cols = 1)
    .dim <- dim(thetaMat)
    if (.dim[1] != .dim[2]) {
      stop("'thetaMat' must be a symmetric, square numeric matrix")
    }
    ## if (all(thetaMat == t(thetaMat))) {
    ##   stop("'thetaMat' must be a symmetric, square numeric matrix (t)")
    ## }
    .matn <- dimnames(thetaMat)[[2]]
    .n <- names(theta)
    if (length(.matn) != length(.n)) {
      stop("dimensions of 'thetaMat' must match the number of elements in 'theta'")
    }
    .theta <- theta[.matn]
    if (any(is.na(.theta))) {
      stop("names in 'theta' do not match names in 'thetaMat' column names")
    }
    return(as.data.frame(rxRmvn(nStud, .theta, thetaMat,
                                lower=thetaLower, upper=thetaUpper,
                                ncores=nCoresRV)))
  }
}


.expandPars <- function(object, params, events, control) {
  ## The event table is needed to get nesting information from the table
  
  ## In common between nested and non-nested models is the expanded
  ## population matrix.
  .et <- .expandTheta(theta=params,
                      thetaMat=control$thetaMat,
                      thetaLower=control$thetaLower,
                      thetaUpper=control$thetaUpper,
                      nStud=control$nStud,
                      nCoresRV=control$nCoresRV)
  .omega <- control$omega;
  if (inherits(.omega, "matrix")) {
    ## Promote to lotri matrix with correct properties
    .names <- dimnames(.omega)[[2]]
    if (is.null(.names)) {
      stop("'omega' must be named for simulations")
    }
    if (control$omegaIsChol) {
      .dimnames <- dimnames(.omega)
      .omega <- .omega %*% t(.omega)
      dimnames(omega) <- .omega
      control$omegaIsChol <- FALSE
    }
    .omega <- lotri::as.lotri(.omega, lower=control$omegaLower,
                              upper=control$omegaUpper, nu=control$omegaDf,
                              default="id")
    control$omegaUpper <-NULL
    control$omegaLower <- NULL
    control$omegaDf <-NULL
  }
  if (inherits(.omega, "lotri")) {
    control$omega <- .omega
    if (control$omegaIsChol) {
      stop("'omegaIsChol' makes no sense with lotri matrices")
    }
    ## This is for hierarchical models.
    .n <- names(events)
    .nl <- tolower(.n)
    .w <- which(.nl == "id")
    if (length(.w) != 1) {
      stop("malformed 'id' column in event data when expanding nested levels and parameters")
    }
    .id <- events[, .w];
    .ni <- .nestingInfo(id=.id, omega=control$omega, data=events)
    .nid <- length(levels(.ni$id))
    if (control$nSub == 1) {
      control$nSub <- .nid
    }
    .en <- rxExpandNesting(obj=object, .ni, compile=TRUE)
    ## now we can expand omega matrices based on the "theta" values above.
    ## We need to determine if the ALL the omega matrix value names
    ## are in the expanded theta matrix.  If so we can use the ijk or
    ## separation strategy
    .allNames <- unlist(lapply(names(control$omega),
                               function(n){
                                 dimnames(control$omega[[n]])[[2]]
                               }))
    .hasDT <- requireNamespace("data.table", quietly = TRUE)
    .method <- control$omegaSeparation
    ## First get the method when .method == "auto"
    if (.method == "auto") {
      if (.hasDT) {
        .all <- all(data.table::`%chin%`(.allNames, names(.expandTheta)))
      } else {
        .all <- all(.allNames %in% names(.expandTheta))
      }
      if (.all && length(.allNames) > 9L) {
        .method <- "separation"
      } else if (.all) {
        .method <- "ijk"
      } else {
        .method <- "invWishart"
      }
    }
    ## Here we choose the type of n needed to generate the above and
    ## below omega matrices.
    if (any(.method == c("separation", "ijk"))) {
      ## In this case, the n is a matrix of the expanded theta
      .n <- as.matrix(.et)
    } else {
      .n <- control$nStud
    }
    .above <- .en$aboveF(.n, type=.method,
                        diagXformType=control$omegaXform)
    .below <- .en$belowF(.n, type=.method,
                         diagXformType=control$omegaXform)
    ## This gives the above ID variability and integrates it into the
    ## study level variability
    if (!is.null(.above)) {
      .et2 <- as.data.frame(rxRmvn(1L, sigma=.above,
                                   lower=.en$aboveUpper(),
                                   upper=.en$aboveLower(),
                                   ncores=control$nCoresRV))
      .et2 <- .et2[, !(names(.et2) %in% names(.et)), drop=FALSE]
      .et <- cbind(.et, .et2)
    }
    .rxModels[".theta"] <- .et
    ## Now simulate the number of subjects for each study
    .ind <- rxRmvn(control$nSub, sigma=.below,
                   lower=.en$belowLower(),
                   upper=.en$belowUpper(),
                   ncores=control$nCoresRV)
    if (rxIs(.below, "list")) {
      .rxModels[".omegaL"] <- .below
    }
    ## object, params, events, control
    .ni2 <- .ni[names(.ni) != "data"]
    .ni2 <- .ni2[names(.ni2) != "id"]
    .rxModels[".nestObj"] <- .en$mod
    .rxModels[".nestEvents"] <- .ni$data
    .rxModels[".nestInfo"] <- .ni2
    return(.Call(`_cbindOme`, .et, .ind,
                             as.integer(control$nSub),
                             PACKAGE='RxODE'))
  } else {
    .rxModels[".nestObj"] <- object
    .rxModels[".nestEvents"] <- events
    .rxModels[".theta"] <- .et
    .rxModels[".nestInfo"] <- NULL
    return(.et)
  }
}
