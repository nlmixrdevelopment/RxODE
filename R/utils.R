##' Unloads all RxODE compiled DLLs
##'
##' @return List of RxODE dlls still loaded
##'
##' @examples
##'
##' print(rxUnloadAll())
##'
##' @export
rxUnloadAll <- function() {
  .Call(`_RxODE_rxUnloadAll_`) # nolint
  .dll <- getLoadedDLLs()
  .n <- names(.dll)
  .n <- .n[regexpr("^rx_[a-f0-9]{32}", .n) != -1]
  .ret <- lapply(.n, function(rx) {
    .cur <- .dll[[rx]]
    class(.cur) <- NULL
    .num <- .rxModels[[.cur$path]]
    if (is.null(.num)) {
      dyn.unload(.cur$path)
      return(NULL)
    } else if (identical(.num, 0L)) {
      dyn.unload(.cur$path)
      return(NULL)
    } else {
      return(.cur$path)
    }
  })
  if (length(.ret) == 0) return(invisible(.ret))
  .ret <- .ret[[!is.null(.ret)]]
  return(invisible(.ret))
}

.normalizePath <- function(path, ...) {
  ifelse(.Platform$OS.type == "windows",
         suppressWarnings(utils::shortPathName(normalizePath(path, ...))),
  ifelse(regexpr("^[/~]", path) != -1,
         suppressWarnings(normalizePath(path, ...)),
         suppressWarnings(normalizePath(file.path(getwd(), path), ...))))
}

##' Require namespace, otherwise throw error.
##'
##' @param pkg Package required for function to work.
##' @return Nothing
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxReq <- function(pkg) {
  ## nocov start
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("package \"%s\" needed for this function to work", pkg),
         call. = FALSE);
  }
  ## nocov end
}
##' Use cat when RxODE.verbose is TRUE
##'
##' @param ... Parameters sent to cat
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxCat <- function(a, ...) {
  ## nocov start
  if (RxODE.verbose) {
    if (is(a, "RxODE")) {
      message(RxODE::rxNorm(a), appendLF = FALSE);
    } else {
      message(a, ..., appendLF = FALSE);
    }
  }
  ## nocov end
}



##' Cleanup anonymous DLLs by unloading them
##'
##' This cleans up any RxODE loaded DLLs
##'
##' @param wd What directory should be cleaned; (DEPRECIATED), this no
##'     longer does anything.
##'
##' This unloads all RxODE anonymous dlls.
##'
##' @return TRUE if successful
##'
##' @author Matthew L. Fidler
##' @export
rxClean <- function(wd) {
  if (!missing(wd)) warning("'wd' is depreciated")
  rxUnloadAll();
}

refresh <- function(derivs=FALSE) {
  ## nocov start
  cat("Dparser Version\n");
  print(dparser::dpVersion())
  if (derivs){
    Sys.setenv(RxODE_derivs=TRUE)
  } else {
    Sys.setenv(RxODE_derivs=FALSE)
  }
  source(devtools::package_file("build/refresh.R"))
  ## nocov end
}

ode.h <- function() {
  ## nocov start
  cat("Generate header string.\n");
  r.files <- list.files(devtools::package_file("R"), "[.]R$", full.names=TRUE);
  r.files <- r.files[regexpr("RxODE_md5.R", r.files, fixed=TRUE) == -1]
  md5 <- digest::digest(c(sapply(list.files(devtools::package_file("src"),
                                            pattern="[.](c|cpp|h|hpp|f|R|in)$", full.names=TRUE),
                                 function(x){digest::digest(x, file=TRUE)}),
                          sapply(r.files, function(x){digest::digest(x, file=TRUE)}),
                          sapply(list.files(devtools::package_file("vignettes"), pattern="[.](Rmd)$",
                                            full.names=TRUE),
                                 function(x){digest::digest(x, file=TRUE)}),
                          sapply(list.files(devtools::package_file("demo"), pattern="(R|index)$",
                                            full.names=TRUE),
                                 function(x){digest::digest(x, file=TRUE)}),
                          sapply(list.files(devtools::package_file(),
                                            pattern="(cleanup.*|configure.*|DESCRIPTION|NAMESPACE)",
                                            full.names=TRUE),
                                 function(x){digest::digest(x, file=TRUE)})))
  hd <- c(sprintf("#define __VER_2__ \"    SET_STRING_ELT(version,2,mkChar(\\\"%s\\\"));\\n\"", md5),
          sprintf("#define __VER_1__ \"    SET_STRING_ELT(version,1,mkChar(\\\"%s\\\"));\\n\"",
                  as.vector(RxODE::rxVersion()["repo"])),
          sprintf("#define __VER_0__ \"    SET_STRING_ELT(version,0,mkChar(\\\"%s\\\"));\\n\"",
                  sessionInfo()$otherPkgs$RxODE$Version),
          sprintf("#define __VER_md5__ \"%s\"", md5),
          sprintf("#define __VER_repo__ \"%s\"", as.vector(RxODE::rxVersion()["repo"])),
          sprintf("#define __VER_ver__ \"%s\"", sessionInfo()$otherPkgs$RxODE$Version))
  writeLines(hd, devtools::package_file("src/ode.h"))
  writeLines(sprintf("RxODE.md5 <- \"%s\"", md5),
             devtools::package_file("R/RxODE_md5.R"));
  ## nocov end
}




##' Defunct setting of sum
##'
##' @param type used to be type of product
##'
##' @export
rxSetSum <- function(type=c("pairwise", "fsum", "kahan", "neumaier", "c")) {
  error("'rxSetSum' has been moved to rxSolve(...,sum=)");
}

##' Defunct setting of product
##'
##' @param type used to be type of product
##'
##' @export
rxSetProd <- function(type=c("long double", "double", "logify")) {
  error("'rxSetProd' has been moved to rxSolve(...,sum=)");
}

##' Set timing for progress bar
##'
##' @param seconds This sets the number of seconds that need to elapse
##'     before drawing the next segment of the progress bar.  When
##'     this is zero or below this turns off the progress bar.
##'
##' @export
##' @author Matthew Fidler
rxSetProgressBar <- function(seconds=1.0) {
  invisible(.Call(`_rxParProgress`, as.double(seconds)))
}

##' Sample a covariance Matrix from the Posterior Inverse Wishart
##' distribution.
##'
##' Note this Inverse wishart rescaled to match the original scale of
##' the covariance matrix.
##'
##' If your covariance matrix is a 1x1 matrix, this uses an scaled
##' inverse chi-squared which is equivalent to the Inverse Wishart
##' distribution in the uni-directional case.
##'
##' @param nu Degrees of Freedom (Number of Observations) for
##'        covariance matrix simulation.
##'
##' @param omega Either the estimate of covariance matrix or the
##'     estimated standard deviations in matrix form each row forming
##'     the standard deviation simulated values
##'
##' @param n Number of Matrices to sample.  By default this is 1.
##'     This is only useful when \code{omega} is a matrix.  Otherwise
##'     it is determined by the number of rows in the input
##'     \code{omega} matrix of standard deviations
##'
##' @param omegaIsChol is an indicator of if the omega matrix is in
##'   the Cholesky decomposition. This is only used when code{type="invWishart"}
##'
##' @param returnChol Return the Cholesky decomposition of the
##'   covariance matrix sample. This is only used when code{type="invWishart"}
##'
##' @param diagXformType Diagonal transformation type.  These could be:
##'
##' \itemize{
##'
##' \item{log} The standard deviations are log transformed, so the
##'   actual standard deviations are exp(omega)
##'
##' \item{identity} The standard deviations are not transformed. The
##' standard deviations are not transformed;  They should be positive.
##'
##' \item{variance} The variances are specified in the \code{omega}
##' matrix; They are transformed into standard deviations.
##'
##' \item{nlmixrSqrt} These standard deviations come from an nlmixr
##' omega matrix where diag(chol(inv(omega))) = x^2
##'
##' \item{nlmixrLog} These standard deviations come from a nlmixr
##' omega matrix omega matrix where diag(chol(solve(omega))) = exp(x)
##'
##' \item{nlmixrIdentity} These standard deviations come from a nlmixr
##' omega matrix omega matrix where diag(chol(solve(omega))) = x
##'
##' }
##'
##'  The nlmixr transformations only make sense when there is no
##'  off-diagonal correlations modeled.
##'
##' @param type The type of covariance posterior that is being
##'     simulated.  This can be:
##'
##' \itemize{
##'
##' \item{invWishart} The posterior is an inverse wishart; This allows
##' for correlations between parameters to be modeled.  All the
##' uncertainty in the parameter is captured in the degrees of freedom
##' parameter.
##'
##' \item{lkj} The posterior separates the standard deviation
##' estimates (modeled outside and provided in the \code{omega}
##' argument) and the correlation estimates. The correlation estimate
##' is simulated with the \code{\link{rLKJ1}}.  This simulation uses
##' the relationship \code{eta=(nu-1)/2}.  This is relationship based
##' on the proof of the relationship between the restricted
##' LKJ-distribution and inverse wishart distribution (XXXXXX).  Once
##' the correlation posterior is calculated, the estimated standard
##' deviations are then combined with the simulated correlation matrix
##' to create the covariance matrix.
##'
##' \item{separation} Like the \code{lkj} option, this separates out
##' the estimation of the correlation and standard deviation.  Instead
##' of using the \code{LKJ} distribution to simulate the correlation,
##' it simulates the inverse wishart of the identity matrix and
##' converts the result to a correlation matrix.  This correlation
##' matrix is then used with the standard deviation to calculate the
##' simulated covariance matrix.
##'
##' }
##'
##' @return a matrix (n=1) or a list of matrices  (n > 1)
##'
##' @details
##'
##' In general, the separation strategy is preferred for diagonal
##' matrices.  If the dimension of the matrix is below 10, \code{lkj}
##' is numerically faster than \code{separation} method.  However, the
##' \code{lkj} method has densities too close to zero (XXXX) when the
##' dimension is above 10.  In that case, though computationally more
##' expensive \code{separation} method performs better.
##'
##' For matrices with modeled covariances, the easiest method to use
##' is the inverse Wishart which allows the simulation of correlation
##' matrices (XXXX).  This method is more well suited for well behaved
##' matrices, that is the variance components are not too low or too
##' high.  When modeling nonlinear mixed effects modeling matrices
##' with too high or low variances are considered sub-optimal in
##' describing a system.  With these rules in mind, it is reasonable
##' to use the inverse Wishart.
##'
##' @author Matthew L.Fidler & Wenping Wang
##'
##' @references
##'
##' Alvarez I, Niemi J and Simpson M. (2014) \emph{Bayesian Inference for a
##' Covariance Matrix}. Conference on Applied Statistics in Agriculture.
##' \url{https://newprairiepress.org/cgi/viewcontent.cgi?article=1004&context=agstatconference}
##'
##'
##' Wang1 Z, Wu Y, and Chu H. (2018) \emph{On Equivalence of the LKJ
##' distribution and the restricted Wishart distribution}. arXiv:1809.04746
##'
##' @examples
##'
##' ## Sample a single covariance.
##' draw1 <- cvPost(3, matrix(c(1,.3,.3,1),2,2))
##'
##' ## Sample 3 covariances
##' set.seed(42)
##' draw3 <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3)
##'
##' ## Sample 3 covariances, but return the cholesky decomposition
##' set.seed(42)
##' draw3c <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3, returnChol=TRUE)
##'
##' ## Sample 3 covariances with lognormal standard deviations via LKJ
##' ## correlation sample
##' cvPost(3,sapply(1:3,function(...){rnorm(10)}), type="lkj")
##'
##' ## or return cholesky decomposition
##' cvPost(3,sapply(1:3,function(...){rnorm(10)}), type="lkj",
##'   returnChol=TRUE)
##'
##' ## Sample 3 covariances with lognormal standard deviations via separation
##' ## strategy using inverse Wishart correlation sample
##' cvPost(3,sapply(1:3,function(...){rnorm(10)}), type="separation")
##'
##' ## or returning the cholesky decomposition
##' cvPost(3,sapply(1:3,function(...){rnorm(10)}), type="separation",
##'   returnChol=TRUE)
##'
##' @export
cvPost <- function(nu, omega, n = 1L, omegaIsChol = FALSE, returnChol = FALSE,
                   type = c("invWishart", "lkj", "separation"),
                   diagXformType = c("log", "identity", "variance", "nlmixrSqrt", "nlmixrLog", "nlmixrIdentity")) {
  if (is.null(nu) && n == 1L) return(omega)
  if (inherits(type, "numeric") || inherits(type, "integer")) {
    .type <- as.integer(type)
  } else {
    .type <- as.vector(c("invWishart"=1L, "lkj"=2L,
                         "separation"=3L)[match.arg(type)])
  }
  if (.type == 1L) {
    .xform <- 1L
  }  else if (inherits(diagXformType, "numeric") || inherits(diagXformType, "integer")) {
    .xform <- as.integer(diagXformType)
  } else {
    .xform <- setNames(c("variance" = 6L, "log" = 5L,
                         "identity" = 4L, "nlmixrSqrt" = 1L,
                         "nlmixrLog" = 2L,
                         "nlmixrIdentity" = 3L)[match.arg(omegaXform)],
                       NULL)
  }
  .ret <- .Call(`_RxODE_cvPost_`, nu, omega, n,
                omegaIsChol, returnChol, .type, .xform,
                PACKAGE = "RxODE")
  return(.ret)
}

##' Simulate from a (truncated) multivariate normal
##'
##' This is simulated with the fast, thread-safe threefry simulator
##' and can use multiple cores to generate the random deviates.
##'
##' @param n Number of random row vectors to be simulated OR the
##'     matrix to use for simulation (faster).
##'
##' @param mu mean vector
##'
##' @param sigma Covariance matrix for multivariate normal or a list
##'   of covariance matrices. If a list of covariance matrix, each
##'   matrix will simulate `n` matrices and combine them to a full
##'   matrix
##'
##' @param lower is a vector of the lower bound for the truncated
##'     multivariate norm
##'
##' @param upper is a vector of the upper bound for the truncated
##'     multivariate norm
##'
##' @param ncores Number of cores used in the simulation
##'
##' @param isChol A boolean indicating if \code{sigma} is a cholesky
##'     decomposition of the covariance matrix.
##'
##' @param keepNames Keep the names from either the mean or covariance
##'     matrix.
##'
##' @return
##'
##' If \code{n==integer} (default) the output is an (n x d) matrix
##' where the i-th row is the i-th simulated vector.
##'
##' If \code{is.matrix(n)} then the random vector are store in \code{n},
##' which is provided by the user, and the function returns
##' \code{NULL} invisibly.
##'
##' @references John K. Salmon, Mark A. Moraes, Ron O. Dror, and David
##'     E. Shaw (2011). Parallel Random Numbers: As Easy as 1, 2, 3.
##'     D. E. Shaw Research, New York, NY 10036, USA.
##'
##' @examples
##'
##' ## From mvnfast
##' ## Unlike mvnfast, uses threefry simulation
##'
##'  d <- 5
##' mu <- 1:d
##'
##' # Creating covariance matrix
##' tmp <- matrix(rnorm(d^2), d, d)
##' mcov <- tcrossprod(tmp, tmp)
##'
##'
##' set.seed(414)
##' rxRmvn(4, 1:d, mcov)
##'
##' set.seed(414)
##' rxRmvn(4, 1:d, mcov)
##'
##' set.seed(414)
##' rxRmvn(4, 1:d, mcov, ncores = 2) # r.v. generated on the second core are different
##'
##' ###### Here we create the matrix that will hold the simulated
##' #  random variables upfront.
##' A <- matrix(NA, 4, d)
##' class(A) <- "numeric" # This is important. We need the elements of A to be of class "numeric".
##'
##' set.seed(414)
##' rxRmvn(A, 1:d, mcov, ncores = 2) # This returns NULL ...
##' A                                # ... but the result is here
##'
##' ## You can also simulate from a truncated normal:
##'
##' rxRmvn(10, 1:d, mcov, lower=1:d-1, upper=1:d+1)
##'
##'
##' # You can also simulate from different matrices (if they match
##' # dimensions) by using a list of matrices.
##'
##' matL <- lapply(1:4,function(...){
##'    tmp <- matrix(rnorm(d^2), d, d)
##'    tcrossprod(tmp, tmp)
##' })
##' 
##'
##' rxRmvn(4, setNames(1:d,paste0("a",1:d)), matL)
##'
##' @author Matthew Fidler, Zdravko Botev and some from Matteo Fasiolo
##' @export
rxRmvn <- function(n, mu=NULL, sigma, lower= -Inf, upper=Inf, ncores=1, isChol=FALSE,
                   keepNames=TRUE, a=0.4, tol=2.05, nlTol=1e-10, nlMaxiter=100L) {
  .ret <- .Call(`_RxODE_rxRmvnSEXP`, n, mu, sigma, lower, upper, ncores,
                isChol, keepNames, a, tol, nlTol, nlMaxiter)
  if (is.matrix(n)) {
    return(invisible())
  }
  return(.ret)
}

##' Collect warnings and just warn once.
##'
##' @param expr R expression
##' @param lst When \code{TRUE} return a list with
##'     list(object,warnings) instead of issuing the warnings.
##'     Otherwise, when \code{FALSE} issue the warnings and return the
##'     object.
##' @return The value of the expression or a list with the value of
##'     the expression and a list of warning messages
##' @author Matthew L. Fidler
##' @noRd
.collectWarnings <- function(expr, lst=FALSE) {
  .ws <- c();
  .thisEnv <- environment()
  .ret <- suppressWarnings(
    withCallingHandlers(expr,
                        warning = function(w) {
                          assign(".ws", unique(c(w$message, .ws)), .thisEnv)
                        }))
  if (lst) {
    return(list(.ret, .ws))
  } else {
    for (.w in .ws) {
      warning(.w)
    }
    return(.ret)
  }
}

##' Convert numeric vector to repeated data.frame
##'
##' @param vec Named input vector
##' @param n Number of columns
##' @return Data frame with repeated vec
##' @author Matthew Fidler
##' @noRd
.vecDf <- function(vec, n) {
  .Call(`_vecDF`, vec, as.integer(n), PACKAGE='RxODE') #nolint
}
##' cbind Ome
##'
##' @param et The theta data frame
##' @param mat The full matrix simulation from omegas
##' @param n number of subject simulated
##' @return data frame with et combined with simulated omega matrix values
##' @author Matthew Fidler
##' @noRd
.cbindOme <- function(et, mat, n) {
  .Call(`_cbindOme`, et, mat, as.integer(n), PACKAGE='RxODE') # nolint
}
