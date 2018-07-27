##' Control Options for FOCEi
##'
##' @param epsilon Precision of estimate for n1qn1 optimization.
##'
##' @param maxInnerIterations Number of iterations for n1qn1
##'     optimization.
##'
##' @param n1qn1nsim Number of function evaluations for n1qn1
##'     optimization.
##'
##' @param maxstepsOde Maximum number of steps for ODE solver.
##'
##' @param tbs Type of transformation used on dependent variable and
##'     prediction for dynamic transform both sides method.
##'
##' @param printInner Integer representing when the inner step is
##'     printed. By default this is 0 or do not print.  1 is print
##'     every function evaluation, 5 is print every 5 evaluations.
##'
##' @param printOuter Integer representing when the outer step is
##'     printed. When this is 0 or do not print the iterations.  1 is
##'     print every function evaluation (default), 5 is print every 5
##'     evaluations.
##'
##' @param scaleTo Scale the initial parameter estimate to this value.
##'     By default this is 1.
##'
##' @param derivMethod indicates the method for calculating
##'     derivatives of the outer problem.  Currently supports
##'     "central" and "forward" difference methods.
##'
##' @param derivEps Central/Forward difference tolerances, which is a vector
##'     of relative difference and absolute difference.  The central/forward
##'     difference step size h is calculated as:
##'
##'         h = abs(x)*derivEps[1]+derivEps[2]
##'
##' @param lbfgsFactr Controls the convergence of the "L-BFGS-B"
##'     method.  Convergence occurs when the reduction in the
##'     objective is within this factor of the machine
##'     tolerance. Default is 1e10, which gives a tolerance of about
##'     \code{2e-6}, approximately 4 sigdigs.  You can check your
##'     exact tolerance by multiplying this value by
##'     \code{.Machine$double.eps}
##'
##'
##' @inheritParams rxSolve
##'
##' @export
##'
##' @author Matthew L. Fidler
foceiControl <- function(epsilon=.Machine$double.eps,
                         maxInnerIterations=10000,
                         maxOuterIterations=50000,
                         n1qn1nsim=NULL,
                         method = c("liblsoda", "lsoda", "dop853"),
                         transitAbs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
                         maxstepsOde = 5000L, hmin = 0L, hmax = NULL, hini = 0, maxordn = 12L, maxords = 5L, cores,
                         covsInterpolation = c("linear", "locf", "nocb", "midpoint"),
                         printInner=0L,
                         printOuter=1L,
                         scaleTo=1.0,
                         derivEps=c(1.0e-6, 1.0e-6),
                         derivMethod=c("forward", "central"),
                         lbfgsLmm=5L,
                         lbfgsPgtol=0,
                         lbfgsFactr=1e9,
                         ## outerOpt=c("lbfgsb", "qnbd"),
                         ..., stiff){
    .xtra <- list(...);
    if (is.null(transitAbs) && !is.null(.xtra$transit_abs)){  # nolint
        transitAbs <- .xtra$transit_abs;  # nolint
    }
    if (missing(covsInterpolation) && !is.null(.xtra$covs_interpolation)){  # nolint
        covsInterpolation <- .xtra$covs_interpolation; # nolint
    }
    if (missing(maxInnerIterations) && !is.null(.xtra$max_iterations)){  # nolint
        maxInnerIterations <- .xtra$max_iterations; # nolint
    }
    if (!missing(stiff) && missing(method)){
        if (rxIs(stiff, "logical")){
            if (stiff){
                method <- "lsoda"
                warning("stiff=TRUE has been replaced with method = \"lsoda\".")
            } else {
                method <- "dop853"
                warning("stiff=FALSE has been replaced with method = \"dop853\".")
            }
        }
    }  else {
        method <- match.arg(method);
    }
    .methodIdx <- c("lsoda"=1L, "dop853"=0L, "liblsoda"=2L);
    method <- as.integer(.methodIdx[method]);
    derivMethod <- match.arg(derivMethod);
    .methodIdx <- c("forward"=0L, "central"=1L);
    derivMethod <- as.integer(.methodIdx[derivMethod])
    if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1];
    covsInterpolation <- tolower(match.arg(covsInterpolation,
                                           c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")))
    if (covsInterpolation == "constant") covsInterpolation <- "locf";
    covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
    if (missing(cores)){
        cores <- RxODE::rxCores();
    }
    if (missing(n1qn1nsim)){
        n1qn1nsim <- 10 * maxInnerIterations + 1;
    }
    ## outerOpt <- match.arg(outerOpt)
    ## .outerIdx <- c("lbfgsb"=0L, "qnbd"=1L)
    ## outerOpt <- as.integer(.outerIdx[outerOpt]);
    .ret <- list(maxOuterIterations=as.integer(maxOuterIterations),
                 maxInnerIterations=as.integer(maxInnerIterations),
                 method=method,
                 transitAbs=transitAbs,
                 atol=atol,
                 rtol=rtol,
                 maxstepsOde=maxstepsOde,
                 hmin=hmin,
                 hmax=hmax,
                 hini=hini,
                 maxordn=maxordn,
                 maxords=maxords,
                 cores=cores,
                 covsInterpolation=covsInterpolation,
                 n1qn1nsim=as.integer(n1qn1nsim),
                 printInner=as.integer(printInner),
                 printOuter=as.integer(printOuter),
                 lbfgsLmm=as.integer(lbfgsLmm),
                 lbfgsPgtol=as.double(lbfgsPgtol),
                 lbfgsFactr=as.double(lbfgsFactr),
                 scaleTo=scaleTo,
                 epsilon=epsilon,
                 derivEps=derivEps,
                 derivMethod=derivMethod,
                 outerOpt=0L)
    class(.ret) <- "foceiControl"
    return(.ret);
}

.foceiSetup <- function(obj, data, theta, thetaFixed = NULL,
                        skipCov=NULL, rxInv = NULL,
                        lower = NULL, upper = NULL, etaMat = NULL,
                        control = RxODE::foceiControl()){
    loadNamespace("n1qn1");
    .Call(`_RxODE_foceiSetup_`, obj, data, theta, thetaFixed, skipCov, rxInv, lower, upper, etaMat, control); # nolint
}


.nearPd <- function(mat){
    if (any(is.na(mat))){
        cat("Bad matrix:\n");
        print(mat);
        return(mat)
    } else {
        return(as.matrix(Matrix::nearPD(mat)$mat));
    }
}
