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
##' @param printInner Print information in the inner optimization
##'
##' @param scaleTo Scale the initial parameter estimate to this value.
##'     By default this is 1.
##'
##' @param centralEps Central difference tolerances, which is a vector
##'     of relative difference and absolute difference.  The central
##'     difference step size h is calculated as:
##'
##'         h = abs(x)*centralEps[1]+centralEps[2]
##'
##' @inheritParams rxSolve
##'
##' @export
##'
##' @author Matthew L. Fidler
foceiControl <- function(epsilon=.Machine$double.eps,
                         maxInnerIterations=100,
                         maxOuterIterations=5000,
                         n1qn1nsim=NULL,
                         method = c("liblsoda", "lsoda", "dop853"),
                         transitAbs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
                         maxstepsOde = 5000L, hmin = 0L, hmax = NULL, hini = 0, maxordn = 12L, maxords = 5L, cores,
                         covsInterpolation = c("linear", "locf", "nocb", "midpoint"),
                         tbs=c("cox-box", "yeo-johnson"),
                         lambda=1.0,
                         estLambda=FALSE,
                         printInner=FALSE,
                         scaleTo=1.0,
                         centralEps=c(0.5e-6, 0.5e-6),
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
    .methodIdx <- c("lsoda"=1, "dop853"=0, "liblsoda"=2);
    method <- as.integer(.methodIdx[method]);
    if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1];
    covsInterpolation <- tolower(match.arg(covsInterpolation,
                                           c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")))
    if (covsInterpolation == "constant") covsInterpolation <- "locf";
    covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
    if (missing(cores)){
        cores <- RxODE::rxCores();
    }
    tbs <- match.arg(tbs);
    .tbsNames <- c("cox-box"=0L, "yeo-johnson"=1L);
    tbs <- as.vector(.tbsNames[tbs]);
    if (missing(n1qn1nsim)){
        n1qn1nsim <- 10 * maxInnerIterations + 1;
    }
    list(maxOuterIterations=as.integer(maxOuterIterations),
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
         tbs=tbs,
         n1qn1nsim=as.integer(n1qn1nsim),
         printInner=as.integer(printInner),
         scaleTo=scaleTo,
         estLambda=as.integer(estLambda),
         lambda=lambda,
         epsilon=epsilon,
         centralEps=centralEps)
}

.foceiSetup <- function(obj, data, theta, thetaFixed = NULL, rxInv = NULL,
                        lower = NULL, upper = NULL, etaMat = NULL,
                        odeOpts = RxODE::foceiControl()){
    loadNamespace("n1qn1");
    .Call(`_RxODE_foceiSetup_`, obj, data, theta, thetaFixed, rxInv, lower, upper, etaMat, odeOpts); # nolint
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
