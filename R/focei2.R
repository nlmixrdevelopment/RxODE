##' Control Options for FOCEi
##'
##' @param epsilon Precision of estimate for n1qn1 optimization.
##'
##' @param maxInnerIterations Number of iterations for n1qn1
##'     optimization.
##'
##' @param maxOuterIterations Maximum number of L-BFGS-B optimization
##'     for outer problem.
##' @param n1qn1nsim Number of function evaluations for n1qn1
##'     optimization.
##'
##' @param maxstepsOde Maximum number of steps for ODE solver.
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
##' @param derivEps Central/Forward difference tolerances, which is a vector
##'     of relative difference and absolute difference.  The central/forward
##'     difference step size h is calculated as:
##'
##'         h = abs(x)*derivEps[1]+derivEps[2]
##'
##' @param derivMethod indicates the method for calculating
##'     derivatives of the outer problem.  Currently supports
##'     "central" and "forward" difference methods.
##'
##' @param covDerivMethod indicates the method for calculating the
##'     derivatives while calculating the covariance components
##'     (Hessian and S).
##'
##' @param covMethod Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of each individual's
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##' \itemize{
##' \item \code{r,s} Uses the sandwich matrix to calculate the covariance, that is: \code{R^-1 * S * R^-1}
##'
##' \item \code{r} Uses the Hessian matrix to calculate the
##'      covariance as \code{2*R^-1}
##' \item \code{s} Uses the crossproduct matrix to calculate the covariance as \code{4*S^-1}
##' }
##'
##' @param lbfgsLmm An integer giving the number of BFGS updates
##'     retained in the "L-BFGS-B" method, It defaults to 40.
##'
##' @param lbfgsPgtol is a double precision variable.
##'
##'     On entry pgtol >= 0 is specified by the user.  The iteration
##'     will stop when:
##'
##'        \code{max{|proj g_i | i = 1, ..., n} <= lbfgsPgtol}
##'
##'     where pg_i is the ith component of the projected gradient.
##'
##'     On exit pgtol is unchanged.  This defaults to zero, when the
##'     check is suppressed.
##'
##' @param lbfgsFactr Controls the convergence of the "L-BFGS-B"
##'     method.  Convergence occurs when the reduction in the
##'     objective is within this factor of the machine
##'     tolerance. Default is 1e10, which gives a tolerance of about
##'     \code{2e-6}, approximately 4 sigdigs.  You can check your
##'     exact tolerance by multiplying this value by
##'     \code{.Machine$double.eps}
##'
##' @param eigen A boolean indicating if eigenvectors are calculated
##'     to include a condition number calculation.
##'
##' @param addPosthoc Boolean indicating if posthoc parameters are
##'     added to the table output.
##'
##' @param diagXform This is the transformation used on the diagonal
##'     of the \code{chol(inv(omega))}. This matrix and values are the
##'     parameters estimated in FOCEi. The possibilities are:
##'
##' \itemize{
##' \item sqrt Estimates the sqrt of the diagonal elements of \code{chol(inv(omega))}.  This is the default method.
##' \item log Estimates the log of the diagonal elements of \code{chol(inv(omega))}
##' \item identity Estimates the diagonal elements without any transformations
##' }
##'
##' c("sqrt", "log", "identity"),
##'
##' @param sumProd Is a boolean indicating if the model should change
##'     multiplication to high precision multiplication and sums to
##'     high precision sums using the PreciseSums package.  By default
##'     this is \code{FALSE}.
##'
##'
##' @param ...
##'
##' @inheritParams rxSolve
##' @inheritParams rxSymPySetupPred
##'
##' @details
##'
##' Note this uses the R's L-BFGS-B in \code{\link{optim}} for the
##' outer problem and the BFGS \code{\link[n1qn1]{n1qn1}} with that
##' allows restoring the prior individual Hessian (for faster
##' optimization speed)
##'
##' @author Matthew L. Fidler
##'
##' @seealso \code{\link{optim}}
##' @seealso \code{\link[n1qn1]{n1qn1}}
##' @export
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
                         covDerivMethod=c("central", "forward"),
                         covMethod=c("r,s", "r", "s"),
                         lbfgsLmm=40L,
                         lbfgsPgtol=0,
                         lbfgsFactr=1e9,
                         eigen=TRUE,
                         addPosthoc=TRUE,
                         diagXform=c("sqrt", "log", "identity"),
                         sumProd=FALSE,
                         optExpression=TRUE,
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
    covDerivMethod <- .methodIdx[match.arg(covDerivMethod)];
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
    .covMethodIdx <- c("r,s" = 0L, "r"=1L, "s"=2L);
    covMethod <- .covMethodIdx[match.arg(covMethod)];
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
                 covDerivMethod=covDerivMethod,
                 covMethod=covMethod,
                 eigen=as.integer(eigen),
                 addPosthoc=as.integer(addPosthoc),
                 diagXform=match.arg(diagXform),
                 sumProd=sumProd,
                 optExpression=optExpression,
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


##' FOCEi fit
##'
##' @param data Data to fit
##' @param inits Initialization list
##' @param PKpars Pk Parameters
##' @param model The RxODE model to use
##' @param pred The Prediction function
##' @param err The Error function
##' @param lower Lower bounds
##' @param upper Upper Bounds
##' @param fixed Boolean vector indicating what parameters should be fixed.
##' @param skipCov Boolean vector indicating what parameters should be fixed when calculating covariances
##' @param control FOCEi options Control list
##' @param thetaNames Names of the thetas to be used in the final object
##' @param etaNames Eta names to be used in the final object
##' @param ... Ignored parameters
##' @return A focei fit object
##' @author Matthew L. Fidler and Wenping Wang
##' @return FOCEi fit object
##' @export
##' @author Matthew L. Fidler & Wenping Wang
foceiFit <- function(data,
                     inits,
                     PKpars,
                     model=NULL,
                     pred=NULL,
                     err=NULL,
                     lower= -Inf,
                     upper= Inf,
                     fixed = NULL,
                     skipCov=NULL,
                     control=foceiControl(),
                     thetaNames=NULL,
                     etaNames=NULL, ...){
    loadNamespace("n1qn1");
    if (!rxIs(control, "foceiControl")){
        control <- do.call(foceiControl, control);
    }
    .ret <- new.env(parent=emptyenv())
    .ret$origData <- data;
    .ret$etaNames <- etaNames;
    .ret$lower <- lower;
    .ret$upper <- upper;
    .ret$thetaFixed <- fixed;
    .ret$skipCov <- fixed;
    .ret$control <- control;
    if(is(model, "RxODE") || is(model, "character")) {
        .ret$ODEmodel <- TRUE
        if (class(pred) != "function"){
            stop("pred must be a function specifying the prediction variables in this model.")
        }
    } else {
        ## Fixme
        .ret$ODEmodel <- TRUE
        model <- constructLinCmt(PKpars);
        pred <- eval(parse(text="function(){return(Central);}"))
    }
    .square <- function(x) x*x
    .ret$diagXformInv = c("sqrt"=".square", "log"="exp", "identity"="identity")[control$diagXformInv]
    if (is.null(err)){
        err <-eval(parse(text=paste0("function(){err",paste(inits$ERROR[[1]],collapse=""),"}")));
    }
    .ret$model <- RxODE::rxSymPySetupPred(model, pred, PKpars, err, grad=(control$derivMethod == 2L),
                                          pred.minus.dv=TRUE, sum.prod=control$sumProd,
                                          theta.derivs=FALSE, optExpression=control$optExpression);

    .covNames <- .parNames <- RxODE::rxParams(.ret$model$pred.only);
    .covNames <- .covNames[regexpr(rex::rex(start, or("THETA", "ETA"), "[", numbers, "]", end), .covNames) == -1];
    colnames(data) <- sapply(names(data), function(x){
        if (any(x == .covNames)){
            return(x)
        } else {
            return(toupper(x))
        }
    })
    .lhs <- c(names(RxODE::rxInits(.ret$model$pred.only)), RxODE::rxLhs(.ret$model$pred.only))
    if (length(.lhs) > 0){
        .covNames <- .covNames[regexpr(rex::rex(start, or(.lhs), end), .covNames) == -1];
    }
    if (length(.covNames) > 0){
        if (!all(.covNames %in% names(data))){
            message("Model:")
            RxODE::rxCat(model$pred.only)
            message("Needed Covariates:")
            nlmixrPrint(.covNames)
            stop("Not all the covariates are in the dataset.")
        }
        message("Needed Covariates:")
        print(.covNames);
    }
    if (is.null(.ret$model$extra.pars)){
        .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)))
    } else {
        .nms <- c(sprintf("THETA[%s]", seq_along(inits$THTA)),
                  sprintf("ERR[%s]", seq_along(.ret$model$extra.pars)))
    }
    if (!is.null(thetaNames) && (length(inits$THTA) + length(.ret$model$extra.pars)) == length(thetaNames)){
        .nms <- thetaNames;
    }
    .ret$thetaNames <- .nms;
    if (length(lower) == 1){
        lower <- rep(lower, length(inits$THTA));
    } else if (length(lower) != length(inits$THTA)){
        print(inits$THTA)
        print(lower)
        stop("Lower must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
    }
    if (length(upper) == 1){
        upper <- rep(upper, length(inits$THTA));
    } else if (length(lower) != length(inits$THTA)){
        stop("Upper must be a single constant for all the THETA lower bounds, or match the dimension of THETA.")
    }

    .extraPars <- c();
    if (!is.null(.ret$model$extra.pars)){
        .ret$model$extra.pars <- eval(call(control$diagXform, .ret$model$extra.pars))
        if (length(.ret$model$extra.pars) > 0){
            inits$THTA <- c(inits$THTA, .ret$model$extra.pars);
            .lowerErr <- rep(control$atol * 10, length(.ret$model$extra.pars));
            .upperErr <- rep(Inf, length(.ret$model$extra.pars));
            lower <-c(lower, .lowerErr);
            upper <- c(upper, .upperErr);
        }
    }

    if (is.null(data$ID)) stop('"ID" not found in data')
    if (is.null(data$DV)) stop('"DV" not found in data')
    if (is.null(data$EVID)) data$EVID = 0
    if (is.null(data$AMT)) data$AMT = 0
    ## Make sure they are all double amounts.
    for (v in c("TIME", "AMT", "DV", cov.names))
        data[[v]] <- as.double(data[[v]]);
    .ret$dataSav = data;
    .ds <- data[data$EVID > 0, c("ID", "TIME", "AMT", .covNames)]
    data <- data[data$EVID == 0, c("ID", "TIME", "DV", .covNames)]
    ## keep the covariate names the same as in the model
    .w <- which(!(names(data.sav) %in% cov.names))
    names(.ret$dataSav)[w] <- tolower(names(.ret$dataSav[w]))         #needed in ev

    .lh = parseOM(inits$OMGA)
    .nlh = sapply(.lh, length)
    .osplt = rep(1:length(.lh), .nlh)
    .lini = list(inits$THTA, unlist(.lh));
    .nlini = sapply(.lini, length)
    .nsplt = rep(1:length(.lini), .nlini)

    .om0 = genOM(.lh)
    .ret$rxInv <- RxODE::rxSymInvCholCreate(mat=.om0, diag.xform=control$diagXform);

    .ret$thetaIni <- inits$THTA

    if (any(.ret$thetaIni == 0)){
        warning("Some of the initial conditions were 0, changing to 0.0001");
        .ret$thetaIni[.ret$thetaIni == 0] <- 0.0001;
    }
    return(.ret)
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
