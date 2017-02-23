rxFoceiEtaSetup <- function(object, ..., dv, eta, omegaInv, nonmem=FALSE, log.det.OMGAinv.5=0){
    args <- list(object=object, ..., eta=eta);
    setup <- do.call(getFromNamespace("rxSolveSetup", "RxODE"), args, envir = parent.frame(1));
    return(with(setup, {
        tmp <- environment(object$.call);
        if (any(names(tmp) == "eta.trans")){
            eta.trans <- tmp$eta.trans;
        } else {
            n <- names(params)
            max.eta <- n[regexpr(regEta, n) != -1]
            max.eta <- max(as.numeric(gsub(regEta, "\\1", max.eta)));
            eta.trans <- sapply(seq(1, max.eta), function(x){return(which(n == sprintf("ETA[%s]", x)) - 1)});
            tmp$eta.trans <- eta.trans;
        }
        setup$nonmem <- as.integer(nonmem)
        setup$DV <- as.double(dv);
        setup$omegaInv <- omegaInv;
        setup$eta <- NULL;
        setup$eta.trans <- as.integer(eta.trans);
        setup$object <- object;
        setup$eta.mat <- matrix(eta, ncol=1);
        setup$log.det.OMGAinv.5 <- 0;
        return(list2env(setup));
    }))
}

##' FOCEI ETA setup,
##'
##' This is basically for testing
##'
##' @param object RxODE object
##' @param ... values sent to rxFoceiEtaSetup
##' @param dv dependent variable
##' @param eta eta values.
##' @param omegaInv Inverse omega
##' @param env Environment to use, instead of rxFoceiEtaSetup environment
##' @param nonmem Match NONMEMs approximation of the hessian.
##' @return environment of solved information.
##' @author Matthew L. Fidler
##' @keywords internal
rxFoceiEta <- function(object, ..., dv, eta, omegaInv, env, nonmem=FALSE){
    UseMethod("rxFoceiEta");
}
##' @rdname rxFoceiEta
rxFoceiEta.RxODE <- function(object, ..., dv, eta, omegaInv, env, nonmem){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiEta.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiEta
rxFoceiEta.rxDll <- function(object, ..., dv, eta, omegaInv, env, nonmem){
    if (missing(env)){
        args <- as.list(match.call(expand.dots=TRUE))[-1];
        env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    }
    return(object$.call(rxTrans(object)["ode_solver_focei_eta"], eta, env));
}

##' Get -LLIK for individual
##'
##' Used for testing.
##'
##' @inheritParams rxFoceiEta
##' @return -llik
##' @author Matthew L. Fidler
##' @keywords internal
rxFoceiLik <- function(object, ..., dv, eta, omegaInv){
    UseMethod("rxFoceiLik");
}
##' @rdname rxFoceiLik
rxFoceiLik.RxODE <- function(object, ..., dv, eta, omegaInv){
    return(rxFoceiLik.rxDll(object=object$cmpMgr$rxDll(), ..., dv=dv, eta=eta, omegaInv=omegaInv));
}
##' @rdname rxFoceiLik
rxFoceiLik.rxDll <- function(object, ..., dv, eta, omegaInv){
    object$.call(rxTrans(object)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    return(RxODE_focei_eta_lik(eta, env))
}

##' Get -Grad for individual
##'
##' For testing
##' @inheritParams rxFoceiEta
##' @return -grad
##' @author Matthew L. Fidler
##' @keywords internal
rxFoceiLp <- function(object, ..., dv, eta, omegaInv){
    UseMethod("rxFoceiLp");
}
##' @rdname rxFoceiLp
rxFoceiLp.RxODE <- function(object, ..., dv, eta, omegaInv){
    return(rxFoceiLp.rxDll(object=object$cmpMgr$rxDll(), ..., dv=dv, eta=eta, omegaInv=omegaInv));
}
##' @rdname rxFoceiLp
rxFoceiLp.rxDll <- function(object, ..., dv, eta, omegaInv){
    object$.call(rxTrans(object)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    return(RxODE_focei_eta_lp(eta, env))
}
##' Solve the FOCEI inner problem
##'
##' @param object RxODE object
##' @param ... values sent to rxFoceiEtaSetup
##' @param dv dependent variable
##' @param eta eta value
##' @param omegaInv inverse omega matrix
##' @param invisible Invisible LBFGS tracing.
##' @param m LBFGS m
##' @param epsilon LBFGS epsilon
##' @param past LBFGS past
##' @param delta LBFGS delta
##' @param max_iterations LBFGS max_iterations
##' @param linesearch_algorithm LBFGS linesearch_algorithm
##' @param max_linesearch LBFGS max_linesearch
##' @param min_step LBFGS min_step
##' @param max_step LBFGS max_step
##' @param ftol LBFGS ftol
##' @param wolfe LBFGS wolfe
##' @param gtol LBFGS gtol
##' @param orthantwise_c LBFGS orthantwise_c
##' @param orthantwise_start LBFGS orthantwise_start
##' @param orthantwise_end LBFGS orthantwise_end
##' @return Loglik with fitted and posthoc attributes
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiInner <- function(object, ..., dv, eta, omegaInv,
                         invisible = 0, m = 6, epsilon = 1e-05, past = 0, delta = 0,
                         max_iterations = 0, linesearch_algorithm = "LBFGS_LINESEARCH_DEFAULT",
                         max_linesearch = 20, min_step = 1e-20, max_step = 1e+20,
                         ftol = 1e-04, wolfe = 0.9, gtol = 0.9, orthantwise_c = 0,
                         orthantwise_start = 0, orthantwise_end = length(eta)){
    UseMethod("rxFoceiInner");
}
##' @rdname rxFoceiInner
##' @export
rxFoceiInner.RxODE <- function(object, ..., dv, eta, omegaInv,
                               invisible = 0, m = 6, epsilon = 1e-05, past = 0, delta = 0,
                               max_iterations = 0, linesearch_algorithm = "LBFGS_LINESEARCH_DEFAULT",
                               max_linesearch = 20, min_step = 1e-20, max_step = 1e+20,
                               ftol = 1e-04, wolfe = 0.9, gtol = 0.9, orthantwise_c = 0,
                               orthantwise_start = 0, orthantwise_end = length(eta)){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiInner.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiInner
##' @export
rxFoceiInner.rxDll <- function(object, ..., dv, eta, omegaInv,
                               invisible = 0, m = 6, epsilon = 1e-05, past = 0, delta = 0,
                               max_iterations = 0, linesearch_algorithm = "LBFGS_LINESEARCH_DEFAULT",
                               max_linesearch = 20, min_step = 1e-20, max_step = 1e+20,
                               ftol = 1e-04, wolfe = 0.9, gtol = 0.9, orthantwise_c = 0,
                               orthantwise_start = 0, orthantwise_end = length(eta)){
    object$.call(rxTrans(object)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    lik <- RxODE_focei_eta("lik");
    lp <- RxODE_focei_eta("lp")
    output <- lbfgs::lbfgs(lik, lp, eta, environment=env,
                           invisible=invisible, m=m, epsilon=epsilon, past=past, delta=delta,
                           max_iterations=max_iterations, linesearch_algorithm=linesearch_algorithm,
                           max_linesearch=max_linesearch, min_step=min_step, max_step=max_step,
                           ftol=ftol, wolfe=wolfe, gtol=gtol, orthantwise_c=orthantwise_c,
                           orthantwise_start=orthantwise_start, orthantwise_end = orthantwise_end);
    return(RxODE_focei_finalize_llik(env));
}
