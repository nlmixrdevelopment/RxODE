rxFoceiEtaSetup <- function(object, ..., dv, eta, theta, nonmem=FALSE, inv.env=parent.frame(1), id= -1,
                            inits.vec=NULL){
    args <- list(object=object, ..., eta=eta, theta=theta);
    setup <- c(do.call(getFromNamespace("rxSolveSetup", "RxODE"), args, envir = parent.frame(1)),
               as.list(inv.env));
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
        setup$eta.trans <- as.integer(eta.trans);
        setup$object <- object;
        setup$eta <- eta;
        setup$eta.mat <- matrix(eta, ncol=1);
        setup$neta <- as.integer(length(eta))
        setup$theta <- theta;
        setup$ntheta <- as.integer(length(theta))
        setup$id <- as.integer(id);
        setup$nearPD <- rxNearPd;
        if (!is.null(inits.vec)){
            setup$inits.vec <- inits.vec;
        }
        return(list2env(setup));
    }))
}

rxNearPd <- function(mat, env){
    if (any(is.na(mat))){
        env$reset <- 1;
        return(mat)
    } else {
        return(as.matrix(Matrix::nearPD(mat)$mat));
    }
}

##' Finalize Likelihood for individual.
##'
##' @param env Environment where likelihood is finalized.
##' @return The likelihood with the fitted values, posthoc eta, and possibly the individual contribution to the gradient.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiFinalizeLlik <- function(env){
    return()
}
##' FOCEI ETA setup,
##'
##' This is basically for testing
##'
##' @param object RxODE object
##' @param ... values sent to rxFoceiEtaSetup
##' @param dv dependent variable
##' @param eta eta values.
##' @param env Environment to use, instead of rxFoceiEtaSetup environment
##' @param nonmem Match NONMEMs approximation of the hessian.
##' @return environment of solved information.
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxFoceiEta <- function(object, ..., dv, eta,  env, nonmem=FALSE){
    UseMethod("rxFoceiEta");
}

##' @rdname rxFoceiEta
##' @export
rxFoceiEta.rxFocei <- function(object, ..., dv, eta, env, nonmem=FALSE){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$inner$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiEta.rxDll", "RxODE"), args, envir=parent.frame(1)));
}

##' @rdname rxFoceiEta
##' @export
rxFoceiEta.RxODE <- function(object, ..., dv, eta,  env, nonmem){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiEta.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiEta
##' @export
rxFoceiEta.rxDll <- function(object, ..., dv, eta, env, nonmem){
    if (missing(env)){
        args <- as.list(match.call(expand.dots=TRUE))[-1];
        env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    }
    object$.call(rxTrans(object)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    rxInner(eta, env);
    return(env);
}

##' Get -LLIK for individual
##'
##' Used for testing.
##'
##' @inheritParams rxFoceiEta
##' @return -llik
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiLik <- function(object, ..., dv, eta){
    UseMethod("rxFoceiLik");
}
##' @rdname rxFoceiLik
##' @export
rxFoceiLik.rxFocei <- function(object, ..., dv, eta){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$inner$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiLik.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiLik
##' @export
rxFoceiLik.RxODE <- function(object, ..., dv, eta){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiLik.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiLik
##' @export
rxFoceiLik.rxDll <- function(object, ..., dv, eta){
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
##' @export
rxFoceiLp <- function(object, ..., dv, eta){
    UseMethod("rxFoceiLp");
}
##' @rdname rxFoceiLp
##' @export
rxFoceiLp.rxFocei <- function(object, ..., dv, eta){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$inner$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiLp.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiLp
##' @export
rxFoceiLp.RxODE <- function(object, ..., dv, eta){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiLp.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiLp
##' @export
rxFoceiLp.rxDll <- function(object, ..., dv, eta){
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
rxFoceiInner <- function(object, ..., dv, eta,
                         invisible = 0, m = 6, epsilon = 1e-05, past = 0, delta = 0,
                         max_iterations = 0, linesearch_algorithm = "LBFGS_LINESEARCH_DEFAULT",
                         max_linesearch = 20, min_step = 1e-20, max_step = 1e+20,
                         ftol = 1e-04, wolfe = 0.9, gtol = 0.9, orthantwise_c = 0,
                         orthantwise_start = 0, orthantwise_end = length(eta)){
    inner.dll <- object$inner$cmpMgr$rxDll();
    inner.dll$.call(rxTrans(inner.dll)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- inner.dll;
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    lik <- RxODE_focei_eta("lik");
    lp <- RxODE_focei_eta("lp")
    output <- lbfgs::lbfgs(lik, lp, eta, environment=env,
                           invisible=invisible, m=m, epsilon=epsilon, past=past, delta=delta,
                           max_iterations=max_iterations, linesearch_algorithm=linesearch_algorithm,
                           max_linesearch=max_linesearch, min_step=min_step, max_step=max_step,
                           ftol=ftol, wolfe=wolfe, gtol=gtol, orthantwise_c=orthantwise_c,
                           orthantwise_start=orthantwise_start, orthantwise_end = orthantwise_end);
    if (any(is.na(env$eta))){
        warning("ETA estimate failed; keeping prior");
        env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
    }
    if (any(env$eta > 1e4)){
        warning("ETA estimate overflow; keeping prior");
        env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
    }
    if (is.null(object$outer)){
        return(tryCatch({RxODE_focei_finalize_llik(env)},
                        error=function(e){
            if (any(is.na(as.vector(env$H)))){
                cat(sprintf("Warning: Underflow/overflow; resetting ETAs to 0 (ID=%s).\n", env$id));
                args$eta <- rep(0, length(env$eta));
                env$eta <- args$eta;
                output <- lbfgs::lbfgs(lik, lp, eta, environment=env,
                                       invisible=invisible, m=m, epsilon=epsilon, past=past, delta=delta,
                                       max_iterations=max_iterations, linesearch_algorithm=linesearch_algorithm,
                                       max_linesearch=max_linesearch, min_step=min_step, max_step=max_step,
                                       ftol=ftol, wolfe=wolfe, gtol=gtol, orthantwise_c=orthantwise_c,
                                       orthantwise_start=orthantwise_start, orthantwise_end = orthantwise_end)
                if (any(is.na(env$eta))){
                    warning("ETA estimate failed; Assume ETA=0");
                    env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
                }
                if (any(env$eta > 1e4)){
                    warning("ETA estimate overflow; Assume ETA");
                    env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
                }
                return(RxODE_focei_finalize_llik(env));
            } else {
                cat("Warning: Hessian not positive definite (correcting with nearPD)\n");
                env$H <- -rxNearPd(-env$H, env);
                return(RxODE_focei_finalize_llik(env))
            }
        }))
    } else {
        args$object <- object;
        args$eta <- env$eta;
        ## print(env$id)
        ## print(args$eta);
        ## print(eta);
        ## print(output);
        env <- do.call(getFromNamespace("rxFoceiTheta", "RxODE"), args, envir = parent.frame(1));
        if (env$reset == 1){
            if (reset){
                stop("Cannot correct ETA underflow/overflow.");
            } else {
                cat(sprintf("Warning: Underflow/overflow; resetting ETAs to 0 (ID=%s).\n", env$id));
                print(env$theta);
                args$eta <- rep(0, length(env$eta));
                args$reset  <- TRUE;
                force(args)
                return(do.call(getFromNamespace("rxFoceiInner", "RxODE"), args, envir = parent.frame(1)));
            }
        }
        return(env$ret);
    }
}


##' FOCEI THETA setup,
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
##' @export
rxFoceiTheta <- function(object, ..., dv, eta, env, nonmem=FALSE){
    UseMethod("rxFoceiTheta");
}

##' @rdname rxFoceiTheta
##' @export
rxFoceiTheta.rxFocei <- function(object, ..., dv, eta, env, nonmem=FALSE){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$outer$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiTheta.rxDll", "RxODE"), args, envir=parent.frame(1)));
}

##' @rdname rxFoceiTheta
##' @export
rxFoceiTheta.RxODE <- function(object, ..., dv, eta, env, nonmem){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- object$cmpMgr$rxDll();
    return(do.call(getFromNamespace("rxFoceiTheta.rxDll", "RxODE"), args, envir=parent.frame(1)));
}
##' @rdname rxFoceiTheta
##' @export
rxFoceiTheta.rxDll <- function(object, ..., dv, eta, env, nonmem){
    object$.call(rxTrans(object)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    if (missing(env)){
        args <- as.list(match.call(expand.dots=TRUE))[-1];
        env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    }
    rxOuter(env);
    return(env);
}
