rxFoceiEtaSetup <- function(object, ..., dv, eta, theta, nonmem=FALSE, inv.env=parent.frame(1), id= -1,
                            inits.vec=NULL, atol.outer=1e-8, rtol.outer=1e-6,
                            switch.solver=FALSE, pred.minus.dv=TRUE){
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
        setup$atol.outer <- atol.outer;
        setup$rtol.outer <- rtol.outer;
        setup$nearPD <- rxNearPd;
        setup$rc <- 0L;
        setup$switch.solver <- as.integer(switch.solver);
        setup$pred.minus.dv <- as.integer(pred.minus.dv);
        if (!is.null(inits.vec)){
            setup$inits.vec <- inits.vec;
        }
        return(list2env(setup));
    }))
}

rxNearPd <- function(mat, env){
    if (any(is.na(mat))){
        cat("Bad matrix:\n");
        print(mat);
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
##' @param ... values sent to rxFoceiEtaSetup and lbfgs
##' @param dv dependent variable
##' @param eta eta value
##' @param estimate Boolean indicating if the optimization should be
##'     perfomed(TRUE), or just keep the eta, and calculate the Loglik
##'     with fitted/posthoc attributes(FALSE).
##' @return Loglik with fitted and posthoc attributes
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFoceiInner <- function(object, ..., dv, eta, eta.bak=NULL,
                         estimate=TRUE){
    inner.dll <- object$inner$cmpMgr$rxDll();
    inner.dll$.call(rxTrans(inner.dll)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    args$object <- inner.dll;
    env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
    lik <- RxODE_focei_eta("lik");
    lp <- RxODE_focei_eta("lp")
    est <- function(){
        if (estimate){
            if (class(args$eta) == "call"){
                args$eta <- eval(args$eta, envir=parent.frame(2))
                args$orthantwise_end <- length(args$eta);
            }
            args$call_eval <- lik;
            args$call_grad <- lp;
            args$vars <- args$eta;
            args$environment <- env
            output <- do.call(getFromNamespace("lbfgs","lbfgs"), args, envir = parent.frame(1))
            rxInner(output$par, env);
        } else {
            rxInner(args$eta, env);
        }
    }
    est();
    if (any(is.na(env$eta))){
        warning("ETA estimate failed; keeping prior");
        env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
    }
    if (any(abs(env$eta) > 1e4)){
        warning("ETA estimate overflow; keeping prior");
        env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
    }
    if (is.null(object$outer)){
        ret <- RxODE_focei_finalize_llik(env);
        if (attr(ret, "corrected") == 1){
            cat(sprintf("Warning: Problem with Hessian or ETA estimate, resetting ETAs to 0 (ID=%s).\n", env$id));
            args$eta <- rep(0, length(env$eta));
            env$eta <- args$eta;
            args$orthantwise_end <- length(args$eta);
            est();
            if (any(is.na(env$eta))){
                warning("ETA estimate failed; Assume ETA=0");
                env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
            }
            if (any(env$eta > 1e4)){
                warning("ETA estimate overflow; Assume ETA");
                env <- do.call(getFromNamespace("rxFoceiEta", "RxODE"), args, envir = parent.frame(1));
            }
            ret <- RxODE_focei_finalize_llik(env);
            attr(ret, "corrected") <- 1L;
            return(ret);
        }
        return(ret);
    } else {
        args$object <- object;
        args$eta <- env$eta;
        args$orthantwise_end <- length(args$eta);
        ## print(env$id)
        ## print(args$eta);
        ## print(eta);
        ## print(output);
        env <- do.call(getFromNamespace("rxFoceiTheta", "RxODE"), args, envir = parent.frame(1));
        if (env$reset == 1){
            cat(sprintf("Warning: Problem with Hessian or ETA estimate, resetting ETAs to 0 (ID=%s).\n", env$id));
            inner.dll$.call(rxTrans(inner.dll)["ode_solver_ptr"]); ## Assign the ODE pointers (and Jacobian Type)
            args$eta <- rep(0, length(env$eta));
            args$orthantwise_end <- length(args$eta);
            args$object <- inner.dll;
            env <- do.call(getFromNamespace("rxFoceiEtaSetup", "RxODE"), args, envir = parent.frame(1));
            est();
            args$object <- object;
            args$eta <- env$eta;
            args$orthantwise_end <- length(args$eta);
            env <- do.call(getFromNamespace("rxFoceiTheta", "RxODE"), args, envir = parent.frame(1));
            if (env$reset == 1){
                print(env$par)
                stop("Cannot reset problem.");
            }
            ret <- env$ret;
            attr(ret, "corrected") <- 1L;
            return(ret);
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
    env$atol.inner <- env$atol;
    env$rtol.inner <- env$rtol;
    env$atol <- env$atol.outer;
    env$rtol <- env$rtol.outer;
    on.exit({
        env$atol <- env$atol.inner;
        env$rtol <- env$rtol.inner;
    })
    rxOuter(env)
    return(env);
}
