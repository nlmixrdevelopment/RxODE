rxFoceiEtaSetup <- function(object, ..., dv, eta, omegaInv){
    args <- list(object=object, ..., eta=eta);
    setup <- do.call("rxSolveSetup", args, envir = parent.frame(1));
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
        setup$DV <- as.double(dv);
        setup$omegaInv <- omegaInv;
        setup$eta <- eta;
        setup$eta.trans <- as.integer(eta.trans);
        setup$object <- object;
        setup$eta.mat <- matrix(eta, ncol=1);
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
##' @return environment of solved information.
##' @author Matthew L. Fidler
##' @keywords internal
rxFoceiEta <- function(object, ..., dv, eta, omegaInv){
    UseMethod("rxFoceiEta");
}
##' @rdname rxFoceiEta
rxFoceiEta.RxODE <- function(object, ..., dv, eta, omegaInv){
    return(rxFoceiEta.rxDll(object=object$cmpMgr$rxDll(), ..., dv=dv, eta=eta, omegaInv=omegaInv));
}
##' @rdname rxFoceiEta
rxFoceiEta.rxDll <- function(object, ..., dv, eta, omegaInv){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call("rxFoceiEtaSetup", args, envir = parent.frame(1));
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
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call("rxFoceiEtaSetup", args, envir = parent.frame(1));
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
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call("rxFoceiEtaSetup", args, envir = parent.frame(1));
    return(RxODE_focei_eta_lp(eta, env))
}

rxFoceiInner <- function(object, ..., dv, eta, omegaInv, epsilon=1e-4, invisible=1){
    UseMethod("rxFoceiInner");
}
rxFoceiInner.RxODE <- function(object, ..., dv, eta, omegaInv, epsilon=1e-4, invisible=1){
    return(rxFoceiInner.rxDll(object=object$cmpMgr$rxDll(), ..., dv=dv, eta=eta, omegaInv=omegaInv,
                              epsilon=epsilon, invisible=invisible));
}
rxFoceiInner.rxDll <- function(object, ..., dv, eta, omegaInv, epsilon=1e-4, invisible=1,
                               slow=TRUE){
    args <- as.list(match.call(expand.dots=TRUE))[-1];
    env <- do.call("rxFoceiEtaSetup", args, envir = parent.frame(1));
    lik <- RxODE_focei_eta("lik");
    lp <- RxODE_focei_eta("lp")
    output <- lbfgs::lbfgs(lik, lp, eta, environment=env, epsilon=epsilon, invisible=invisible);
    return(output)
}
