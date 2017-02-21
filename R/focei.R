rxFoceiEta <- function(object, ..., dv, eta, omegaInv){
    UseMethod("rxFoceiEta");
}
rxFoceiEta.RxODE <- function(object, ..., dv, eta, omegaInv){
    return(rxFoceiEta.rxDll(object=object$cmpMgr$rxDll(), ..., dv=dv, eta=eta, omegaInv=omegaInv));
}
rxFoceiEta.rxDll <- function(object, ..., dv, eta, omegaInv){
    args <- list(object=object, ..., eta=eta);
    setup <- do.call("rxSolveSetup", args, envir = parent.frame(1));
    return(with(setup,{
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
        env <- list2env(list(DV=as.double(dv),
                             omegaInv=omegaInv,
                             eta=eta,
                             eta.trans=as.integer(eta.trans),
                             params=params,
                             inits=inits,
                             lhs_vars=lhs_vars,
                             ## events
                             time=time,
                             evid=evid,
                             amt=amt,
                             ## Covariates
                             pcov=pcov,
                             cov=cov,
                             isLocf=isLocf,
                             ## Solver options (double)
                             atol=atol,
                             rtol=rtol,
                             hmin=hmin,
                             hmax=hmax,
                             hini=hini,
                             ## Solver options ()
                             maxordn=maxordn,
                             maxords=maxords,
                             maxsteps=maxsteps,
                             stiff=stiff,
                             transit_abs=transit_abs,
                             ## Passed to build solver object.
                             object=object,
                             extra.args=extra.args))
        object$.call(rxTrans(object)["ode_solver_focei_eta"],
                     eta,
                     env)}))
}
