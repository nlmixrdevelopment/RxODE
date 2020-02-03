##' Parameters specified by the model
##'
##' This returns the model's parameters that are required to solve the
##' ODE system, and can be used to pipe parameters into an RxODE solve
##'
##' @inheritParams rxModelVars
##'
##' @param constants is a boolean indicting if constants should be
##'     included in the list of parameters. Currently RxODE parses
##'     constants into variables in case you wish to change them
##'     without recompiling the RxODE model.
##'
##' @inheritParams rxControl
##'
##' @return When extracting the parameters from an RxODE model, a
##'     character vector listing the parameters in the model.
##'
##' @author Matthew L.Fidler
##' @export
rxParams <- function(obj, ...){
    UseMethod("rxParams");
}

##' @rdname rxParams
##' @export
rxParams.RxODE <- function(obj, constants=TRUE, ...,
                           params=NULL, inits=NULL, iCov=NULL,
                           keep=NULL,
                           thetaMat = NULL,
                           omega = NULL, iov = NULL, dfSub = NULL,
                           sigma=NULL, dfObs = NULL, dfOcc=NULL,
                           nSub=NULL, nStud=NULL){
    .ret <- list(params=params, inits=inits, iCov=iCov, keep=keep,
                 thetaMat = thetaMat,
                 omega = omega, iov=iov, dfSub = dfSub,
                 sigma=sigma, dfObs = dfObs, dfOcc = dfOcc,
                 nSub=nSub, nStud=nStud);
    if (all(sapply(seq_along(.ret), function(x){ is.null(.ret[[x]]) }))){
        if (length(list(...)) > 0){
            .clearPipe();
            stop("unknown arguments in 'rxParams'");
        }
        return(rxParams.default(obj, constants=constants));
    } else {
        .lst <- list(...);
        if (length(.lst) > 0){
            .clearPipe();
            stop(sprintf("unknown arguments in 'rxParams': %s\ntry piping to 'rxSolve'",paste(names(.lst), collapse=", ")));
        }
        ## Most likely
        ## RxODE() %>% rxParams() %>%
        assignInMyNamespace(".pipelineRx", obj);
        assignInMyNamespace(".pipelineInits", NULL)
        assignInMyNamespace(".pipelineEvents", NULL)
        assignInMyNamespace(".pipelineParams", NULL)
        assignInMyNamespace(".pipelineICov", NULL)
        assignInMyNamespace(".pipelineKeep", NULL)
        assignInMyNamespace(".pipelineThetaMat", NULL)
        assignInMyNamespace(".pipelineOmega", NULL)
        assignInMyNamespace(".pipelineIov", NULL)
        assignInMyNamespace(".pipelineSigma", NULL)
        assignInMyNamespace(".pipelineDfObs", NULL)
        assignInMyNamespace(".pipelineDfOcc", NULL)
        assignInMyNamespace(".pipelineDfSub", NULL)
        assignInMyNamespace(".pipelineNSub", NULL)
        assignInMyNamespace(".pipelineNStud", NULL)


        class(.ret) <- "rxParams"
        return(.ret)
    }
}

##' @rdname rxParams
##' @export
rxParams.rxSolve <- function(obj, constants=TRUE, ...,
                             params=NULL, inits=NULL, iCov=NULL,
                             keep=NULL,
                             thetaMat = NULL,
                             omega = NULL, iov=NULL, dfSub = NULL,
                             sigma=NULL, dfObs = NULL, dfOcc = NULL,
                             nSub=NULL, nStud=NULL){
    .ret <- list(params=params, inits=inits, iCov=iCov, keep=keep,
                 thetaMat = thetaMat,
                 omega = omega, iov=iov, dfSub = dfSub,
                 sigma=sigma, dfObs = dfObs, dfOcc = dfOcc,
                 nSub=nSub, nStud=nStud);
    if (all(sapply(seq_along(.ret), function(x){ is.null(.ret[[x]]) }))){
        if (length(list(...)) > 0){
            .clearPipe();
            stop("unknown arguments in 'rxParams'");
        }
        return(rxParams.default(obj, constants=constants));
    } else {
        .lst <- list(...);
        if (length(.lst) > 0){
            .clearPipe();
            stop(sprintf("unknown arguments in 'rxParams': %s\ntry piping to 'rxSolve'",paste(names(.lst), collapse=", ")));
        }
        ## Most likely
        ## solveObject %>% rxParams() %>%
        .x <- obj;
        ## Assign prior information
        ## Need to extract:
        ## 1. RxODE model
        assignInMyNamespace(".pipelineRx", .x$.args.object)
        ## Events
        assignInMyNamespace(".pipelineEvents", .x$.args.events)
        ## 2. RxODE parameters
        assignInMyNamespace(".pipelineParams", .x$.args.par0);
        assignInMyNamespace(".pipelineICov", .x$.args$iCov);
        ## 3. RxODE inits
        assignInMyNamespace(".pipelineInits", .x$.args.inits);
        ## 4. RxODE thetaMat
        assignInMyNamespace(".pipelineThetaMat", .x$.args$thetaMat);
        ## 5. RxODE omega
        assignInMyNamespace(".pipelineOmega", .x$.args$omega);
        assignInMyNamespace(".pipelineIov", .x$.args$iov);
        ## 6. RxODE sigma
        assignInMyNamespace(".pipelineSigma", .x$.args$sigma);
        ## 7. RxODE dfObs
        assignInMyNamespace(".pipelineDfObs", .x$env$.args$dfObs)
        assignInMyNamespace(".pipelineDfOcc", .x$env$.args$dfOcc)
        ## 8. RxODE dfSub
        assignInMyNamespace(".pipelineDfSub", .x$env$.args$dfSub)
        class(.ret) <- "rxParams"
        return(.ret)
    }
}

##' @rdname rxParams
##' @export
rxParams.rxEt <- function(obj, ...,
                          params=NULL, inits=NULL, iCov=NULL,
                          keep=NULL,
                          thetaMat = NULL,
                          omega = NULL, iov=NULL, dfSub = NULL,
                          sigma=NULL, dfObs = NULL, dfOcc = NULL,
                          nSub=NULL, nStud=NULL){
    # et() %>% rxParams() %>%
    assignInMyNamespace(".pipelineEvents", obj);
    .lst <- list(...);
    if (length(.lst) > 0){
        .clearPipe();
        stop(sprintf("unknown arguments in 'rxParams': %s\ntry piping to 'rxSolve'",paste(names(.lst), collapse=", ")));
    }
    .ret <- list(params=params, inits=inits, iCov=iCov, keep=keep,
                 thetaMat = thetaMat, omega = omega, iov=iov, dfSub = dfSub,
                 sigma=sigma, dfObs = dfObs, dfOcc=dfOcc,
                 nSub=nSub, nStud=nStud);
    class(.ret) <- "rxParams"
    return(.ret)
}

.rxParams <- function(obj, constants=TRUE){
    .ret <- rxParams_(obj)
    if (!constants){
        .init <- RxODE::rxInit(obj);
        .ret <- .ret[!(.ret %in% names(.init))]
    }
    return(.ret);

}
##' @export
rxParams.default <- function(obj,..., constants=TRUE){
    if (!missing(obj)){
        return(.rxParams(obj, constants));
    } else {
        .lst <- list(...);
        .nm <- c("cov", "params", "inits", "iCov", "keep",
                 "thetaMat",
                 "omega", "iov", "dfSub",
                 "sigma", "dfObs", "dfOcc", "nSub", "nStud");
        .ret <- lapply(.nm, function(x){ return(.lst[[x]]) })
        names(.ret) <- .nm
        class(.ret) <- "rxParams"
        return(.ret)
    }

}

##' @rdname rxParams
##' @export
rxParam <- rxParams
