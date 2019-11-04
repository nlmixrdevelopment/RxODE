## Refactor Jacobian calculation
##' Faster expand.grid
##'
##' Only support x and y as characters right now
##'
##' @param x first element (must be character)
##' @param y second element (must be character)
##' @param type Internal type=0L is traditional expand gridn and
##'     type=1L is jacobian expand grid (adds symbols)
##' @return Expand grid
##' @author Matthew Fidler
##' @keywords internal
##' @examples
##'
##' ##
##' rxExpandGrid(letters, letters)
##'
##' ## Another fast method; See
##' ## https://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid
##'
##' expand.grid.jc <- function(seq1,seq2) {
##'  cbind(Var1 = rep.int(seq1, length(seq2)),
##'   Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
##' }
##'
##' \dontrun{
##'  microbenchmark::microbenchmark(rxExpandGrid(letters, letters), expand.grid.jc(letters, letters))
##' }
##' @export
rxExpandGrid <- function(x, y, type=0L){
    rxExpandGrid_(x, y, type)
}

## Assumes model is loaded.
.rxJacobian <- function(model, vars=TRUE){
    if (rxIs(vars,"logical")){
        if (vars){
            .pars  <- .rxParams(model, TRUE);
            if (any(.pars=="ETA[1]")){
                .pars  <- .pars[regexpr(rex::rex(start,"ETA[",any_numbers,"]"), .pars) != -1]
            }
            .jac <- rxExpandGrid(rxState(model),
                                 c(rxState(model), .pars),
                                 1L);
        } else  {
            .pars <- c();
            .jac <- rxExpandGrid(rxState(model),
                                 rxState(model),
                                 1L);
        }
    } else if (rxIs(vars,"character")){
        .pars <- vars;
        .jac <- rxExpandGrid(rxState(model),
                             c(rxState(model), vars),
                             1L)
    }
    assign("..vars", .pars, envir=model);
    rxCat("Calculate Jacobian\n");
    rxProgress(dim(.jac)[1]);
    on.exit({rxProgressAbort()});
    .ret <- apply(.jac, 1, function(x){
        .l <- x["line"];
        .l <- eval(parse(text=.l));
        rxTick();
        paste0(x["rx"], "=", rxFromSE(.l));
    })
    assign("..jacobian", .ret, envir=model);
    rxProgressStop();
    return(.ret)
}

## Assumes .rxJacobian called on model c(state,vars)
.rxSens <- function(model, vars, vars2){
    .state <- rxState(model);
    if (length(.state) > 0L){
        if (missing(vars)) vars <- get("..vars", envir=model);
        if (!missing(vars2)){
            .grd <- rxExpandSens2_(.state, vars, vars2);
        } else {
            .grd <- rxExpandSens_(.state, vars);
        }
        message("calculate sensitivities");
        rxProgress(dim(.grd)[1]);
        on.exit({rxProgressAbort()});
        lapply(c(.grd$ddtS, .grd$ddS2), function(x){
            assign(x, symengine::Symbol(x), envir=model)
        })
        .ret <- apply(.grd, 1, function(x){
            .l <- x["line"];
            .l <- eval(parse(text=.l));
            .ret <- paste0(x["ddt"], "=", rxFromSE(.l));
            if (exists(x["s0"], envir=model)){
                .l <- x["s0D"];
                .l <- eval(parse(text=.l));
                if (paste(.l) != "0"){
                    .ret <- paste0(.ret, "\n", x["s0r"], "=", rxFromSE(.l),
                                   "+0.0");
                }
            }
            rxTick();
            return(.ret);
        })
        if (missing(vars2)){
            assign("..sens", .ret, envir=model);
        } else {
            assign("..sens2", .ret, envir=model);
        }
        rxProgressStop();
    } else {
        assign("..sens", NULL, envir=model)
        .ret <- NULL
    }
    return(.ret)
}


## Check for good functions for predfn and pkpars and error functions
.goodFns <- c(".GlobalEnv", "package:RxODE", "package:nlmixr")
.checkGood <- function(x){
    .tmp <- suppressWarnings({find(deparse(substitute(x)))})
    if (!identical(.tmp, character())){
        if (!any(.tmp == .goodFns)){
            stop(sprintf("'%s' is from '%s' and cannot be used in this context", deparse(substitute(x)), .tmp))
        }
    }
}

##' Get the state information for the current model
##'
##' @param obj RxODE model that can take rxModelVars
##' @return character vector of initial preserved states (state),
##'     extra states (statef) and dvid translation information
##'     (dvid). This is used in generating the final RxODE model.
##' @author Matthew Fidler
.rxGenFunState <- function(obj){
    .mv0 <- rxModelVars(obj);
    .curDvid <- .mv0$dvid;
    .state0  <- .mv0$state;
    if (length(.state0) > 0){
        .state0 <- paste(paste0("cmt(",.mv0$state,");\n"),collapse="");
    } else {
        .state0 <- "";
    }
    .statef  <- .mv0$stateExtra;
    if (length(.statef) > 0){
        .statef <- paste0(paste(paste0("\ncmt(",.mv0$stateExtra,");"),collapse=""),"\n");
    } else {
        .statef <- "";
    }
    if (length(.curDvid) > 1){
        .dvidF <- paste0("\ndvid(",
                         paste(.curDvid, collapse=","),");\n");
    } else {
        .dvidF <- "";
    }
    return(c(state=.state0, statef=.statef, dvid=.dvidF))
}

##' Generate a RxODE block augmented by a pkpars function
##'
##' @inheritParams rxSEinner
##' @return New model augmented by pred function
##' @author Matthew L. Fidler
.rxGenPkpars <- function(rx, pkpars, init){
    if (!is.null(pkpars)){
        .txt <- as.vector(unlist(strsplit(rxParsePk(pkpars, init=init), "\n")));
        .re <- rex::rex(start, "init",
                        or("_", ".", ""),
                        or("C", "c"), "ond",
                        or("ition", ""),
                        or("", "s"), any_spaces,
                        or("=", "<-", "~"), any_spaces,
                        "c(", capture(anything), ")", any_spaces, ";",
                        any_spaces, end);
        .w <- which(regexpr(.re, .txt) != -1);
        if (length(.w) == 1L){
            .inis <- txt[.w];
            .inis <- strsplit(gsub(.re, "\\1", .inis), " *, *")[[1]];
            if (length(rxState(rx)) == length(inis)){
                .inis <- paste(paste0(rxState(rx), "(0)=", inis, "+0.0;"), collapse="\n");
                .txt[w] <- .inis;
            } else {
                stop("specified %s initial conditions when there are only %s states", length(.inis), length(rxState(rx)));
            }
        } else if (length(.w) > 1){
            stop("only one 'initCondition=' supported");
        }
        .newmod <- rxGetModel(paste0(paste(.txt, collapse="\n"), "\n", rxNorm(rx)));
        return(.newmod);
    } else {
        return(rx)
    }
}
##' Generate Augmented Pred/err RxODE model
##'
##' @inheritParams rxSEinner
##' @return A list of (1) RxODE model variables augmented with
##'     pred/error information and (2) extra error variables created.
##' @author Matthew L Fidler
.rxGenPred <- function(obj, predfn, errfn, init){
    .extraPars <- c();
    if (is.null(errfn)){
        errfn <- function(){add(0.1)};
    }
    ## Get maximum theta.
    .txt <- rxParsePred(predfn, init=init);
    .predMod <- rxGetModel(.txt);
    .pars <- .rxParams(rxGetModel(paste0(rxNorm(obj), "\n",
                                        rxNorm(.predMod))), FALSE);
    .w <- which(sapply(.pars, function(x){
        return(substr(x, 0, 2) == "TH")
    }))
    .pars <- .pars[.w];
    .mtheta <- max(as.numeric(gsub(rex::rex("THETA[", capture(numbers), "]"), "\\1", .pars)))
    .err <- rxParseErr(errfn, base.theta=.mtheta + 1, init=init)
    .predMod <- rxParsePred(predfn, init=init, .err);
    .extraPars <- attr(.predMod, "ini");
    .predMod <- rxGetModel(.predMod)
    .full <- paste0(rxNorm(obj), "\n", rxNorm(.predMod));
    .full <- rxGetModel(.full);
    return(list(.full, .extraPars));
}

.rxLoadPrune <- function(mod, doConst=TRUE, promoteLinSens=TRUE, fullModel=FALSE){
    if (fullModel){
        message(sprintf("pruning branches of full model...",
            appendLF=FALSE));
    } else {
        message(sprintf("pruning branches...",
            appendLF=FALSE))
    }
    .newmod <-rxGetModel(rxPrune(mod));
    message("done")
    ## message("Loading into symengine environment...", appendLF=FALSE)
    if (fullModel){
        message(sprintf("loading full model into symengine environment"));
    } else {
        message(sprintf("loading into symengine environment"))
    }
    .newmod <- rxS(.newmod, doConst, promoteLinSens=promoteLinSens);
    ## message("done.")
    return(.newmod)
}

##' Generate and load RxODE function into sympy environment
##'
##' @inheritParams rxSEinner
##' @return Sympy environment
##' @author Matthew Fidler
.rxGenFun <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                      init=NULL,
                      promoteLinSens=TRUE){
    rxSolveFree();
    rxTempDir();
    .checkGood(predfn);
    .checkGood(pkpars);
    .checkGood(errfn);
    ## Probably need assignInMyNamespace...
    message("creating full model...", appendLF=FALSE)
    .stateInfo <- .rxGenFunState(obj);
    .newmod <- .rxGenPkpars(obj, pkpars, init);
    .newmod <- .rxGenPred(.newmod, predfn, errfn, init);
    .extraPars <- .newmod[[2]]
    .newmod <- .newmod[[1]]
    .newmod <- .rxLoadPrune(.newmod, promoteLinSens=promoteLinSens, )
    message("done")
    .newmod$..stateInfo <- .stateInfo
    .newmod$..extraPars <- .extraPars
    return(.newmod)
}
##' Generate the ETA sensitivities for FO related methods
##'
##' @inheritParams rxSEinner
##' @return RxODE/symengine environment
##' @author Matthew L. Fidler
.rxGenEtaS <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                       init=NULL, promoteLinSens=TRUE,
                       theta=FALSE){
    .s <- .rxGenFun(obj, predfn, pkpars, errfn, init,
                    promoteLinSens=promoteLinSens)
    .etaVars <- c()
    if (theta && exists("..maxTheta", .s)){
        .etaVars <- paste0("THETA_", seq(1, .s$..maxTheta), "_")
    } else if (exists("..maxEta", .s)) {
        .etaVars <- paste0("ETA_", seq(1, .s$..maxEta), "_")
    }
    if (length(.etaVars) == 0L){
        stop("cannot identify parameters for sensitivity analysis");
    }
    .stateVars <- rxState(.s)
    .s <- .rxGenFun(obj, predfn, pkpars, errfn, init,
                    promoteLinSens=promoteLinSens)
    .rxJacobian(.s, c(.stateVars, .etaVars))
    .rxSens(.s, .etaVars)
    return(.s)
}

##' Generate the d(err)/d(eta) values for FO related methods
##'
##' @inheritParams rxSEinner
##' @return RxODE/symengine environment
##' @author Matthew L. Fidler
.rxGenHdEta <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                        init=NULL, pred.minus.dv=TRUE,
                        promoteLinSens=TRUE,
                        theta=FALSE){
    ## Equation 19 in Almquist
    .s <- .rxGenEtaS(obj, predfn, pkpars, errfn, init,
                     promoteLinSens=promoteLinSens,
                     theta=theta)
    .stateVars <- rxState(.s)
    .grd <- rxExpandFEta_(.stateVars, .s$..maxEta,
                          ifelse(pred.minus.dv, 1L, 2L))
    if (.useUtf()){
        rxCat("Calculate \u2202(f)/\u2202(\u03B7)\n");
    } else {
        rxCat("Calculate d(f)/d(eta)\n");
    }
    rxProgress(dim(.grd)[1]);
    on.exit({rxProgressAbort()});
    .any.zero <- FALSE
    .all.zero <- TRUE
    .ret <- apply(.grd, 1, function(x){
        .l <- x["calc"];
        .l <- eval(parse(text=.l));
        .ret <- paste0(x["dfe"], "=", rxFromSE(.l));
        .zErr <- suppressWarnings(try(as.numeric(get(x["dfe"], .s)),silent=TRUE))
        if (identical(.zErr, 0)){
            .any.zero <<- TRUE
        } else if (.all.zero){
            .all.zero <<- FALSE
        }
        rxTick();
        return(.ret);
    })
    if (.all.zero){
        stop("none of the predictions depend on 'ETA'")
    }
    if (.any.zero){
        warning("some of the predictions do not depend on 'ETA'")
    }
    .s$..HdEta <- .ret;
    .s$..pred.minus.dv <- pred.minus.dv;
    rxProgressStop();
    return(.s)
}
##' Finalize RxODE pred based on symengine saved info
##'
##' @param .s Symengine/RxODE object
##' @inheritParams rxSEinner
##' @return Nothing
##' @author Matthew L Fidler
.rxFinalizePred <- function(.s, sum.prod=FALSE,
                             optExpression=TRUE){
    .prd <- get("rx_pred_", envir=.s);
    .prd <- paste0("rx_pred_=", rxFromSE(.prd))
    .r <- get("rx_r_", envir=.s);
    .r <- paste0("rx_r_=", rxFromSE(.r))
    .yj <- paste(get("rx_yj_", envir=.s));
    .yj <- paste0("rx_yj_~", rxFromSE(.yj))
    .lambda <- paste(get("rx_lambda_", envir=.s));
    .lambda <- paste0("rx_lambda_~", rxFromSE(.lambda))
    .lhs0 <- .s$..lhs0
    if (is.null(.lhs0)) .lhs0 <- ""
    .lhs <- .s$..lhs
    if (is.null(.lhs)) .lhs <- ""
    .ddt <- .s$..ddt
    if (is.null(.ddt)) .ddt <- ""
    .s$..pred <- paste(c(.s$..stateInfo["state"],
                         .lhs0,
                         .ddt,
                         .yj,
                         .lambda,
                         .prd,
                         .r,
                         .lhs,
                         .s$..stateInfo["statef"],
                         .s$..stateInfo["dvid"],
                         ""), collapse="\n")
    .s$..pred.nolhs <- paste(c(.s$..stateInfo["state"],
                               .lhs0,
                               .ddt,
                               .yj,
                               .lambda,
                               .prd,
                               .r,
                               .s$..stateInfo["statef"],
                               .s$..stateInfo["dvid"],
                               ""), collapse="\n")
    if (sum.prod){
        message("stabilizing round off errors in predictions or EBE model...", appendLF=FALSE);
        .s$..pred <- rxSumProdModel(.s$..pred);
        message("done");
    }
    if (optExpression){
        .s$..pred <- rxOptExpr(.s$..pred, "EBE model")
    }
}
##' Finalize inner RxODE based on symengine saved info
##'
##' @param .s Symengine/RxODE object
##' @inheritParams rxSEinner
##' @return Nothing
##' @author Matthew L Fidler
.rxFinalizeInner <- function(.s, sum.prod=FALSE,
                             optExpression=TRUE){
    .prd <- get("rx_pred_", envir=.s);
    .prd <- paste0("rx_pred_=", rxFromSE(.prd))
    .r <- get("rx_r_", envir=.s);
    .r <- paste0("rx_r_=", rxFromSE(.r))
    .yj <- paste(get("rx_yj_", envir=.s));
    .yj <- paste0("rx_yj_~", rxFromSE(.yj))
    .lambda <- paste(get("rx_lambda_", envir=.s));
    .lambda <- paste0("rx_lambda_~", rxFromSE(.lambda))
    .ddt <- .s$..ddt
    if (is.null(.ddt)) .ddt <- character(0)
    .sens <- .s$..sens
    if (is.null(.sens)) .sens <- character(0)
    .s$..inner <- paste(c(.ddt,
                          .sens,
                          .yj,
                          .lambda,
                          .prd,
                          .s$..HdEta,
                          .r,
                          .s$..REta,
                          .s$..stateInfo["statef"],
                          .s$..stateInfo["dvid"],
                          ""), collapse="\n")
    if (sum.prod){
        message("stabilizing round off errors in inner problem...", appendLF=FALSE);
        .s$..inner <- rxSumProdModel(.s$..inner);
        message("done");
    }
    if (optExpression){
        .s$..inner <- rxOptExpr(.s$..inner, "inner model")
    }
    .s$..inner <- paste0(.s$..stateInfo["state"],"\n", .s$..inner);
}
##' Generate FOCE without interaction
##'
##' @inheritParams rxSEinner
##' @return RxODE/symengine environment
##' @author Matthew Fidler
.rxGenFoce <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                       init=NULL, pred.minus.dv=TRUE,
                       sum.prod=FALSE,
                       optExpression=TRUE,
                       promoteLinSens=TRUE,
                       theta=FALSE){
    .s <- .rxGenHdEta(obj, predfn, pkpars, errfn, init, pred.minus.dv,
                      promoteLinSens=promoteLinSens, theta=theta)
    .s$..REta <- NULL;
    ## Take etas from rx_r
    eval(parse(text=rxRepR0_(.s$..maxEta)))
    .rxFinalizeInner(.s, sum.prod, optExpression)
    .rxFinalizePred(.s, sum.prod, optExpression)
    .s$..outer <- NULL
    return(.s)
}

##' Generate peices for FOCEi inner problem
##'
##' @inheritParams rxSEinner
##' @return RxODE/symengine environment
##' @author Matthew L. Fidler
.rxGenFocei <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                        init=NULL, pred.minus.dv=TRUE,
                        sum.prod=FALSE,
                        optExpression=TRUE,
                        promoteLinSens=TRUE, theta=FALSE){
    .s <- .rxGenHdEta(obj, predfn, pkpars, errfn, init, pred.minus.dv,
                      promoteLinSens=promoteLinSens, theta=theta)
    .stateVars <- rxState(.s)
    .grd <- rxExpandFEta_(.stateVars, .s$..maxEta, FALSE)
    if (.useUtf()){
        rxCat("Calculate \u2202(R\u00B2)/\u2202(\u03B7)\n");
    } else {
        rxCat("Calculate d(R^2)/d(eta)\n");
    }
    rxProgress(dim(.grd)[1]);
    on.exit({rxProgressAbort()});
    .ret <- apply(.grd, 1, function(x){
        .l <- x["calc"];
        .l <- eval(parse(text=.l));
        .ret <- paste0(x["dfe"], "=", rxFromSE(.l));
        rxTick();
        return(.ret);
    })
    .s$..REta <- .ret;
    rxProgressStop();
    .rxFinalizeInner(.s, sum.prod, optExpression)
    .rxFinalizePred(.s, sum.prod, optExpression)
    .s$..outer <- NULL
    return(.s)
}

##' Generate peices for EBE only problem
##'
##' @inheritParams rxSEinner
##' @return RxODE/symengine environment
##' @author Matthew L. Fidler
.rxGenEBE <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                      init=NULL, pred.minus.dv=TRUE,
                      sum.prod=FALSE,
                      optExpression=TRUE,
                      theta=FALSE){
    .s <- .rxGenFun(obj, predfn, pkpars, errfn, init)
    .s$..inner <- NULL
    .s$..outer <- NULL
    .rxFinalizePred(.s, sum.prod, optExpression)
    return(.s);
}

rxSymPyAbsLog <- FALSE
rxSymPyLogSign <- c()

rxSymPyExpThetas <- c()
rxSymPyExpEtas <- c()


##' Setup Pred function based on RxODE object.
##'
##' This is for the so-called inner problem.
##'
##' @param obj RxODE object
##' @param predfn Prediction function
##' @param pkpars Pk Pars function
##' @param errfn Error function
##' @param init Initialization parameters for scaling.
##' @param grad Boolaen indicated if the the equations for the
##'     gradient be calculated
##' @param sum.prod A boolean determining if RxODE should use more
##'     numerically stable sums/products.
##' @param pred.minus.dv Boolean stating if the FOCEi objective
##'     function is based on PRED-DV (like NONMEM).  Default TRUE.
##' @param only.numeric Instead of setting up the sensitivities for
##'     the inner problem, modify the RxODE to use numeric
##'     differentiation for the numeric inner problem only.
##' @param optExpression Optimize the model text for computer
##'     evaluation.
##' @param interaction Boolean to determine if dR^2/deta is calculated
##'     for FOCEi (not needed for FOCE)
##' @param theta Calculate THETA derivatives instead of ETA
##'     derivatives.  By default FALSE
##' @return RxODE object expanded with predfn and with calculated
##'     sensitivities.
##' @inheritParams rxS
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
##' @importFrom utils find
rxSEinner <- function(obj, predfn, pkpars=NULL, errfn=NULL, init=NULL,
                      grad=FALSE, sum.prod=FALSE, pred.minus.dv=TRUE,
                      only.numeric=FALSE,
                      optExpression=TRUE,
                      interaction=TRUE, ...,
                      promoteLinSens=TRUE,
                      theta=FALSE){
    assignInMyNamespace("rxErrEnv.lambda", NULL);
    assignInMyNamespace("rxErrEnv.yj", NULL);
    assignInMyNamespace("rxSymPyExpThetas", c());
    assignInMyNamespace("rxSymPyExpEtas", c());
    if (only.numeric){
        .s <- .rxGenEBE(obj, predfn, pkpars, errfn, init, pred.minus.dv,
                        sum.prod, optExpression, theta=theta)
    } else if (interaction){
        .s <- .rxGenFocei(obj, predfn, pkpars, errfn, init, pred.minus.dv,
                          sum.prod, optExpression,
                          promoteLinSens=promoteLinSens, theta=theta);
    } else {
        .s <- .rxGenFoce(obj, predfn, pkpars, errfn, init, pred.minus.dv,
                         sum.prod, optExpression,
                         promoteLinSens=promoteLinSens, theta=theta);
    }
    .toRx <- function(x, msg){
        if (is.null(x)) return(NULL)
        message(msg, appendLF=FALSE)
        .ret <- RxODE(x);
        message("done")
        return(.ret);
    }
    if (exists("..maxTheta", .s)){
        .eventTheta <- rep(0L,.s$..maxTheta)
    } else {
        .eventTheta <- integer();
    }
    if (exists("..maxEta", .s)){
        .eventEta <- rep(0L,.s$..maxEta)
    } else {
        .eventEta <- integer();
    }
    for (.v in .s$..eventVars){
        .vars <- rxGetModel(paste0("rx_lhs=", as.character(get(.v, envir=.s))))$params
        for (.v2 in .vars){
            .reg <- rex::rex(start, "ETA_", capture(any_numbers), "_", end)
            if (regexpr(.reg, .v2) != -1){
                .num <- as.numeric(sub(.reg, "\\1", .v2))
                .eventEta[.num] <- 1L;
            }
            .reg <- rex::rex(start, "THETA_", capture(any_numbers), "_", end)
            if (regexpr(.reg, .v2) != -1){
                .num <- as.numeric(sub(.reg, "\\1", .v2))
                .eventTheta[.num] <- 1L;
            }
        }
    }
    pred.opt <- NULL
    inner <- .toRx(.s$..inner, "compiling inner model...");
    if (any(.eventEta == 1L) && !is.null(inner)){
        if (sum.prod){
            message("stabilizing round off errors in events FD model...", appendLF=FALSE);
            .s$..pred.nolhs <- rxSumProdModel(.s$..pred.nolhs);
            message("done");
        }
        if (optExpression){
            .s$..pred.nolhs <- rxOptExpr(.s$..pred.nolhs, "events FD model")
        }
        .s$..pred.nolhs <- paste(c(paste0("rx_dum_", seq_along(inner$params), "~", inner$params),
                                   .s$..pred.nolhs), collapse="\n");
        pred.opt <- .s$..pred.nolhs
    }
    .ret <- list(obj=obj,
                 inner=inner,
                 pred.only=.toRx(.s$..pred, "compiling EBE model..."),
                 extra.pars=.s$..extraPars,
                 outer=.toRx(.s$..outer),
                 pred.nolhs=.toRx(pred.opt, "compiling events FD model..."),
                 theta=NULL,
                 ## warn=.zeroSens,
                 pred.minus.dv=pred.minus.dv,
                 log.thetas=.s$..extraTheta[["exp"]],
                 log.etas=.s$..extraEta[["exp"]],
                 extraProps=.s$..extraTheta,
                 eventTheta=.eventTheta,
                 eventEta=.eventEta
                 ## ,
                 ## cache.file=cache.file
                 )
    class(.ret) <- "rxFocei";
    return(.ret)
}
##'@rdname rxSEinner
##'@export
rxSymPySetupPred <- rxSEinner

## norm <- RxODE("
## d/dt(y)  = dy
## d/dt(dy) = mu*(1-y^2)*dy - y
## ## Initial conditions
## y(0) = 2
## dy(0) = 0
## ## mu
## mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
## ")