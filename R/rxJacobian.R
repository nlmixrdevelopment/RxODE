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
    .symengine <- rxIs(model, "rxS");
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
.rxSens <- function(model, vars, vars2, msg="Calculate Sensitivites"){
    .state <- rxState(model);
    if (missing(vars)) vars <- get("..vars", envir=model);
    if (!missing(vars2)){
        .grd <- rxExpandSens2_(.state, vars, vars2);
    } else {
        .grd <- rxExpandSens_(.state, vars);
    }
    message(msg);
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
            if (.l != "0"){
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
    return(.ret)
}


## Check for good functions for predfn and pkpars and error functions
.goodFns <- c(".GlobalEnv", "package:RxODE", "package:nlmixr")
.checkGood <- function(x){
    .tmp <- suppressWarnings({find(deparse(substitute(x)))})
    if (!identical(.tmp, character())){
        if (!any(.tmp == .goodFns)){
            stop(sprintf("%s is from %s and can't be used in this context.", deparse(substitute(x)), .tmp))
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
                stop("Specified %s initial conditions when there are only %s states.", length(.inis), length(rxState(rx)));
            }
        } else if (length(.w) > 1){
            stop("Multiple initCondition= found.  Please choose one.");
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
##' @return RxODE model variables augmented with pred/error information
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
    return(.full);
}

##' Generate and load RxODE function into sympy environment
##'
##' @inheritParams rxSEinner
##' @return Sympy environment
##' @author Matthew Fidler
.rxGenFun <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                      init=NULL){
    rxSolveFree();
    rxTempDir();
    .checkGood(predfn);
    .checkGood(pkpars);
    .checkGood(errfn);
    ## Probably need assignInMyNamespace...
    message("Combining pred/pk/error...", appendLF=FALSE)
    .stateInfo <- .rxGenFunState(obj);
    .newmod <- .rxGenPkpars(obj, pkpars, init);
    .newmod <- .rxGenPred(.newmod, predfn, errfn, init);
    message("done.")
    message("Pruning branches...", appendLF=FALSE)
    .newmod <-rxGetModel(rxPrune(.newmod));
    message("done.")
    message("Loading into symengine environment...", appendLF=FALSE)
    .newmod <- rxS(.newmod);
    .newmod$..stateInfo <- .stateInfo
    message("done.")
    return(.newmod)
}
##' Generate the ETA sensitivities for FO related methods
##'
##' @inheritParams rxSEinner
##' @return RxODE/symengine environment
##' @author Matthew L. Fidler
.rxGenEtaS <- function(obj, predfn, pkpars=NULL, errfn=NULL,
                        init=NULL){
    .s <- .rxGenFun(obj, predfn, pkpars, errfn, init)
    .etaVars <- paste0("ETA_", seq(1, .s$..maxEta), "_")
    .stateVars <- rxState(.s)
    .rxJacobian(.s, c(.stateVars, .etaVars))
    .rxSens(.s, .etaVars)
    return(.s)
}


.rxGenHdEta <- function(){
    ## Equation 19
}



rxSEinner <- function(obj, predfn, pkpars=NULL, errfn=NULL, init=NULL, grad=FALSE, sum.prod=FALSE, pred.minus.dv=TRUE,
                      theta.derivs=FALSE,only.numeric=FALSE,
                      grad.internal=FALSE, theta.internal=FALSE,
                      optExpression=TRUE,
                      interaction=TRUE){

}

## norm <- RxODE("
## d/dt(y)  = dy
## d/dt(dy) = mu*(1-y^2)*dy - y
## ## Initial conditions
## y(0) = 2
## dy(0) = 0
## ## mu
## mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
## ")
