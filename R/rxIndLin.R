.rxIndLinStrategy <- "curState";
##' This sets the inductive linearization strategy for matrix building
##'
##' When there is more than one state in a ODE that cannot be
##' separated this specifies how it is incorporated into the matrix
##' exponential.
##'
##' @param strategy The strategy for inductive linearization matrix building
##'
##' \itemize{
##'
##' \item[curState] Prefer parameterizing in terms of the current
##'    state, followed by the first state observed in the term.
##'
##' \item[split] Split the paramterization between all states in the
##' term by dividing each by the number of states in the term and then
##' adding a matrix term for each state.
##' }
##'
##' @return Nothing
##' @author Matthew L. Fidler
##' @export
rxIndLinStrategy <- function(strategy=c("curState", "split")) {
    assignInMyNamespace(".rxIndLinStrategy", match.arg(strategy));
}

.rxIndLinState <- NULL
##' Set the preferred factoring by state
##'
##' @param prefered A list of each state's prefered factorization
##' @return Nothing
##' @author Matthew Fidler
##' @export
rxIndLinState <- function(prefered=NULL) {
    if (is.null(prefered)) {
        return(assignInMyNamespace(".rxIndLinState", prefered))
    }
    checkmate::assertList(prefered, names="unique")
    lapply(seq_along(prefered), function(x) {
        if (!checkmate::checkCharacter(prefered[[x]],
                                       names="unnamed")) {
            stop(sprintf("'rxIndLinState' list element '%s' must be a unnamed character vector", names(prefered)[x]))
        }
    })
    assignInMyNamespace(".rxIndLinState", prefered)
}

.rxIndLinLine <- function(line, states, state0){
    .tmp <- symengine::expand(line) ## Expand line
    .tmp <- rxFromSE(.tmp); ## Convert SE->RxODE; Changes things like X^1 -> X

    .ret <- eval(parse(text=paste0("rxSplitPlusQ(quote(", .tmp, "))")));
    .lst <- list()
    .fullIndLin <- FALSE
    ## This adds the x to the name environment in the .lst list above.
    .addIt <- function(x, name){
        .cur <- .lst
        if (any(names(.lst) == name)){
            .cur[[name]] <- c(.cur[[name]], x)
        } else {
            .cur[[name]] <- x;
        }
        .lst <<- .cur;
    }
    ## This fixes the multiplication so that -1*x becomes -x and *1* becomes *
    ## Also changes y^2*z*1/y to y^1 to y
    ## This collapses a character vector with "*" between it
    ## If we have *1/x this becomes simply /x
    .multCollapse <- function(x){
        as.character(symengine::S(paste(x, collapse="*")))
    }
    sapply(.ret, function(x){
        .mult <- eval(parse(text=paste0("rxSplitPlusQ(quote(", x, "),mult=TRUE)")))
        .mult <- unlist(lapply(.mult, function(x){
            if (x == "-1") return(x);
            .min <- substr(x, 0, 1);
            if (.min == "-"){
                c("-1", substr(x, 2, nchar(x)))
            } else {
                return(x);
            }
        }))

        .curStates <- unlist(lapply(.mult, function(x){
            .pars <- rxModelVars(paste0("rx_expr=", x))$params
            return(.pars[.pars %in% states]);
       }))
        .addState <- function(.state, .mult){
            .num <- which(.mult == .state)
            if (length(.num) == 0){
                .fullIndLin <<- TRUE
                .rest <- c(.mult, paste0("1/", .state))
                .expr <- .multCollapse(.rest)
                .addIt(.expr, .state);
            } else {
                if (length(.num) > 1) .fullIndLin <<- TRUE
                .num <- .num[1];
                .rest <- .mult[-.num];
                if (length(.rest) == 0){
                    .addIt("1", .state);
                } else if (length(.rest) == 1 && .rest[1] == "-1"){
                    .addIt("-1", .state);
                } else {
                    if (!.fullIndLin){
                        ## Check to see if this is inductive or
                        ## depends on other states.
                        .otherStates <- unlist(lapply(.rest, function(x){
                            .pars <- rxModelVars(paste0("rx_expr=", x))$params
                            return(.pars[.pars %in% states]);
                        }))
                        if (length(.otherStates) > 0) .fullIndLin <<- TRUE;
                    }
                    .expr <- .multCollapse(.rest)
                    .addIt(.expr, .state);
                }
            }
        }
        if (length(.curStates) == 1){
            .addState(.curStates, .mult)
        } else if (length(.curStates) > 1) {
            if (.rxIndLinStrategy == "split") {
                ## Use strategy #3, split between all the compartments
                .extra <- paste0("1/", length(.curStates))
                for (.s in .curStates){
                    .addState(.s, c(.mult, .extra));
                }
            } else {
                .pref <- .rxIndLinState[[state0]]
                .addPref <- FALSE
                if (!is.null(.pref)) {
                    .pref <- intersect(.pref, .curStates)
                    if (length(.pref) > 0){
                        .addState(.pref[1], .mult)
                        .addPref <- TRUE
                    }
                }
                if (!.addPref) {
                    if (any(.curStates == state0)) {
                        ## If there is d/dt(state1) = ... (state1*state2*state3) ...
                        ## or some other complex expression prefer expressing
                        .addState(state0, .mult);
                    } else {
                        ## Otherwise just use the first state identified.
                        ## Strategy #1 just add the first
                        .addState(.curStates[1], .mult);
                    }
                }
            }
        } else {
            ## This is some "constant" or time-base experssion that
            ## does not depend on state.  Hence it is the extra constant
            .expr <- .multCollapse(.mult);
            .addIt(.expr, "_rxF")
        }
    })
    .ret <- sapply(c(states, "_rxF"), function(x){
        if (any(names(.lst) == x)) return(gsub("[+][-]", "-", paste(.lst[[x]], collapse="+")))
        return("0");
    });
    rxTick();
    return(c(.ret, paste0(.fullIndLin)))
}

##' This creates the inductive linearization pieces to integrate into
##' a RxODE model.
##'
##' @param model RxODE model type of object
##'
##' @param doConst Replace constants with values; By default this is
##'     \code{FALSE}.
##'
##' @return List:
##' \itemize{
##'
##' \item Matrix Exponential initial matrix A for exp(t*A)
##'
##' \item Inductive Linerization vector for F
##'
##' \item Extra RxODE code for model generation
##'
##' \item Generated C code for model variables; With ME only this will
##' be a list of size 1, otherwise it is a list of size 2.
##'
##' }
##' @author Matthew Fidler
##' @noRd
.rxIndLin <- function(model, doConst=FALSE){
    .env <- .rxLoadPrune(model, doConst=doConst)
    .states <- rxState(.env);
    rxProgress(length(.states));
    message("create inductive linearization matrices")
    on.exit({rxProgressAbort()});
    .ret <- eval(parse(text=rxIndLin_(.states)))
    .w <- setNames(which(.ret[, "indLin"] == "TRUE") - 1, NULL);
    .fullIndLin <- length(.w) > 0
    .ret0 <- .ret[.states, .states, drop = FALSE]
    .ret1 <- .ret[, "_rxF", drop = FALSE]
    .code <- c(paste0("_rxM=", as.vector(.ret0), ";"));
    .codeSave <- c("SEXP matLst = PROTECT(allocVector(VECSXP, 4));pro++;",
                   "SEXP matNme = PROTECT(allocVector(STRSXP, 4));pro++;",
                   "SET_STRING_ELT(matNme, 0, mkChar(\"A\"));",
                   "SET_STRING_ELT(matNme, 1, mkChar(\"f\"));",
                   "SET_STRING_ELT(matNme, 2, mkChar(\"fullIndLin\"));",
                   "SET_STRING_ELT(matNme, 3, mkChar(\"wIndLin\"));",
                   "setAttrib(matLst, R_NamesSymbol, matNme);",
                   ## mat0 setup the nstate x nstate matrix for
                   ## inductive linearization
                   sprintf("SEXP mat0 = PROTECT(allocMatrix(STRSXP, %s, %s));pro++;", dim(.ret0)[1], dim(.ret0)[2]),
                   paste(paste0("SET_STRING_ELT(mat0, ", seq_along(as.vector(.ret0)) - 1, ",mkChar(\"", as.vector(.ret0), "\"));"), collapse="\n"),
                   sprintf("SEXP mat0nn = PROTECT(allocVector(STRSXP, %s));pro++;", dim(.ret0)[[1]]),
                   paste(paste0("SET_STRING_ELT(mat0nn, ", seq_along(dimnames(.ret0)[[1]]) - 1,
                                ",mkChar(\"", dimnames(.ret0)[[1]], "\"));"), collapse="\n"),
                   ## Assumes state is the SEXP that we can use for
                   ## state names.
                   "SEXP mat0n = PROTECT(allocVector(VECSXP, 2));pro++;",
                   "SET_VECTOR_ELT(mat0n, 0, mat0nn);",
                   "SET_VECTOR_ELT(mat0n, 1, mat0nn);",
                   "setAttrib(mat0, R_DimNamesSymbol, mat0n);");
    if (!all(.ret1 == "0")){
        .code <- c(.code, paste("_rxF=", as.vector(.ret1)))
        .codeSave <- c(.codeSave,
                       ## Now setup matF or the F matrix from inductive linearization
                       sprintf("SEXP matF = PROTECT(allocMatrix(STRSXP, %s, %s));pro++;", dim(.ret1)[1], dim(.ret1)[2]),
                       paste(paste0("SET_STRING_ELT(matF, ", seq_along(as.vector(.ret1)) - 1, ",mkChar(\"", as.vector(.ret1), "\"));"), collapse="\n"),
                       "SEXP matFd = PROTECT(allocVector(INTSXP, 2));pro++;",
                       ## Assumes state is the SEXP that we can use for state names.
                       "SEXP matFn = PROTECT(allocVector(VECSXP, 2));pro++;",
                       "SET_VECTOR_ELT(matFn, 0, state);",
                       "SEXP matrxF = PROTECT(allocVector(STRSXP,1)); pro++;",
                       "SET_STRING_ELT(matrxF, 0, mkChar(\"_rxF\"));",
                       "SET_VECTOR_ELT(matFn, 1, matrxF);",
                       "setAttrib(matF,R_DimNamesSymbol, matFn);",
                       "SET_VECTOR_ELT(matLst, 0, mat0);",
                       "SET_VECTOR_ELT(matLst, 1, matF);");
    } else {
        .codeSave <- c(.codeSave,
                       ## Now setup matF or the F matrix from inductive linearization
                       "SET_VECTOR_ELT(matLst, 0, mat0);",
                       "SET_VECTOR_ELT(matLst, 1, R_NilValue);"
                       );
    }
    .codeSave <- c(.codeSave,
                   "SEXP fullIndLin = PROTECT(allocVector(LGLSXP, 1)); pro++;",
                   paste0("LOGICAL(fullIndLin)[0] = ", ifelse(.fullIndLin, "1", "0"), ";"),
                   "SET_VECTOR_ELT(matLst, 2, fullIndLin);",
                   sprintf("SEXP intLin = PROTECT(allocVector(INTSXP, %s));pro++;", length(.w)),
                   paste(sprintf("INTEGER(intLin)[%s]=%s;", seq_along(.w) - 1, .w), collapse="\n"),
                   "SET_VECTOR_ELT(matLst, 3, intLin);",
                   "SET_VECTOR_ELT(lst,  17, matLst);"
                   );
    ## Generate extra RxODE code for generating right functions.
    .code <- paste(.code, collapse="\n");
    ## Generate C code for .ret0 and .ret1
    .ret <- list(.ret0, .ret1, .code,
                 paste(.codeSave, collapse="\n"));
    rxProgressStop();
    return(.ret);
}
